#include "ising.hpp"
#include "nlohmann/json.hpp"
#include <random>
#include <unsupported/Eigen/FFT>

using json = nlohmann::json;

int main(int argc, char *argv[])
{
    double hbar = 1.0;
    int N = 5;
    double J = 1.0;
    double hx = 1.0;
    double hz = 1.0;
    double Tmax = 1.0;
    double dt = 0.1;
    int seed = 1;
    double beta = 1.0;
    int averages = 1;
    int krylovDim = 20;

    if (argc != 2)
    {
        throw std::invalid_argument("usage: ./_cyclicPerturbationFast.out <maxTime>");
    }

    std::ifstream dataFile("_params.json");
    json params = json::parse(dataFile);

    hbar = params["hbar"];
    N = params["N"];
    J = params["J"];
    hx = params["hx"];
    hz = params["hz"];
    dt = params["dt"];
    seed = params["seed"];
    beta = params["beta"];
    std::string saveFolder = params["saveFolder"];
    krylovDim = params["krylovDim"];
    int dataPoints = params["dataPoints"];

    double dbeta = dt;
    int Betaiters = static_cast<int>(std::round(beta / (2.0 *dbeta)));
    if (Betaiters <= 0) Betaiters = 1;
    dbeta = beta / (Betaiters * 2.0);

    Tmax = std::stof(argv[1]);
    double mu = 2.0 * pi / Tmax;

    if (N < 0)
    {
        throw std::invalid_argument("N must be non-negative");
    }
    if (hbar <= 0.0)
    {
        throw std::invalid_argument("hbar must be positive");
    }
    if (Tmax <= 0.0 || dt <= 0.0)
    {
        throw std::invalid_argument("T and dt must be positive");
    }

    int iters = static_cast<int>(std::round(Tmax / dt));
    if (iters <= 0) iters = 1;
    dt = Tmax / iters;

    int itersPerDataPoint = iters / dataPoints;
    if (itersPerDataPoint <= 0) itersPerDataPoint = 1;
    dataPoints = iters / itersPerDataPoint;

    int dim = 1 << N;
    std::cout << "Hilbert space dimension: " << dim << std::endl;

    VectorReal correlations(dataPoints);
    correlations.setZero();
    double WeightSum = 0;

    std::string fileName = saveFolder + "specFn_N" + std::to_string(N) + "_T" + std::to_string(Tmax) + "_B" + std::to_string(beta) + ".csv";

    double numerator = 0.0;
    double denominator = 0.0;

    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    for (int s = 0; s < dim; s++)
    {
        Vector state1(dim);
        state1.setZero();
        state1(s) = 1.0;

        Vector state2(dim);
        state2.setZero();
        state2(s) = 1.0;

        VectorReal corrLocal(dataPoints);
        corrLocal.setZero();

        double weightLocal = 0.0;

        for (int b = 0; b < Betaiters; b++)
            evolveKrylovImaginary(state1, N, J, hz, hx, dbeta, krylovDim);
        
        applyHIsing(state1, state2, N, 0, 0, 1);

        weightLocal = state1.squaredNorm();
        
        Vector Astate2(dim);
        Astate2.setZero();
        
        for (int iter = 0; iter < iters; iter++)
        {
            //std::cout << "s: " << s << ", iter: " << iter << "/" << iters << "\r" << std::flush;
            if (iter % itersPerDataPoint == 0)
            {
                applyHIsing(state2, Astate2, N, 0.0, 0.0, 1.0);
                corrLocal(iter / itersPerDataPoint) += std::real(state1.dot(Astate2));
            }
            evolveKrylov(state1, N, hbar, J, hz, hx, 0, 0, iter*dt, dt, krylovDim);
            evolveKrylov(state2, N, hbar, J, hz, hx, 0, 0, iter*dt, dt, krylovDim);
            
        }

        #ifdef _OPENMP
        #pragma omp critical
        #endif
        {
            correlations += corrLocal;
            WeightSum += weightLocal;
        }
        
    }
    correlations /= WeightSum;

    Eigen::FFT<double> fft;
    Eigen::VectorXcd currSpectrum;

    Vector frequencies = 2.0 * pi * VectorReal::LinSpaced(dataPoints, 0.0, dataPoints - 1.0) / (dataPoints * dt);
    fft.fwd(currSpectrum, correlations);

    std::ofstream file;
    file.open(fileName);
    if (!file)
    {
        throw std::runtime_error("could not open output file: " + fileName);
    }
    file << std::setprecision(15);

    file << "frequency,spectralFunction";
    for (int i = 1; i < int(dataPoints/2.0); i++)
    {
        file << "\n" << std::real(frequencies(i)) << "," << dt * std::real(currSpectrum(i));
    }

    file.close();
    
}