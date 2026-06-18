#include "ising.hpp"
#include "nlohmann/json.hpp"
#include <random>

using json = nlohmann::json;

int main(int argc, char *argv[])
{
    double hbar = 1.0;
    int N = 5;
    double J = 1.0;
    double hx = 1.0;
    double hz = 1.0;
    double Tmax = 1.0;
    double lambda = 0.1;
    double dt = 0.1;
    int seed = 1;
    double beta = 1.0;
    int averages = 1;

    if (argc != 3)
    {
        throw std::invalid_argument("usage: ./build/bin/cyclicPerturbationFast.exe <maxTime> <lambda>");
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
    averages = params["averages"];
    std::string saveFolder = params["saveFolder"];

    double dbeta = dt;
    int Betaiters = static_cast<int>(std::round(beta / (dbeta * 2.0)));
    if (Betaiters <= 0) Betaiters = 1;
    dbeta = beta / (Betaiters * 2.0);

    Tmax = std::stof(argv[1]);
    lambda = std::stof(argv[2]);
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

    int dim = 1 << N;
    std::cout << "Hilbert space dimension: " << dim << std::endl;

    VectorReal Ei(averages);
    VectorReal Ef(averages);
    VectorReal Weights(averages);

    std::string fileName = saveFolder + "isingFast_N" + std::to_string(N) + "_lambda" + std::to_string(lambda) + "_T" + std::to_string(Tmax) + "_B" + std::to_string(beta) + ".csv";

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for(int average = 0; average < averages; average++)
    {
        std::mt19937                        generator(seed + average);
        std::normal_distribution<double>  distr(0.0,1.0);

        Vector state(dim);
        state.setZero();
        for (int s = 0; s < dim; s++)
        {
            state(s) = std::complex<double>(distr(generator),distr(generator));
        }

        state = state / std::sqrt(std::real(state.dot(state)));

        for(int Betaiter = 0; Betaiter < Betaiters; Betaiter++)
        {
            evolveTrotterImaginary(state, N, J, hz, hx, dbeta);
        }

        Weights(average) = std::real(state.dot(state));
        Ei(average) = (expectationHz(state, N, J, hz) + expectationHx(state, N, hx));

        for(int iter = 0; iter < iters; iter++)
        {
            evolveTrotter(state, N, hbar, J, hz, hx, lambda, mu, dt * iter, dt);
        }

        Ef(average) = (expectationHz(state, N, J, hz) + expectationHx(state, N, hx)) / std::real(state.dot(state)) * Weights(average);
    }

    std::cout << "Work: "<< (Ef-Ei).sum()/Weights.sum();

    std::ofstream file;
    file.open(fileName);
    if (!file)
    {
        throw std::runtime_error("could not open output file: " + fileName);
    }
    file << std::setprecision(15);

    file << "Ei,Ef,deltaE,weight";
    for (int i = 0; i < averages; i++)
    {
        file << "\n" << Ei(i) << "," << Ef(i) << "," << Ef(i) - Ei(i) << "," << Weights(i);
    }

    file.close();
    
}
