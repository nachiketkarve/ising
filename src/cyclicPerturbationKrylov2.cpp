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
    int krylovDim = 20;

    if (argc != 3)
    {
        throw std::invalid_argument("usage: ./_cyclicPerturbationFast.out <maxTime> <lambda>");
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
    krylovDim = params["krylovDim"];

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

    VectorReal Eis(dim);
    VectorReal Efs(dim);
    VectorReal Weights(dim);

    std::string fileName = saveFolder + "isingFast_N" + std::to_string(N) + "_lambda" + std::to_string(lambda) + "_T" + std::to_string(Tmax) + "_B" + std::to_string(beta) + ".csv";

    double numerator = 0.0;
    double denominator = 0.0;

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int s = 0; s < dim; s++)
    {
        Vector state(dim);
        state.setZero();
        state(s) = 1.0;

        for (int b = 0; b < Betaiters; b++)
            evolveKrylovImaginary(state, N, J, hz, hx, dbeta, krylovDim);

        double weight = state.squaredNorm();
        double Ei = expectationHz(state, N, J, hz) + expectationHx(state, N, hx);

        for (int iter = 0; iter < iters; iter++)
            evolveKrylov(state, N, hbar, J, hz, hx, lambda, mu, iter*dt, dt, krylovDim);

        double normf = state.squaredNorm();
        double EfRaw = expectationHz(state, N, J, hz) + expectationHx(state, N, hx);
        double Ef = EfRaw * weight / normf;

        Weights(s) = weight;

        numerator += Ef - Ei;
        denominator += weight;

        Eis(s) = Ei;
        Efs(s) = Ef;
    }

    double Work = numerator / denominator;

    std::cout << "Work: "<< Work;

    std::ofstream file;
    file.open(fileName);
    if (!file)
    {
        throw std::runtime_error("could not open output file: " + fileName);
    }
    file << std::setprecision(15);

    file << "Ei,Ef,deltaE,weight";
    for (int i = 0; i < dim; i++)
    {
        file << "\n" << Eis(i) << "," << Efs(i) << "," << Efs(i) - Eis(i) << "," << Weights(i);
    }

    file.close();
    
}