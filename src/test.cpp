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
    int s = 1;
    double Tmax = 1.0;
    double lambda = 0.1;
    double dt = 0.1;

    std::ifstream dataFile("_params.json");
    json params = json::parse(dataFile);

    hbar = params["hbar"];
    N = params["N"];
    J = params["J"];
    hx = params["hx"];
    hz = params["hz"];
    dt = params["dt"];

    s = std::stoi(argv[1]);
    Tmax = std::stof(argv[2]);
    lambda = std::stof(argv[3]);
    double mu = 2.0 * pi / Tmax;

    int dim = 1 << N;
    std::cout << "Hilbert space dimension: " << dim << std::endl;

    Vector state1(dim);
    state1.setZero();
    state1(s) = 1.0;

    int iters = static_cast<int>(std::round(Tmax / dt));
    if (iters <= 0) iters = 1;
    dt = Tmax / iters;

    for(int iter = 0; iter < iters; iter++)
    {
        evolveTrotter4(state1, N, hbar, J, hz, hx, lambda, mu, dt * iter, dt);
    }

    Vector state2(dim);
    state2.setZero();
    state2(s) = 1.0;

    MatrixReal H(dim, dim);
    MatrixReal V(dim, dim);
    buildHamiltonian(H, N, J, hx, hz);
    buildPerturbation(V, N);

    Matrix U = Matrix::Identity(dim, dim);

    for(int iter = 0; iter < iters; iter++)
    {
        UpdateUnitary(H,V,U,hbar,lambda,mu,dt*iter,dt);
    }
    state2 = U * state2;

    double norm = std::real((state1-state2).norm());

    std::cout << std::setprecision(17);
    std::cout << "Error: " << norm << std::endl;
    
}