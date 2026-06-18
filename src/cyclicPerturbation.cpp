#include "ising.hpp"
#include "nlohmann/json.hpp"
#include <cstdlib>

using json = nlohmann::json;

int main(int argc, char *argv[])
{
    int N = 5;
    double hbar = 1.0;
    double J = 1.0;
    double hx = 1.0;
    double hz = 1.0;
    double dt = 0.1;
    double lambda = 0.1;
    double T = 1.0;
    std::string saveFolder = "./";

    if (argc != 3)
    {
        throw std::invalid_argument("usage: ./build/bin/cyclicPerturbation.exe <maxTime> <lambda>");
    }

    std::ifstream dataFile("_params.json");
    if (!dataFile)
    {
        throw std::runtime_error("could not open _params.json");
    }

    json params = json::parse(dataFile);
    N = params["N"];
    hbar = params["hbar"];
    J = params["J"];
    hx = params["hx"];
    hz = params["hz"];
    dt = params["dt"];
    saveFolder = params["saveFolder"];

    T = std::atof(argv[1]);
    lambda = std::atof(argv[2]);
    if (N < 0)
    {
        throw std::invalid_argument("S must be non-negative");
    }
    if (hbar <= 0.0)
    {
        throw std::invalid_argument("hbar must be positive");
    }
    if (T <= 0.0 || dt <= 0.0)
    {
        throw std::invalid_argument("T and dt must be positive");
    }

    int Nsteps = static_cast<int>(T / dt);
    if (Nsteps <= 0)
    {
        throw std::invalid_argument("T must be at least one dt");
    }
    dt = T / Nsteps; // Adjust dt to fit exactly into T
    double mu = 2.0 * pi / T;

    std::string fileName = saveFolder + "ising_N" + std::to_string(N) + "_lambda" + std::to_string(lambda) + "_T" + std::to_string(T) + ".csv";

    int dim = std::pow(2,N);
    MatrixReal H(dim, dim);
    MatrixReal V(dim, dim);
    buildHamiltonian(H, N, J, hx, hz);
    buildPerturbation(V, N);

    Eigen::SelfAdjointEigenSolver<MatrixReal> eigSolver(H);
    VectorReal eigVals = eigSolver.eigenvalues();
    Matrix eigVecs = eigSolver.eigenvectors();

    Matrix U = Matrix::Identity(dim, dim);
    for (int step = 0; step < Nsteps; step++)
    {
        
        UpdateUnitary(H, V, U, hbar, lambda, mu, step * dt, dt);
        
    }

    Matrix I = Matrix::Identity(dim, dim);
    std::cout << "Unitarity error: " << (U.adjoint() * U - I).norm() << std::endl;

    Matrix H_eig = eigVecs.adjoint() * H * eigVecs;
    
    Matrix D = Matrix::Zero(dim, dim);
    for (int i = 0; i < dim; i++)
    {
        D(i, i) = std::complex<double>(eigVals(i), 0.0);
    }

    double diagonalization_error = (H_eig - D).norm();
    std::cout << "diagonalization_error = " << diagonalization_error << std::endl;

    VectorReal Es(dim);
    for (int i = 0; i < dim; i++)
    {
        
        Vector psi_i = eigVecs.col(i);
        Vector psi_f = U * psi_i;

        double norm = std::real(psi_f.dot(psi_f));
        double Ef = std::real(psi_f.dot(H * psi_f)) / norm;

        Es(i) = Ef;
    }

    Vector psi0 = eigVecs.col(0);
    Vector psi0f = U * psi0;

    double norm0 = std::real(psi0f.dot(psi0f));
    double Ef0_direct = std::real(psi0f.dot(H * psi0f)) / norm0;
    double ground_delta_direct = Ef0_direct - eigVals(0);

    std::cout << std::setprecision(17);

    std::cout << "E0 = " << eigVals(0) << std::endl;
    std::cout << "Ef0 direct = " << Ef0_direct << std::endl;
    std::cout << "Ef0 - E0 direct = " << ground_delta_direct << std::endl;
    std::cout << "norm0 direct = " << norm0 << std::endl;


    std::ofstream file;
    file.open(fileName);
    if (!file)
    {
        throw std::runtime_error("could not open output file: " + fileName);
    }
    file << std::setprecision(15);

    file << "Ei,Ef,deltaE";
    for (int i = 0; i < dim; i++)
    {
        file << "\n" << eigVals(i) << "," << Es(i) << "," << Es(i) - eigVals(i);
    }

    file.close();





}
