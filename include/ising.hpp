#ifndef HEADER_ISING
#define HEADER_ISING

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>
#include "eigenClasses.hpp"

const double TOL = 1e-12;

// Return the ith bit from the right in the state s
int bit_i(int s, int i);

// Return eigenvalue of sigmaz at the ith site 
int sigmaz_i(int s, int i);

// Flip spin at ith site
int flip_spin_i(int s, int i);

// Build the Ising Hamiltonian
void buildHamiltonian(MatrixReal &H, int N, double J, double hx, double hz);

// Build the perturbation matrix
void buildPerturbation(MatrixReal &V, int N);

// Update the unitary matrix
void UpdateUnitary(const MatrixReal &H, const MatrixReal &V, Matrix &U, double hbar, double lambda, double mu, double t, double dt);

// Apply the exponential of the z-component of the Ising Hamiltonian
void expHz(Vector &state, int N, double J, double hz, double hbar, double dt);

// Apply the exponential of the z-component of the Ising Hamiltonian (imaginary time)
void expHzReal(Vector &state, int N, double J, double hz, double dbeta);

// Apply the exponential of the x-component of the Ising Hamiltonian
void expHx(Vector &state, int N, double hx, double hbar, double dt);

// Apply the exponential of the x-component of the Ising Hamiltonian (imaginary time)
void expHxReal(Vector &state, int N, double hx, double dbeta);

// Evolve the state using the Trotter decomposition
void evolveTrotter(Vector &state, int N, double hbar, double J, double hz, double hx, double lambda, double mu, double time, double dt);

// Evolve the state using the 4th-order Trotter decomposition
void evolveTrotter4(Vector &state, int N, double hbar, double J, double hz, double hx, double lambda, double mu, double time, double dt);

// Evolve the state using the imaginary-time Trotter decomposition
void evolveTrotterImaginary(Vector &state, int N, double J, double hz, double hx, double dbeta);

// Calculate the expectation value of the z-component of the Ising Hamiltonian
double expectationHz(const Vector &state, int N, double J, double hz);

// Calculate the expectation value of the x-component of the Ising Hamiltonian
double expectationHx(const Vector &state, int N, double hx);

// Apply the Ising Hamiltonian to a state
void applyHIsing(const Vector &psi, Vector &out, int N, double J, double hz, double hxEff);

// Evolve the state using the Krylov method
void evolveKrylov(Vector &psi, int N, double hbar, double J, double hz, double hx, double lambda, double mu, double time, double dt, int krylovDim);

// Evolve the state using the Krylov method (imaginary time)
void evolveKrylovImaginary(Vector &psi, int N, double J, double hz, double hx, double dbeta, int krylovDim);

#endif