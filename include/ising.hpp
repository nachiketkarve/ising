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


void buildHamiltonian(MatrixReal &H, int N, double J, double hx, double hz);
void buildPerturbation(MatrixReal &V, int N);
void UpdateUnitary(const MatrixReal &H, const MatrixReal &V, Matrix &U, double hbar, double lambda, double mu, double t, double dt);
void expHz(Vector &state, int N, double J, double hz, double hbar, double dt);
void expHzReal(Vector &state, int N, double J, double hz, double dbeta);
void expHx(Vector &state, int N, double hx, double hbar, double dt);
void expHzReal(Vector &state, int N, double J, double hz, double dbeta);
void evolveTrotter(Vector &state, int N, double hbar, double J, double hz, double hx, double lambda, double mu, double time, double dt);
void evolveTrotter4(Vector &state, int N, double hbar, double J, double hz, double hx, double lambda, double mu, double time, double dt);
void evolveTrotterImaginary(Vector &state, int N, double J, double hz, double hx, double dbeta);
double expectationHz(const Vector &state, int N, double J, double hz);
double expectationHx(const Vector &state, int N, double hx);
void applyHIsing(const Vector &psi, Vector &out, int N, double J, double hz, double hxEff);
void evolveKrylov(Vector &psi, int N, double hbar, double J, double hz, double hx, double lambda, double mu, double time, double dt, int krylovDim);
void evolveKrylovImaginary(Vector &psi, int N, double J, double hz, double hx, double dbeta, int krylovDim);

#endif