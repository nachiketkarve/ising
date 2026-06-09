#include "ising.hpp"

int bit_i(int s, int i)
{
    return (s >> i) & 1;
}

int sigmaz_i(int s, int i)
{
    int si = bit_i(s,i);
    if (si == 1) return 1;
    else return -1;
}

int flip_spin_i(int s, int i)
{
    return s ^ (1 << i);
}

void buildHamiltonian(MatrixReal &H, int N, double J, double hx, double hz)
{
    int dim = std::pow(2,N);

    if (H.rows() != dim || H.cols() != dim) return;

    H.setZero();

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for(int s = 0; s < dim; s++)
    {
        for(int i = 0; i < N; i++)
        {
            int iNext = (i+1)%N;
            H(s,s) = H(s,s) + J * float(sigmaz_i(s,i) * sigmaz_i(s,iNext)) + hz * float(sigmaz_i(s,i));
            int sFlip = flip_spin_i(s,i);
            H(s,sFlip) = H(s,sFlip) + hx;
        }
    }
}

void buildPerturbation(MatrixReal &V, int N)
{
    int dim = std::pow(2,N);

    if (V.rows() != dim || V.cols() != dim) return;

    V.setZero();

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for(int s = 0; s < dim; s++)
    {
        for(int i = 0; i < N; i++)
        {
            int sFlip = flip_spin_i(s,i);
            V(s,sFlip) = V(s,sFlip) + 1.0;
        }
    }
}

void UpdateUnitary(const MatrixReal &H, const MatrixReal &V, Matrix &U, double hbar, double lambda, double mu, double t, double dt)
{
    double tm = t + 0.5 * dt;
    MatrixReal Hmid = H + lambda * std::pow(std::sin(mu * tm), 2) * V;

    Eigen::SelfAdjointEigenSolver<MatrixReal> eig(Hmid);
    VectorReal evals = eig.eigenvalues();
    Matrix evecs = eig.eigenvectors();

    Matrix phase = Matrix::Zero(H.rows(), H.cols());
    for (int i = 0; i < evals.size(); ++i)
    {
        phase(i, i) = std::exp(-std::complex<double>(0, 1) * evals(i) * dt / hbar);
    }

    Matrix Ustep = evecs * phase * evecs.adjoint();
    U = Ustep * U;
}


void expHz(Vector &state, int N, double J, double hz, double hbar, double dt)
{
    const int dim = 1 << N;
    if (state.size() != dim) return;

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int s = 0; s < dim; s++)
    {
        double phase = 0;
        for(int i = 0; i < N; i++)
        {
            int iNext = (i+1)%N;
            phase = phase + J * float(sigmaz_i(s,i) * sigmaz_i(s,iNext)) + hz * float(sigmaz_i(s,i));
        }
        phase = - phase * dt / hbar;
        state(s) = state(s) * std::complex<double>(std::cos(phase),std::sin(phase));
    }
}

void expHzReal(Vector &state, int N, double J, double hz, double dbeta)
{
    const int dim = 1 << N;
    if (state.size() != dim) return;

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int s = 0; s < dim; s++)
    {
        double phase = 0;
        for(int i = 0; i < N; i++)
        {
            int iNext = (i+1)%N;
            phase = phase + J * float(sigmaz_i(s,i) * sigmaz_i(s,iNext)) + hz * float(sigmaz_i(s,i));
        }
        phase = - phase * dbeta;
        state(s) = state(s) * std::exp(phase);
    }
}

void expHx(Vector &state, int N, double hx, double hbar, double dt)
{
    const int dim = 1 << N;
    if (state.size() != dim) return;

    const double theta = hx * dt / hbar;
    const double c = std::cos(theta);
    const std::complex<double> minus_i_s(0.0, -std::sin(theta));

    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    for (int i = 0; i < N; i++)
    {
        const int mask = 1 << i;

        #ifdef _OPENMP
        #pragma omp for schedule(static)
        #endif
        for (int s = 0; s < dim; s++)
        {
            if ((s & mask) == 0)
            {
                int sFlip = s | mask;

                std::complex<double> a = state(s);
                std::complex<double> b = state(sFlip);

                state(s)     = c * a + minus_i_s * b;
                state(sFlip) = minus_i_s * a + c * b;
            }
        }
    }
}

void expHxReal(Vector &state, int N, double hx, double dbeta)
{
    const int dim = 1 << N;
    if (state.size() != dim) return;

    const double theta = hx * dbeta;
    const double ch = std::cosh(theta);
    const double sh = -std::sinh(theta);

    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    for (int i = 0; i < N; i++)
    {
        const int mask = 1 << i;

        #ifdef _OPENMP
        #pragma omp for schedule(static)
        #endif
        for (int s = 0; s < dim; s++)
        {
            if ((s & mask) == 0)
            {
                int sFlip = s | mask;

                std::complex<double> a = state(s);
                std::complex<double> b = state(sFlip);

                state(s)     = ch * a + sh * b;
                state(sFlip) = sh * a + ch * b;
            }
        }
    }
}

void evolveTrotter(Vector &state, int N, double hbar, double J, double hz, double hx, double lambda, double mu, double time, double dt)
{
    const int dim = 1 << N;
    if (state.size() != dim) return;

    if ((std::abs(lambda) < TOL) && (std::abs(hx) < TOL))
    {
        expHz(state, N, J, hz, hbar, dt);
        return;
    }

    double tmid = time + dt / 2.0;
    double hxmid = hx + lambda * std::sin(mu * tmid) * std::sin(mu * tmid);

    expHz(state, N, J, hz, hbar, dt / 2.0);
    expHx(state, N, hxmid, hbar, dt);
    expHz(state, N, J, hz, hbar, dt / 2.0);
}

void evolveTrotter4(Vector &state, int N, double hbar, double J, double hz, double hx, double lambda, double mu, double time, double dt)
{
    const int dim = 1 << N;
    if (state.size() != dim) return;

    const double p = 1.0 / (4.0 - std::cbrt(4.0));
    const double q = 1.0 - 4.0 * p;

    double t0 = time;

    evolveTrotter(state, N, hbar, J, hz, hx, lambda, mu, t0, p * dt);
    t0 += p * dt;

    evolveTrotter(state, N, hbar, J, hz, hx, lambda, mu, t0, p * dt);
    t0 += p * dt;

    evolveTrotter(state, N, hbar, J, hz, hx, lambda, mu, t0, q * dt);
    t0 += q * dt;

    evolveTrotter(state, N, hbar, J, hz, hx, lambda, mu, t0, p * dt);
    t0 += p * dt;

    evolveTrotter(state, N, hbar, J, hz, hx, lambda, mu, t0, p * dt);
}

void evolveTrotterImaginary(Vector &state, int N, double J, double hz, double hx, double dbeta)
{
    const int dim = 1 << N;
    if (state.size() != dim) return;

    expHzReal(state, N, J, hz, dbeta / 2.0);
    expHxReal(state, N, hx, dbeta);
    expHzReal(state, N, J, hz, dbeta / 2.0);
}

double expectationHz(const Vector &state, int N, double J, double hz)
{
    const int dim = 1 << N;

    if (state.size() != dim)
    {
        throw std::invalid_argument("state has wrong dimension");
    }

    double expHz = 0.0;

    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:expHz) schedule(static)
    #endif
    for (int s = 0; s < dim; s++)
    {
        double Ez = 0.0;

        for (int i = 0; i < N; i++)
        {
            int iNext = (i + 1) % N;

            Ez += J * sigmaz_i(s, i) * sigmaz_i(s, iNext);
            Ez += hz * sigmaz_i(s, i);
        }

        expHz += std::norm(state(s)) * Ez;
    }

    return expHz;
}

double expectationHx(const Vector &state, int N, double hx)
{
    const int dim = 1 << N;

    if (state.size() != dim)
    {
        throw std::invalid_argument("state has wrong dimension");
    }

    double expHx = 0.0;

    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:expHx) schedule(static)
    #endif
    for (int s = 0; s < dim; s++)
    {
        double local = 0.0;

        for (int i = 0; i < N; i++)
        {
            int sFlip = s ^ (1 << i);
            local += std::real(std::conj(state(s)) * state(sFlip));
        }

        expHx += local;
    }
    return hx * expHx;
}


////////////////////////////////////////////////////////////////////////////////////////////

void applyHIsing(const Vector &psi, Vector &out, int N, double J, double hz, double hxEff)
{
    const int dim = 1 << N;

    if (psi.size() != dim || out.size() != dim)
    {
        throw std::invalid_argument("wrong vector size in applyHIsing");
    }

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int s = 0; s < dim; s++)
    {
        double Ez = 0.0;

        for (int i = 0; i < N; i++)
        {
            int iNext = (i + 1) % N;

            int zi = sigmaz_i(s, i);
            int zj = sigmaz_i(s, iNext);

            Ez += J * zi * zj;
            Ez += hz * zi;
        }

        std::complex<double> val = Ez * psi(s);

        for (int i = 0; i < N; i++)
        {
            int sFlip = s ^ (1 << i);
            val += hxEff * psi(sFlip);
        }

        out(s) = val;
    }
}

void evolveKrylov(Vector &psi, int N, double hbar, double J, double hz, double hx, double lambda, double mu, double time, double dt, int krylovDim)
{
    const int dim = 1 << N;

    if (psi.size() != dim)
    {
        throw std::invalid_argument("wrong vector size in evolveKrylov");
    }

    if (krylovDim <= 0)
    {
        throw std::invalid_argument("krylovDim must be positive");
    }

    const double tmid = time + 0.5 * dt;
    const double sinmt = std::sin(mu * tmid);
    const double hxEff = hx + lambda * sinmt * sinmt;

    const double psiNorm = psi.norm();

    if (psiNorm == 0.0)
    {
        throw std::runtime_error("zero vector in evolveKrylov");
    }

    const int mMax = std::min(krylovDim, dim);

    Matrix Q(dim, mMax);
    VectorReal alpha = VectorReal::Zero(mMax);
    VectorReal beta = VectorReal::Zero(mMax);

    Vector qPrev = Vector::Zero(dim);
    Vector q = psi / psiNorm;
    Vector w(dim);

    int mActual = 0;

    const double lanczosTol = 1e-13;

    for (int j = 0; j < mMax; j++)
    {
        Q.col(j) = q;

        applyHIsing(q, w, N, J, hz, hxEff);

        if (j > 0)
        {
            w -= beta(j - 1) * qPrev;
        }

        alpha(j) = std::real(q.dot(w));

        w -= alpha(j) * q;

        // Full reorthogonalization for numerical stability.
        // This is useful for moderate Krylov dimensions, e.g. 20--80.
        for (int k = 0; k <= j; k++)
        {
            std::complex<double> coeff = Q.col(k).dot(w);
            w -= coeff * Q.col(k);
        }

        beta(j) = w.norm();

        mActual = j + 1;

        if (beta(j) < lanczosTol)
        {
            break;
        }

        qPrev = q;
        q = w / beta(j);
    }

    MatrixReal T = MatrixReal::Zero(mActual, mActual);

    for (int j = 0; j < mActual; j++)
    {
        T(j, j) = alpha(j);

        if (j + 1 < mActual)
        {
            T(j, j + 1) = beta(j);
            T(j + 1, j) = beta(j);
        }
    }

    Eigen::SelfAdjointEigenSolver<MatrixReal> eig(T);

    if (eig.info() != Eigen::Success)
    {
        throw std::runtime_error("Krylov tridiagonal diagonalization failed");
    }

    VectorReal evals = eig.eigenvalues();
    MatrixReal evecsReal = eig.eigenvectors();

    Vector e0 = Vector::Zero(mActual);
    e0(0) = 1.0;

    Matrix evecs = evecsReal.cast<std::complex<double>>();

    Vector coeffs = evecs.adjoint() * e0;

    for (int j = 0; j < mActual; j++)
    {
        coeffs(j) *= std::exp(
            std::complex<double>(0.0, -evals(j) * dt / hbar)
        );
    }

    Vector y = evecs * coeffs;

    psi = psiNorm * Q.leftCols(mActual) * y;
}

void evolveKrylovImaginary(Vector &psi, int N, double J, double hz, double hx, double dbeta, int krylovDim)
{
    const int dim = 1 << N;

    if (psi.size() != dim)
    {
        throw std::invalid_argument("wrong vector size in evolveKrylovImaginary");
    }

    if (krylovDim <= 0)
    {
        throw std::invalid_argument("krylovDim must be positive");
    }

    if (dbeta < 0.0)
    {
        throw std::invalid_argument("dbeta must be non-negative for imaginary time");
    }

    const double psiNorm = psi.norm();

    if (psiNorm == 0.0)
    {
        throw std::runtime_error("zero vector in evolveKrylovImaginary");
    }

    const int mMax = std::min(krylovDim, dim);

    Matrix Q(dim, mMax);
    VectorReal alpha = VectorReal::Zero(mMax);
    VectorReal beta = VectorReal::Zero(mMax);

    Vector qPrev = Vector::Zero(dim);
    Vector q = psi / psiNorm;
    Vector w(dim);

    int mActual = 0;

    const double lanczosTol = 1e-13;

    for (int j = 0; j < mMax; j++)
    {
        Q.col(j) = q;

        // Static Hamiltonian H0: lambda = 0, hxEff = hx
        applyHIsing(q, w, N, J, hz, hx);

        if (j > 0)
        {
            w -= beta(j - 1) * qPrev;
        }

        alpha(j) = std::real(q.dot(w));

        w -= alpha(j) * q;

        // Full reorthogonalization for stability
        for (int k = 0; k <= j; k++)
        {
            std::complex<double> coeff = Q.col(k).dot(w);
            w -= coeff * Q.col(k);
        }

        beta(j) = w.norm();

        mActual = j + 1;

        if (beta(j) < lanczosTol)
        {
            break;
        }

        qPrev = q;
        q = w / beta(j);
    }

    MatrixReal T = MatrixReal::Zero(mActual, mActual);

    for (int j = 0; j < mActual; j++)
    {
        T(j, j) = alpha(j);

        if (j + 1 < mActual)
        {
            T(j, j + 1) = beta(j);
            T(j + 1, j) = beta(j);
        }
    }

    Eigen::SelfAdjointEigenSolver<MatrixReal> eig(T);

    if (eig.info() != Eigen::Success)
    {
        throw std::runtime_error("Krylov tridiagonal diagonalization failed");
    }

    VectorReal evals = eig.eigenvalues();
    MatrixReal evecsReal = eig.eigenvectors();

    Matrix evecs = evecsReal.cast<std::complex<double>>();

    Vector e0 = Vector::Zero(mActual);
    e0(0) = 1.0;

    Vector coeffs = evecs.adjoint() * e0;

    for (int j = 0; j < mActual; j++)
    {
        coeffs(j) *= std::exp(-evals(j) * dbeta);
    }

    Vector y = evecs * coeffs;

    psi = psiNorm * (Q.leftCols(mActual) * y);
}