//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define the matrix and vector classes
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef HEADERFILE_EIGENCLASSES
#define HEADERFILE_EIGENCLASSES
#include "Eigen/Dense"
#include "unsupported/Eigen/CXX11/Tensor"
#include "unsupported/Eigen/Polynomials"
#include "Eigen/LU"
#include "Eigen/Eigenvalues"

// Define Matrix and Vector classes
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::Vector<std::complex<double>, Eigen::Dynamic> Vector;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixReal;
typedef Eigen::Vector<double, Eigen::Dynamic> VectorReal;

// Define pi
const double pi = 3.14159265358979323846;

#endif