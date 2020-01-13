#ifndef MODALANALYSIS_H
#define MODALANALYSIS_H

#include <iostream>
#include <cmath>
#include <tuple>
#include <vector>
#include <complex>
#include <algorithm>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/Eigenvalues"
#include "Util.h"

typedef Eigen::Triplet<double> Triplet;
typedef Eigen::Matrix<double, 5, Eigen::Dynamic> ElementsMatrix;
typedef std::tuple<std::complex<double>, Eigen::VectorXcd> EigTuple;

class ModalAnalysis
{
  private:
    const size_t numEl;  // Number of elements
    const size_t numDof; // Number of degress of freedom
    float omega;         // Angular velocity [rad/s]

    // Global matrices
    Eigen::MatrixXd M;
    Eigen::MatrixXd G;
    Eigen::MatrixXd K;

    // State matrices
    Eigen::SparseMatrix<double> A, B;

    void buildShaftMatrices(const ElementsMatrix& elements);

    void buildStateSpace();

  public:
    ModalAnalysis(double omega, const ElementsMatrix& elements);

    ~ModalAnalysis();

    void solve();

    void printInfo() const;

    // List of tuples containing eigenvector and eigenvalue pairs
    std::vector<EigTuple>* eigenSolution;
};

#endif
