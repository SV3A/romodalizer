#ifndef MODALANALYSIS_H
#define MODALANALYSIS_H

#include <iostream>
#include <cmath>
#include <tuple>
#include <vector>
#include <complex>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/Eigenvalues"

typedef Eigen::Triplet<double> Triplet;
typedef Eigen::Matrix<double, 5, Eigen::Dynamic> ElementsMatrix;

class ModalAnalysis
{
  private:
    size_t numEl;  // Number of elements
    size_t numDof; // Number of degress of freedom
    float omega;   // Angular velocity [rad/s]

    // Global matrices
    Eigen::MatrixXd M;
    Eigen::MatrixXd G;
    Eigen::MatrixXd K;

    // State matrices
    Eigen::SparseMatrix<double> A, B;

    // List of tuples containing eigenvector and eigenvalue pairs
    std::vector<std::tuple<std::complex<double>, Eigen::VectorXcd>> eigenSolution;

    void buildShaftMatrices(const ElementsMatrix &elements);

    void buildStateSpace();

  public:
    ModalAnalysis(double omega, const ElementsMatrix &elements);

    void solve();

    void printInfo();
};

#endif
