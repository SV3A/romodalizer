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
#include "NodeComponents.h"

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

    // Defines the local matrices for each shaft element and inserts the local
    // matrices into the global M, G, and K matrices.
    void buildShaftMatrices(const ElementsMatrix& elements);

    // Buids state matrices A and B as:
    //     | omega*G -K |      | M  0 |
    // A = |    K     0 |, B = | 0  M |
    void buildStateSpace();

  public:
    ModalAnalysis(double omega, const ElementsMatrix& elements);

    ~ModalAnalysis();

    // List of tuples containing eigenvector and eigenvalue pairs
    std::vector<EigTuple>* eigenSolution;

    // Solve the generalized EVP: A*phi_i = lambda_i*B*phi_i
    void solve();

    // Prints summary
    void printInfo() const;

    // Overloads for adding external components to the rotor

    // Method for adding disc elements, i.e. an inertial contribution
    void addNodeComponent(size_t node, Disc& disc);
    // Method for adding a support, i.e. a stiffness contribution
    void addNodeComponent(size_t node, Bearing& bearingElement);
};

#endif
