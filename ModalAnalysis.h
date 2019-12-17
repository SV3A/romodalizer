#ifndef MODALANALYSISHEADDEF
#define MODALANALYSISHEADDEF

#include <iostream>
#include <cmath>
#include "Eigen/Dense"

typedef Eigen::Matrix<double, 5, Eigen::Dynamic> ElementsMatrix;

class ModalAnalysis
{
  private:
    unsigned int numEl;    // Number of elements
    unsigned int numDof;   // Number of degress of freedom
    //Eigen::Matrix<double,5,Eigen::Dynamic> *mesh;  // Pointer to mesh

    // Matrices
    Eigen::MatrixXd M;
    Eigen::MatrixXd G;
    Eigen::MatrixXd K;

    void buildShaftMatrices(const ElementsMatrix &elements);

  public:
    ModalAnalysis(const ElementsMatrix &elements);

    void printInfo();
};

#endif
