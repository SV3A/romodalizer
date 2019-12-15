#ifndef MODALANALYSISHEADDEF
#define MODALANALYSISHEADDEF

#include <iostream>
#include "Eigen/Eigenvalues"

class ModalAnalysis
{
  private:
    unsigned int numEl;    // Number of elements
    unsigned int numDof;   // Number of degress of freedom
    Eigen::MatrixXf *mesh; // Pointer to mesh

    // Matrices
    Eigen::MatrixXf mass;
    Eigen::MatrixXf damp;
    Eigen::MatrixXf stiff;

    void buildMassMatrix();

  public:
    ModalAnalysis(Eigen::MatrixXf &mesh);

    void printInfo();
};

#endif
