#include "ModalAnalysis.h"

ModalAnalysis::ModalAnalysis(Eigen::MatrixXf &mesh)
{
  this->mesh  = &mesh;
  this->numEl = mesh.cols();

  // Calculate number of degress of freedom
  numDof = (numEl + 1)*4;

  buildMassMatrix();
}

void ModalAnalysis::buildMassMatrix()
{
}

void ModalAnalysis::printInfo()
{
  std::cout << "Number of elements: " << numEl << std::endl;
  std::cout << "Number of DOFs: " << numDof << std::endl;
  std::cout << "Mesh:" << std::endl;
  std::cout << *mesh << std::endl;
}
