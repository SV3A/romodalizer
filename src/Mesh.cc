#include "Mesh.h"

Mesh::Mesh(Eigen::MatrixXd & mesh)
{
  numEl = mesh.row(3).sum();
  elements.resize(5, numEl);
  setGeometry(mesh);
}


void Mesh::setDensity(double rho, unsigned int element){
  if (element != 0)
    elements(3, element) = rho;
  else
    elements.row(3).setConstant(rho);
}


void Mesh::setEmod(double eMod, unsigned int element){
  if (element != 0)
    elements(4, element) = eMod;
  else
    elements.row(4).setConstant(eMod);
}


void Mesh::setGeometry(Eigen::MatrixXd &mesh){
  double elLength, elOutRadius, elInRadius;
  unsigned int startIdx = 0, endIdx;

  for (int i = 0; i < mesh.cols(); ++i) {
    endIdx = mesh(3,i);

    elLength    = mesh(0, i)/mesh(3, i);
    elOutRadius = mesh(1, i);
    elInRadius  = mesh(2, i);

    elements.row(0).segment(startIdx, endIdx).setConstant(elLength);
    elements.row(1).segment(startIdx, endIdx).setConstant(elOutRadius);
    elements.row(2).segment(startIdx, endIdx).setConstant(elInRadius);

    startIdx += endIdx;
  }
}
