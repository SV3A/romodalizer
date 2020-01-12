#include <iostream>
#include "Mesh.h"
#include "ModalAnalysis.h"
#include "Eigen/Dense"

int main(int argc, char *argv[])
{
  // Define mesh
  double ro = 0.025*0.5;
  Eigen::MatrixXd m(4, 1); m << 200, // Length
                                 ro, // Outer radius
                                0.0, // Inner radius
                                  4; // Partition num

  m.row(0) = m.row(0)*0.001; // Convert from [mm] to [m]

  Mesh *mesh = new Mesh(m);
  mesh->setDensity(2770.0);
  mesh->setEmod((double) 69.0e9);

  ModalAnalysis *modAnalysis = new ModalAnalysis(0.0f, mesh->elements);

  modAnalysis->printInfo();

  modAnalysis->solve();

  // Print eigenvalues
  std::vector<EigTuple> * eigSol = & modAnalysis->eigenSolution;

  std::cout << "Eigenvalues:" << std::endl;
  for (auto const & i:*eigSol){
    std::cout << "   " << std::get<0>(i) << std::endl;
  }

  // Cleanup
  delete mesh;
  delete modAnalysis;
}
