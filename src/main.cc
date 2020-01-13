#include <iostream>
#include "Mesh.h"
#include "ModalAnalysis.h"
#include "Eigen/Dense"

int main(int argc, char *argv[])
{
  // Define mesh
  double ro = 0.025*0.5;
  Eigen::MatrixXd m(4, 1); m << 1000, // Length [mm]
                                  ro, // Outer radius
                                 0.0, // Inner radius
                                   4; // Partition num

  m.row(0) = m.row(0)*0.001; // Convert from [mm] to [m]

  Mesh *mesh = new Mesh(m);
  mesh->setDensity(2770.0);
  mesh->setEmod((double) 69.0e9);

  ModalAnalysis* modAnalysis = new ModalAnalysis(0.0f, mesh->elements);

  modAnalysis->printInfo();

  modAnalysis->solve();

  // Print eigenvalues and natural frequencies
  std::vector<EigTuple>* eigSol = modAnalysis->eigenSolution;

  //std::cout << "Eigenvalues:" << std::endl;
  //for (auto const & i:*eigSol){
    //printf("\t%.9f\n", imag(std::get<0>(i)));
  //}

  std::cout << "Natural frequencies:" << std::endl;
  for (auto const & i:*eigSol){
    printf("   %10.3f Hz\n", imag(std::get<0>(i))/2.0/M_PI);
  }

  // Cleanup
  delete mesh;
  delete modAnalysis;
}
