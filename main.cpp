#include <iostream>
#include "Mesh.h"
#include "ModalAnalysis.h"
#include "Eigen/Dense"

int main(int argc, char *argv[])
{
  // Define mesh
  double ro = 0.025*0.5;
  Eigen::MatrixXd m(4, 5); m << 162,  88, 120,  52,  63, // Length
                                 ro,  ro,  ro,  ro,  ro, // Outer radius
                                0.0, 0.0, 0.0, 0.0, 0.0, // Inner radius
                                  3,   2,   3,   1,   2; // Partition num

  m.row(0) = m.row(0)*0.001; // Convert from [mm] to [m]

  Mesh *mesh = new Mesh(m);
  mesh->setDensity(2770.0);
  mesh->setEmod((double) 69.0e9);

  ModalAnalysis *modAnalysis = new ModalAnalysis(mesh->elements);

  modAnalysis->printInfo();

  // Cleanup
  delete mesh;
  delete modAnalysis;
}
