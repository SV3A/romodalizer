#include "ModalAnalysis.h"

// Constructor
ModalAnalysis::ModalAnalysis(const ElementsMatrix &elements)
{
  // Get number of elements and derrive from it the number of D.O.F.s
  numEl  = elements.cols();
  numDof = (numEl + 1)*4;

  // Set size of system matrices
  M.resize(numDof, numDof);
  D.resize(numDof, numDof);
  K.resize(numDof, numDof);

  // Build matrices
  buildMassMatrix(elements);
  //buildDampMatrix(l);
  //buildStiffMatrix(l);
}


// Defines the local mass matrices (linear and rotational) for each element and
// inserts the these local matrices into the global mass matrix M.
void ModalAnalysis::buildMassMatrix(const ElementsMatrix &elements)
{
  // Start and end indices
  unsigned int a = 0, b = 7;

  double l, lsq, rho, transArea, ro, ri;

  // Temp matrices
  Eigen::MatrixXd linMassAux(8, 8), rotMassAux(8, 8);

  for (int e = 0; e < numEl; e++) {
    l   = elements(0, e);
    ro  = elements(1, e);
    ri  = elements(2, e);
    rho = elements(3, e);
    lsq = l*l;
    transArea = M_PI*(pow(ro,2) - pow(ri,2));

    // Linear inertia matrix
    linMassAux <<
    156,   0,     0,      22*l,   54,    0,     0,     -13*l,
    0,     156,  -22*l,   0,      0,     54,    13*l,   0,
    0,    -22*l,  4*lsq,  0,      0,    -13*l, -3*lsq,  0,
    22*l,  0,     0,      4*lsq,  13*l,  0,     0,     -3*l*l,
    54,    0,     0,      13*l,   156,   0,     0,     -22*l,
    0,     54,   -13*l,   0,      0,     156,   22*l,   0,
    0,     13*l, -3*lsq,  0,      0,     22*l,  4*lsq,  0,
   -13*l,  0,     0,     -3*lsq, -22*l,  0,     0,      4*lsq;

    linMassAux *= (rho*transArea*l) / 420.0;

    // Angular inertia matrix
    rotMassAux <<
    36,   0,    0,      3*l,   -36,   0,    0,      3*l,
    0,    36,  -3*l,    0,      0,   -36,  -3*l,    0,
    0,   -3*l,  4*lsq,  0,      0,    3*l, -lsq,    0,
    3*l,  0,    0,      4*lsq, -3*l,  0,    0,     -lsq,
   -36,   0,    0,     -3*l,    36,   0,    0,     -3*l,
    0,   -36,   3*l,    0,      0,    36,   3*l,    0,
    0,   -3*l, -lsq,    0,      0,    3*l,  4*lsq,  0,
    3*l,  0,    0,     -lsq,   -3*l,  0,    0,      4*lsq;

    rotMassAux *= (rho*transArea*(ro*ro - ri*ri)) / (120.0*l);

    // Add the matrices
    linMassAux += rotMassAux;

    // Construct the global mass matrix (of size numDofxnumDof)
    for (int i = a; i <= b; ++i) {
      for (int j = a; j <= b; ++j) {
        M(i, j) += linMassAux(i - e*4, j - e*4);
      }
    }

    a += 4; b += 4;
  }
}

void ModalAnalysis::printInfo()
{
  std::cout << "Number of elements: " << numEl << std::endl;
  std::cout << "Number of DOFs: " << numDof << std::endl;
}
