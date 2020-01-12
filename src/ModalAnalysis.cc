#include "ModalAnalysis.h"

// Constructor
ModalAnalysis::ModalAnalysis(double omega, const ElementsMatrix &elements)
{
  // Set angular velocity
  this->omega = omega;

  // Get number of elements and derrive from it the number of D.O.F.s
  numEl  = elements.cols();
  numDof = (numEl + 1)*4;

  // Set size of system matrices and state matrices
  M.resize(numDof, numDof);
  G.resize(numDof, numDof);
  K.resize(numDof, numDof);
  A.resize(numDof*2, numDof*2);
  B.resize(numDof*2, numDof*2);

  // Build shaft matrices
  buildShaftMatrices(elements);
}


// Defines the local matrices for each shaft element and inserts the local
// matrices into the global M, G, and K matrices.
void ModalAnalysis::buildShaftMatrices(const ElementsMatrix &elements)
{
  // Start- and end indices
  size_t a = 0, b = 7;

  double l, lsq, rho, eMod, transArea, momInert, ro, ri;
  double momInertFact = M_PI/4.0;

  // Temp matrices
  Eigen::Matrix<double, 8, 8> localMLin, localMRot, localG, localK;

  for (size_t e = 0; e < numEl; e++) {
    l    = elements(0, e);
    ro   = elements(1, e);
    ri   = elements(2, e);
    rho  = elements(3, e);
    eMod = elements(4, e);

    lsq       = l*l;
    transArea = M_PI*(pow(ro,2) - pow(ri,2));
    momInert  = momInertFact*(pow(ro,4) - pow(ri,4));

    // Mass matrices
    // Linear inertia matrix
    localMLin <<
    156,   0,     0,      22*l,   54,    0,     0,     -13*l,
    0,     156,  -22*l,   0,      0,     54,    13*l,   0,
    0,    -22*l,  4*lsq,  0,      0,    -13*l, -3*lsq,  0,
    22*l,  0,     0,      4*lsq,  13*l,  0,     0,     -3*l*l,
    54,    0,     0,      13*l,   156,   0,     0,     -22*l,
    0,     54,   -13*l,   0,      0,     156,   22*l,   0,
    0,     13*l, -3*lsq,  0,      0,     22*l,  4*lsq,  0,
   -13*l,  0,     0,     -3*lsq, -22*l,  0,     0,      4*lsq;

    localMLin *= (rho*transArea*l) / 420.0;

    // Angular inertia matrix
    localMRot <<
    36,   0,    0,      3*l,   -36,   0,    0,      3*l,
    0,    36,  -3*l,    0,      0,   -36,  -3*l,    0,
    0,   -3*l,  4*lsq,  0,      0,    3*l, -lsq,    0,
    3*l,  0,    0,      4*lsq, -3*l,  0,    0,     -lsq,
   -36,   0,    0,     -3*l,    36,   0,    0,     -3*l,
    0,   -36,   3*l,    0,      0,    36,   3*l,    0,
    0,   -3*l, -lsq,    0,      0,    3*l,  4*lsq,  0,
    3*l,  0,    0,     -lsq,   -3*l,  0,    0,      4*lsq;

    localMRot *= (rho*transArea*(ro*ro - ri*ri)) / (120.0*l);

    // Add the inertia matrices
    localMLin += localMRot;


    // Gyro matrix
    localG <<
    0,   -36,   3*l,    0,     0,    36,   3*l,    0,
    36,   0,    0,      3*l,  -36,   0,    0,      3*l,
   -3*l,  0,    0,     -4*lsq, 3*l,  0,    0,      lsq,
    0,   -3*l,  4*lsq,  0,     0,    3*l, -lsq,    0,
    0,    36,  -3*l,    0,     0,   -36,  -3*l,    0,
   -36,   0,    0,     -3*l,   36,   0,    0,     -3*l,
   -3*l,  0,    0,      lsq,   3*l,  0,    0,     -4*lsq,
    0,   -3*l, -lsq,    0,     0,    3*l,  4*lsq,  0;

    localG *= 2.0*( rho*transArea*(ro*ro + ri*ri) / (120.0*l) );


    // Stiffness matrix
    localK <<
    12,   0,    0,      6*l,   -12,   0,    0,      6*l,
    0,    12,  -6*l,    0,      0,   -12,  -6*l,    0,
    0,   -6*l,  4*lsq,  0,      0,    6*l,  2*lsq,  0,
    6*l,  0,    0,      4*lsq, -6*l,  0,    0,      2*lsq,
   -12,   0,    0,     -6*l,    12,   0,    0,     -6*l,
    0,   -12,   6*l,    0,      0,    12,   6*l,    0,
    0,   -6*l,  2*lsq,  0,      0,    6*l,  4*lsq,  0,
    6*l,  0,    0,      2*lsq, -6*l,  0,    0,      4*lsq;

    localK *= (eMod * momInert) / pow(l,3);


    // Construct the global mass- and gyro matrix (of size numDof x numDof)
    for (size_t i = a; i <= b; ++i) {
      for (size_t j = a; j <= b; ++j) {
        M(i, j) += localMLin(i - e*4, j - e*4);
        G(i, j) +=    localG(i - e*4, j - e*4);
        K(i, j) +=    localK(i - e*4, j - e*4);
      }
    }

    a += 4; b += 4;
  }
}


void ModalAnalysis::buildStateSpace(){

  std::vector<Triplet> aTripList, bTripList;

  for (size_t i = 0; i < numDof; i++) {
    for (size_t j = 0; j < numDof; j++) {
      if (M(i,j) != 0.0)  bTripList.push_back(Triplet(i, j, M(i,j)));
      if (G(i,j) != 0.0)  aTripList.push_back(Triplet(i, j, omega*G(i,j)));
      if (K(i,j) != 0.0) {
        aTripList.push_back(Triplet(i       , j+numDof, -K(i,j)));
        aTripList.push_back(Triplet(i+numDof, j       ,  K(i,j)));
        bTripList.push_back(Triplet(i+numDof, j+numDof,  K(i,j)));
      }
    }
  }

  // Populate sparse matrices
  A.setFromTriplets(aTripList.begin(), aTripList.end());
  B.setFromTriplets(bTripList.begin(), bTripList.end());
}


// Generalized EVP: A*phi_i = lambda_i*B*phi_i
void ModalAnalysis::solve()
{
  Eigen::VectorXcd eigVals;
  Eigen::MatrixXcd eigVects;
  Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges;

  // Create the A and B matrices in the EVP
  buildStateSpace();

  // Compute solution
  ges.compute(A, B);

  // Copy solution to tmp. variables
  eigVals  = ges.eigenvalues();
  eigVects = ges.eigenvectors();

  // Collect eigen- values and vectors in global "eigenSolution" container
  for (size_t i = 0; i < eigVals.size(); i++) {
    // Create tuple
    std::tuple<std::complex<double>,Eigen::VectorXcd> eigPair(eigVals[i],
                                                              eigVects.row(i));
    // and append it to container
    eigenSolution.push_back(eigPair);
  }

  // Sort the eigen solution
  std::sort(eigenSolution.begin(), eigenSolution.end(), util::eigSort);
}


void ModalAnalysis::printInfo()
{
  std::cout << "Number of elements: " << numEl << std::endl;
  std::cout << "Number of DOFs: " << numDof << std::endl;
}
