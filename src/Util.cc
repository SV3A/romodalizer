#include "Util.h"

bool util::complexSort(std::complex<double> &a, std::complex<double> &b) {
  // If real part equal sort by imag part, else sort by real part
  if (real(a) == real(b))
    return imag(a) < imag(b);
  return real(a) < real(b);
}


bool util::eigSort(EigTuple &a, EigTuple &b) {
  return abs(imag(std::get<0>( a ))) < abs(imag(std::get<0>( b )));
}


void util::printEigenvalues(std::vector<EigTuple>* const eigSol,
                            std::string const& output,
                            unsigned int incr)
{
  const char* formatSpec;
  std::vector<float> eigs;

  // Extract every incr'th eigenvalue
  for (std::vector<EigTuple>::iterator it = eigSol->begin();
      it != eigSol->end(); it += incr) {
    eigs.push_back(abs(imag(std::get<0>(*it))));

    // rad/s -> Hz
    if (output == "eigenfrequencies")
      eigs.back() *= 1.0/(2.0*M_PI);
  }

  if (output == "eigenvalues") {
    std::cout << "\nEigenvalues:" << std::endl;
    formatSpec = "%14.4f rad/s\n";

  } else if (output == "eigenfrequencies") {
    std::cout << "\nNatural Frequencies:" << std::endl;
    formatSpec = "%14.4f Hz\n";
  }

  for (auto const& eig:eigs) {
      printf(formatSpec, eig);
  }
}
