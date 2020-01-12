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
