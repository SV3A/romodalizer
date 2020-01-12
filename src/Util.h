#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <complex>
#include "Eigen/Dense"

typedef std::tuple<std::complex<double>, Eigen::VectorXcd> EigTuple;

namespace util
{
// Function to evaluate two complex numbers, to be used by std::sort.  Returns
// 1 if a is "smaller" than b
bool complexSort(std::complex<double> &a, std::complex<double> &b);

bool eigSort(EigTuple &a, EigTuple &b);
}

#endif
