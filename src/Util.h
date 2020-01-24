#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <vector>
#include <complex>
#include <chrono>
#include "Eigen/Dense"
#include "Eigen/Sparse"

typedef std::tuple<std::complex<double>, Eigen::VectorXcd> EigTuple;
typedef std::chrono::time_point<std::chrono::high_resolution_clock> HPTimeMark;

namespace util
{
  // Function to evaluate two complex numbers, to be used by std::sort.  It 
  // returns 1 if a is "smaller" than b
  bool complexSort(std::complex<double> &a, std::complex<double> &b);

  bool eigSort(EigTuple &a, EigTuple &b);

  // Print eigenvalues natural frequencies
  void printEigenvalues(std::vector<EigTuple>* const eigSol,
      std::string const& output,
      unsigned int incr = 1);

  // Timing functions
  void tic();
  void toc();

  // Class for scope timer
  class Timer
  {
    private:
      HPTimeMark startMark;
      std::string timerName;
      void stop();
    public:
      Timer();
      Timer(std::string name);
      ~Timer();
  };

  // Calculate condition number of matrix
  template <typename T>
  double cond(T& A) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A);

    return svd.singularValues()(0) /
           svd.singularValues()(svd.singularValues().size()-1);
  }
}

#endif
