#include "Util.h"

// Time points used by tic/toc
static HPTimeMark* start;
static HPTimeMark* end;


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


void util::tic(){
  start = new HPTimeMark;
  end   = new HPTimeMark;

  *start = std::chrono::high_resolution_clock::now();
}


void util::toc(){
  *end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = *end - *start;

  std::cout << "Execution time: " << elapsed.count() << " s" << std::endl;

  delete start;
  delete end;
}


util::Timer::Timer() {
  startMark = std::chrono::high_resolution_clock::now();
}


util::Timer::Timer(std::string name) {
  startMark = std::chrono::high_resolution_clock::now();
  timerName = name;
}


util::Timer::~Timer() {
  stop();
}


void util::Timer::stop()
{
  auto endMark = std::chrono::high_resolution_clock::now();

  // Cast to microseconds
  auto start_us = std::chrono::time_point_cast
    <std::chrono::microseconds>(startMark).time_since_epoch().count();

  auto end_us   = std::chrono::time_point_cast
    <std::chrono::microseconds>(endMark).time_since_epoch().count();

  auto elapsed = end_us - start_us;

  if (timerName.size() > 1)
    std::cout << "Execution time (" << timerName << "): ";
  else
    std::cout << "Execution time: ";

  if (elapsed > 1e6)      std::cout << elapsed*1e-6  << " s";
  else if (elapsed > 1e3) std::cout << elapsed*0.001 << " ms";
  else                    std::cout << elapsed       << " μs";

  std::cout << std::endl;
}
