/*
  ToyEndpoint.cpp

  Makes a sample endpoint spectrum for a given neutrino mass
*/

#include "TBetaGenerator.hpp"

#include <span>

int main(int argc, char *argv[])
{
  double mBeta{0}; // Effective neutrino mass

  auto args = std::span(argv, size_t(argc));
  mBeta = std::stod(args[1]);

  // No steriles please
  const double mN{0};
  const double eta{0};

  return 0;
}