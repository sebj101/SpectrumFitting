/*
  ToyEndpoint.cpp

  Makes a sample endpoint spectrum for a given neutrino mass
*/

#include <iostream>
#include <memory>
#include <random>

#include "TBetaGenerator.hpp"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"

using std::cout;
using std::endl;
using std::make_unique;
using std::unique_ptr;

// build generation functor
// operator calls eactly one distribution function.
// Make different functor for alternative distribution.
class betaGenerator {
  // data parameters
 private:
  bool order_;
  double munu_;
  double mNu_;
  double mixing_;

 public:
  // fix parameter at construction
  betaGenerator(bool o, double ml, double mn, double eta)
      : order_(o), munu_(ml), mNu_(mn), mixing_(eta) {}

  double operator()(double x) {
    return TBeta::dGammadE(order_, munu_, mNu_, mixing_, x);
  }
};

double SpectrumIntegral(double eMin, double eMax, int nBins, bool order,
                        double mBeta) {
  double integral{0};
  const double rectangleWidth{(eMax - eMin) / double(nBins)};
  double eEval{eMin + rectangleWidth / 2};
  for (int n{0}; n < nBins; n++) {
    integral += rectangleWidth * TBeta::dGammadE(order, mBeta, 0, 0, eEval);
    eEval += rectangleWidth;
  }
  return integral;
}

// Main
int main(int argc, char *argv[]) {
  typedef std::piecewise_linear_distribution<double> pld;

  double mBeta{std::stod(argv[1])};  // Effective neutrino mass
  TString outputFile{argv[2]};
  auto fout = make_unique<TFile>(outputFile, "RECREATE");

  // No steriles please
  const double mN{0};
  const double eta{0};

  bool ordering{true};  // Normal ordering

  // Distribution parameters
  const int nBins{1000};
  const double minEnergy{0.01};                    // keV
  const double maxEnergy{TBeta::endAt(mBeta, 1)};  // Maximum energy
  // Generate the distribution
  pld spectrum(nBins, minEnergy, maxEnergy,
               betaGenerator(ordering, mBeta, mN, eta));

  // Need to calculate how many times to draw from spectrum for a hypothetical
  // exposure of tritium
  const double totalSpectrumIntegral{
      SpectrumIntegral(0, maxEnergy, 40000, ordering, mBeta)};
  const double windowSpectrumIntegral{SpectrumIntegral(
      maxEnergy - minEnergy, maxEnergy, 20000, ordering, mBeta)};

  cout << "Window fraction = " << windowSpectrumIntegral / totalSpectrumIntegral
       << endl;

  fout->Close();
  return 0;
}