/*
  Spectrum.h

  Contains a class which allows for generation of arbitrary spectra

  S. Jones 08-03-23
*/

#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <memory>

#include "TH1.h"

namespace spec {
class Spectrum {
 private:
  bool normalOrdering;  // If true, normally ordered
  double mBeta;         // Effective neutrino mass [eV/c^2]
  double runningTime;   // Running time [years]
  double nAtoms;        // Number of atoms
  double spectrumSize;  // Spectrum size [eV]
  double endpoint;      // Spectrum endpoint [eV]
  double background;    // Background rate [eV^-1 s^-1]
  double windowFrac;    // Fraction of events within window

  // Masses of the neutrino mass eigenstates
  double m1;
  double m2;
  double m3;

  TH1D hSpec;

  /// @brief Calculate the effective neutrino mass from the mass eigenstates
  /// @param m1 Mass of m1 [eV]
  /// @param m2 Mass of m2 [eV]
  /// @param m3 Mass of m3 [eV]
  /// @param no Is the neutrino mass hierarchy normal?
  /// @return Effective neutrino mass [eV]
  double CalcMBetaFromStates(double m1, double m2, double m3, bool no);

  /// @brief Calculate differential decay rate for a given electron KE
  /// @param electronT Electron kinetic energy [eV]
  /// @return Returns differential decay rate [s^-1 eV^-1]
  double dGammadE(double electronT);

  /// @brief Integrates the differential decay spectrum up to the endpoint,
  /// using the rectangle method
  /// @param eMin Lower bound to integrate from
  /// @param eMax Upper bound to integrate to
  /// @param nPnts Number of bins to use for the integration
  /// @return Spectrum integral
  double SpectrumIntegral(double eMin, double eMax, int nPnts);

 public:
  Spectrum(bool NO, double nuMass, double time, double atoms, double specSize,
           double endE = 18575, double bkg = 1e-6);

  double GetSpecSize() { return spectrumSize; }
  double GetBkgRate() { return background; }
  double GetDecayFrac() { return windowFrac; }
  double GetMBeta() { return mBeta; }

  std::unique_ptr<TH1D> GetSpectrum() { return std::make_unique<TH1D>(hSpec); }
};
}  // namespace spec

#endif