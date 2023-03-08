/*
  Spectrum.h

  Contains a class which allows for generation of arbitrary spectra

  S. Jones 08-03-23
*/

#ifndef SPECTRUM_H
#define SPECTRUM_H

namespace spec {
class Spectrum {
 private:
  bool normalOrdering;  // If true, normally ordered
  double mBeta;         // Effective neutrino mass [eV/c^2]
  double exposure;      // Exposure [atom years]
  double spectrumSize;  // Spectrum size [eV]
  double endpoint;      // Spectrum endpoint [eV]
  double background;    // Background rate [eV^-1 s^-1]
  double windowFrac;    // Fraction of events within window

  // Masses of the neutrino mass eigenstates
  double m1;
  double m2;
  double m3;

  double CalcMBetaFromStates(double m1, double m2, double m3, bool no);

  /// @brief Calculate differential decay rate for a given electron KE
  /// @param electronT Electron kinetic energy [eV]
  /// @return Returns differential decay rate [s^-1 eV^-1]
  double dGammadE(double electronT);

 public:
  Spectrum(bool NO, double nuMass, double exp, double specSize,
           double endE = 18575, double bkg = 1e-6);

  double GetSpecSize() { return spectrumSize; }
  double GetBkgRate() { return background; }
  double GetDecayFrac() { return windowFrac; }
  double GetMBeta() { return mBeta; }
};
}  // namespace spec

#endif