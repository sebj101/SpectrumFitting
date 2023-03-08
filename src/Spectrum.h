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
  double eta;           // Fraction of events within window

  // Masses of the neutrino mass eigenstates
  double m1;
  double m2;
  double m3;

  double CalcMBetaFromStates(double m1, double m2, double m3, bool no);

 public:
  Spectrum(bool NO, double nuMass, double exp, double specSize,
           double endE = 18575, double bkg = 1e-6);

  double GetSpecSize() { return spectrumSize; }
  double GetBkgRate() { return background; }
  double GetDecayFrac() { return eta; }
  double GetMBeta() { return mBeta; }
};
}  // namespace spec

#endif