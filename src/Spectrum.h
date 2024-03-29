/*
  Spectrum.h

  Contains a class which allows for generation of arbitrary spectra

  S. Jones 08-03-23
*/

#ifndef SPECTRUM_H
#define SPECTRUM_H

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
  double dE;            // Energy interval size [eV]

  // Masses of the neutrino mass eigenstates
  double m1;
  double m2;
  double m3;

  // Histogram containing the spectrum
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

  /// @brief Returns the minimum possible value of the effective neutrino mass
  /// @return The minimum effective neutrino mass in eV
  double GetMBetaMin();

  /// @brief Calculates the fraction of the total decays happening in window
  /// @return Fraction of decays in this bit of the beta decay spectrum
  double CalcWindowFrac();

  /// @brief Actually fills the fake spectrum according to the class parameters
  void FillSpectrum();

 public:
  /// @brief Constructor where neutrino mass and endpoint are specified
  /// @return Spectrum object
  Spectrum(bool NO, double nuMass, double time, double atoms, double specSize,
           double endE, double bkg, double deltaE = 0);

  /// @brief Constructor for random spectrum
  /// @return Spectrum object
  Spectrum(double time, double atoms, double specSize, double bkg,
           double deltaE = 0);

  // Member data getters
  double GetSpecSize() { return spectrumSize; }

  double GetBkgRate() { return background; }

  double GetDecayFrac() { return windowFrac; }

  double GetMBeta() { return mBeta; }

  double GetQValue() { return endpoint; }

  bool IsNormallyOrdered() { return normalOrdering; }

  TH1D GetSpectrum() { return hSpec; }

  /// @brief Return a normalised spectrum (units of eV^-1 s^-1 on y axis)
  /// @return The normalised histogram
  TH1D GetSpectrumNorm();

  /// @brief Calculates (analytically) the standard deviation on m_beta^2
  /// @return Standard deviation on m_beta^2 [eV^2]
  double GetSigmaMBetaSq();
};
}  // namespace spec

#endif