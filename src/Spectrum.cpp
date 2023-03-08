/*
  Spectrum.cxx
*/

#include "Spectrum.h"

#include <iostream>

#include "src/include/FundamentalConstants.h"
#include "src/include/NuFitParams.h"

spec::Spectrum::Spectrum(bool NO, double nuMass, double exp, double specSize,
                         double endE, double bkg)
    : normalOrdering(NO),
      mBeta(nuMass),
      exposure(exp),
      spectrumSize(specSize),
      endpoint(endE),
      background(bkg) {
  // First calculate the masses of the neutrino eigenstates
  double mBetaMin{0};
  if (normalOrdering) {
    const double m1Min{0};
    const double m2Min{sqrt(m1Min * m1Min + kNuFitDmsq21NH)};
    const double m3Min{sqrt(m1Min * m1Min + kNuFitDmsq31NH)};
    mBetaMin = CalcMBetaFromStates(m1Min, m2Min, m3Min, true);
  } else {
    const double m3Min{0};
    const double m2Min{sqrt(m3Min * m3Min - kNuFitDmsq32IH)};
    const double m1Min{sqrt(m2Min * m2Min + kNuFitDmsq21NH)};
    mBetaMin = CalcMBetaFromStates(m1Min, m2Min, m3Min, false);
  }

  // Check that the neutrino mass is not too low
  // If it is, set neutrino to be massless
  if (mBeta < mBetaMin) {
    m1 = 0;
    m2 = 0;
    m3 = 0;
  } else {
    // Calculate the masses of the eigenstates
    if (normalOrdering) {
      m1 = sqrt(mBeta * mBeta - Ue2SqNH * kNuFitDmsq21NH -
                Ue3SqNH * (kNuFitDmsq21NH + kNuFitDmsq32NH));
      m1 /= sqrt(Ue1SqNH + Ue2SqNH + Ue3SqNH);
      m2 = sqrt(m1 * m1 + kNuFitDmsq21NH);
      m3 = sqrt(m1 * m1 + kNuFitDmsq31NH);
    } else {
      m3 = sqrt(mBeta * mBeta + kNuFitDmsq32IH * (Ue2SqIH + Ue1SqIH) -
                Ue1SqIH * kNuFitDmsq21NH);
      m3 /= sqrt(Ue1SqIH + Ue2SqIH + Ue3SqIH);
      m2 = sqrt(m3 * m3 - kNuFitDmsq32IH);
      m1 = sqrt(m2 * m2 + kNuFitDmsq21NH);
    }

    // Some output for checking we are doing this correctly
    std::cout << "Inputted m_beta = " << mBeta << " eV\t Calculated m_beta = "
              << CalcMBetaFromStates(m1, m2, m3, normalOrdering) << " eV\n";

    // Now time to calculate the fraction of events in our window
    const double totalSpecInt{SpectrumIntegral(0, 40000)};
    const double windowSpecInt{
        SpectrumIntegral(endpoint - spectrumSize, 20000)};
    windowFrac = windowSpecInt / totalSpecInt;
  }
}

double spec::Spectrum::dGammadE(double electronT) {
  double beta{sqrt(1.0 - pow(EMASS_EV / (electronT + EMASS_EV), 2))};
  double p{sqrt(electronT * electronT + 2 * EMASS_EV * electronT)};
  double eta{ALPHA * 2 / beta};
  double fermiFunction{2 * M_PI * eta / (1 - exp(-2 * M_PI * eta))};
  double nuE{endpoint - electronT};
  double rate{G_F * G_F * pow(0.97425, 2) * fermiFunction *
              (1 + 3.0 * pow(-1.2646, 2)) * p * (electronT + EMASS_EV) /
              (2 * pow(M_PI, 3))};

  double neutrinoPhaseSpc{0};
  if (normalOrdering) {
    double neutrinoPhaseSpc1 =
        (m1 <= nuE) ? Ue1SqNH * nuE * sqrt(nuE * nuE - m1 * m1) : 0.0;
    double neutrinoPhaseSpc2 =
        (m2 <= nuE) ? Ue2SqNH * nuE * sqrt(nuE * nuE - m2 * m2) : 0.0;
    double neutrinoPhaseSpc3 =
        (m3 <= nuE) ? Ue3SqNH * nuE * sqrt(nuE * nuE - m3 * m3) : 0.0;
    neutrinoPhaseSpc =
        neutrinoPhaseSpc1 + neutrinoPhaseSpc2 + neutrinoPhaseSpc3;
  } else {
    double neutrinoPhaseSpc1 =
        (m1 <= nuE) ? Ue1SqIH * nuE * sqrt(nuE * nuE - m1 * m1) : 0.0;
    double neutrinoPhaseSpc2 =
        (m2 <= nuE) ? Ue2SqIH * nuE * sqrt(nuE * nuE - m2 * m2) : 0.0;
    double neutrinoPhaseSpc3 =
        (m3 <= nuE) ? Ue3SqIH * nuE * sqrt(nuE * nuE - m3 * m3) : 0.0;
    neutrinoPhaseSpc =
        neutrinoPhaseSpc1 + neutrinoPhaseSpc2 + neutrinoPhaseSpc3;
  }
  rate *= neutrinoPhaseSpc;
  rate *= 1.0 / 6.58e-16;  // Account for natural units
  return rate;
}

double spec::Spectrum::CalcMBetaFromStates(double m1, double m2, double m3,
                                           bool no) {
  if (no) {
    return sqrt(Ue1SqNH * m1 * m1 + Ue2SqNH * m2 * m2 + Ue3SqNH * m3 * m3);
  } else {
    return sqrt(Ue1SqIH * m1 * m1 + Ue2SqIH * m2 * m2 + Ue3SqIH * m3 * m3);
  }
}

double spec::Spectrum::SpectrumIntegral(double eMin, int nBins) {
  double integral{0};
  double rectangleWidth{(endpoint - eMin) / double(nBins)};
  double eEval{eMin + rectangleWidth / 2};
  for (int n{0}; n < nBins; n++) {
    integral += rectangleWidth * dGammadE(eEval);
    eEval += rectangleWidth;
  }
  return integral;
}