/*
  Spectrum.cxx
*/

#include "Spectrum.h"

#include <iostream>

#include "src/include/NuFitParams.h"

spec::Spectrum::Spectrum(bool NO, double nuMass, double time, double specSize,
                         double endE, double bkg)
    : normalOrdering(NO),
      mBeta(nuMass),
      runTime(time),
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
  }
}

double spec::Spectrum::CalcMBetaFromStates(double m1, double m2, double m3,
                                           bool no) {
  if (no) {
    return sqrt(Ue1SqNH * m1 * m1 + Ue2SqNH * m2 * m2 + Ue3SqNH * m3 * m3);
  } else {
    return sqrt(Ue1SqIH * m1 * m1 + Ue2SqIH * m2 * m2 + Ue3SqIH * m3 * m3);
  }
}