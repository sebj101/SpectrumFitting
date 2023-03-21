/*
  Spectrum.cxx
*/

#include "Spectrum.h"

#include <chrono>
#include <iostream>
#include <random>
#include <vector>

#include "src/include/CommonNumbers.h"
#include "src/include/FundamentalConstants.h"
#include "src/include/NuFitParams.h"

spec::Spectrum::Spectrum(bool NO, double nuMass, double time, double atoms,
                         double specSize, double endE, double bkg,
                         double deltaE)
    : normalOrdering(NO),
      mBeta(nuMass),
      runningTime(time),
      nAtoms(atoms),
      spectrumSize(specSize),
      endpoint(endE),
      background(bkg),
      dE(deltaE) {
  // First calculate the minimum value of the neutrino mass
  double mBetaMin{GetMBetaMin()};

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
  }

  windowFrac = CalcWindowFrac();

  // Actually fill the spectrum
  FillSpectrum();
}

spec::Spectrum::Spectrum(double time, double atoms, double specSize, double bkg,
                         double deltaE)
    : runningTime(time),
      nAtoms(atoms),
      spectrumSize(specSize),
      background(bkg),
      dE(deltaE) {
  // Seed generator
  long int seed{std::chrono::system_clock::now().time_since_epoch().count()};
  std::mt19937 rng(seed);
  // Generate the random mass hierarchy
  std::uniform_real_distribution<double> dist(0, 1);
  if (dist(rng) > 0.5) {
    normalOrdering = true;
  } else {
    normalOrdering = false;
  }

  // Generate the lightest neutrino mass
  const double mLeastMax{0.8};  // eV
  const double mLeastMin{0};    // eV
  // Flat throw between these ranges
  std::uniform_real_distribution<double> mLeastDist(mLeastMin, mLeastMax);
  if (normalOrdering) {
    m1 = mLeastDist(rng);
    m2 = sqrt(m1 * m1 + kNuFitDmsq21NH);
    m3 = sqrt(m1 * m1 + kNuFitDmsq31NH);
  } else {
    m3 = mLeastDist(rng);
    m2 = sqrt(m3 * m3 - kNuFitDmsq32IH);
    m1 = sqrt(m2 * m2 + kNuFitDmsq21NH);
  }
  // Calculate effective neutrino mass
  mBeta = CalcMBetaFromStates(m1, m2, m3, normalOrdering);

  // Use a gaussian throw for the endpoint
  // Q value and uncertainty from
  // Myers, Wagner, Kracke, & Wesson. PRL 114, 013003 (2015)
  const double litQ{18575.72};  // eV
  const double litQUnc{0.07};   // eV
  std::normal_distribution<double> QDist(litQ, litQUnc);
  endpoint = QDist(rng);

  // Calculate fraction of events in window
  windowFrac = CalcWindowFrac();

  // Fill the spectrum
  FillSpectrum();
}

double spec::Spectrum::CalcWindowFrac() {
  const double totalSpecInt{SpectrumIntegral(0, endpoint, 400000)};
  const double windowSpecInt{
      SpectrumIntegral(endpoint - spectrumSize, endpoint, 20000)};
  return windowSpecInt / totalSpecInt;
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

double spec::Spectrum::SpectrumIntegral(double eMin, double eMax, int nBins) {
  double integral{0};
  const double rectangleWidth{(eMax - eMin) / double(nBins)};
  const double eEvalInit{eMin + rectangleWidth / 2};
  for (int n{0}; n < nBins; n++) {
    double eEval{eEvalInit + double(n) * rectangleWidth};
    integral += rectangleWidth * dGammadE(eEval);
  }
  return integral;
}

double spec::Spectrum::GetMBetaMin() {
  if (normalOrdering) {
    const double m1Min{0};
    const double m2Min{sqrt(m1Min * m1Min + kNuFitDmsq21NH)};
    const double m3Min{sqrt(m1Min * m1Min + kNuFitDmsq31NH)};
    return CalcMBetaFromStates(m1Min, m2Min, m3Min, normalOrdering);
  } else {
    const double m3Min{0};
    const double m2Min{sqrt(m3Min * m3Min - kNuFitDmsq32IH)};
    const double m1Min{sqrt(m2Min * m2Min + kNuFitDmsq21NH)};
    return CalcMBetaFromStates(m1Min, m2Min, m3Min, normalOrdering);
  }
}

void spec::Spectrum::FillSpectrum() {
  const double lasteVFrac{2.9e-13};
  // Calculate rate in last eV in absence of mass
  const double r{nAtoms * lasteVFrac / T_MEAN_LIFETIME};
  // Calculate optimum energy window
  const double deltaEOpt{sqrt(background / r)};
  // If left as default, use this for energy binning
  if (dE == 0) dE = deltaEOpt;

  // Number of bins
  const double distanceBeyondE0{10};
  int nDistBins{int(std::round((spectrumSize + distanceBeyondE0) / dE))};

  hSpec = TH1D("", "; Electron energy [eV]; N_{electrons}", nDistBins,
               endpoint - spectrumSize, endpoint + distanceBeyondE0);
  hSpec.SetLineColor(kBlack);
  hSpec.GetXaxis()->SetTitleSize(.05);
  hSpec.GetYaxis()->SetTitleSize(.05);
  hSpec.GetXaxis()->SetLabelSize(.05);
  hSpec.GetYaxis()->SetLabelSize(.05);

  // Fill the histogram
  // For each bin, calculate the mean number of decays and then use a Poisson
  // distribution to get the acutal number of events
  long int seed{std::chrono::system_clock::now().time_since_epoch().count()};
  std::mt19937 rng(seed);

  // Do the background events at the same time
  // Poisson dist with mean of the calculated number of events per bin
  // Calculate the number of background throws
  const double avgBkgEvents{background * hSpec.GetBinWidth(1) * ONE_YEAR *
                            runningTime};
  std::poisson_distribution<int> pBkg(avgBkgEvents);

  for (int iBin{1}; iBin <= hSpec.GetNbinsX(); iBin++) {
    // Integrate over bin
    double binDecayRate{SpectrumIntegral(
        hSpec.GetBinLowEdge(iBin),
        hSpec.GetBinLowEdge(iBin) + hSpec.GetBinWidth(iBin), 10)};
    double meanDecays{binDecayRate * nAtoms * runningTime * ONE_YEAR};
    std::poisson_distribution<long> pSignal(meanDecays);
    int nBkg{pBkg(rng)};
    long nSignal{pSignal(rng)};
    double totalBinEvents{double(nSignal) + double(nBkg)};
    hSpec.SetBinContent(iBin, totalBinEvents);
  }
}

double spec::Spectrum::GetSigmaMBetaSq() {
  const double lasteVFrac{2.9e-13};
  const double t{ONE_YEAR * runningTime};  // seconds
  // Calculate rate in last eV in absence of mass
  const double r{nAtoms * lasteVFrac / T_MEAN_LIFETIME};
  // Calculate optimum energy window
  double sigma{sqrt(r * t * dE + background * t / dE)};
  return 2 / (3 * r * t) * sigma;
}

TH1D spec::Spectrum::GetSpectrumNorm() {
  TH1D hSpecNorm{hSpec};  // Copy the original spectrum
  hSpecNorm.GetYaxis()->SetTitle("Events [eV^{-1} s^{-1} atom^{-1}]");
  // Scale by the bin width
  hSpecNorm.Scale(1, "width");
  // Divide by the number of seconds and the number of atoms
  hSpecNorm.Scale(1 / (runningTime * ONE_YEAR * nAtoms));
  return hSpecNorm;
}