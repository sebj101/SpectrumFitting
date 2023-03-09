/*
  ToyEndpoint.cpp

  Makes a sample endpoint spectrum for a given neutrino mass
*/

#include <iostream>
#include <memory>
#include <random>

#include "Spectrum.h"
#include "TAxis.h"
#include "TBetaGenerator.hpp"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TMath.h"
#include "TString.h"

using namespace spec;

using std::cout;
using std::endl;
using std::make_unique;
using std::unique_ptr;

const double ME{9.1093837015e-31};
const double ME_EV{ME * TMath::C() * TMath::C() / TMath::Qe()};
const double G_F{1.1663787e-23};
const double ALPHA{7.2973525698e-3};

// PMNS matrix elements
const double kNuFitTh12NH{33.44 * TMath::Pi() / 180.0};
const double kNuFitTh13NH{8.57 * TMath::Pi() / 180.0};
const double Ue1Sq{pow(cos(kNuFitTh12NH) * cos(kNuFitTh13NH), 2)};
const double Ue2Sq{pow(sin(kNuFitTh12NH) * cos(kNuFitTh13NH), 2)};
const double Ue3Sq{pow(sin(kNuFitTh13NH), 2)};

// Mass splittings
const double kNuFitDmsq21{7.42e-5};
const double kNuFitDmsq32{23.2e-4};

double M1FromMBeta(double mBeta) {
  double m1sq{
      (mBeta * mBeta - Ue2Sq * kNuFitDmsq21 -
       Ue3Sq * (kNuFitDmsq21 + kNuFitDmsq32) / (Ue1Sq + Ue2Sq + Ue3Sq))};
  return sqrt(m1sq);
}

double SebDecayRate(double electronT, double m1, double m2, double m3,
                    double endpointE = 18575) {
  double beta{sqrt(1.0 - pow(ME_EV / (electronT + ME_EV), 2))};
  double p{sqrt(electronT * electronT + 2 * ME_EV * electronT)};
  double eta{ALPHA * 2 / beta};
  double fermiFunction{2 * TMath::Pi() * eta /
                       (1 - exp(-2 * TMath::Pi() * eta))};
  double nuE{endpointE - electronT};
  double rate{G_F * G_F * pow(0.97425, 2) * fermiFunction *
              (1 + 3.0 * pow(-1.2646, 2)) * p * (electronT + ME_EV) /
              (2 * pow(TMath::Pi(), 3))};

  double neutrinoPhaseSpc1 =
      (m1 <= nuE) ? Ue1Sq * nuE * sqrt(nuE * nuE - m1 * m1) : 0.0;
  double neutrinoPhaseSpc2 =
      (m2 <= nuE) ? Ue2Sq * nuE * sqrt(nuE * nuE - m2 * m2) : 0.0;
  double neutrinoPhaseSpc3 =
      (m3 <= nuE) ? Ue3Sq * nuE * sqrt(nuE * nuE - m3 * m3) : 0.0;
  double neutrinoPhaseSpc =
      neutrinoPhaseSpc1 + neutrinoPhaseSpc2 + neutrinoPhaseSpc3;
  rate *= neutrinoPhaseSpc;
  rate *= 1.0 / 6.58e-16;  // Account for natural units

  return rate;
}

double SebDecayRate(double electronT, double mBeta, double endpointE = 18575) {
  if (mBeta == 0) {
    return SebDecayRate(electronT, 0, 0, 0, endpointE);
  } else {
    double m1{M1FromMBeta(mBeta)};
    double m2{sqrt(kNuFitDmsq21 + m1 * m1)};
    double m3{sqrt(kNuFitDmsq32 + m2 * m2)};
    return SebDecayRate(electronT, m1, m2, m3, endpointE);
  }
}

double SpectrumIntegral(double eMin, double eMax, int nBins, double mBeta) {
  double integral{0};
  const double rectangleWidth{(eMax - eMin) / double(nBins)};
  double eEval{eMin + rectangleWidth / 2};
  for (int n{0}; n < nBins; n++) {
    integral += rectangleWidth * SebDecayRate(eEval, mBeta);
    eEval += rectangleWidth;
  }
  return integral;
}

unique_ptr<TGraph> MakeSpectrumGraph(double endpointDist, double mBeta,
                                     int nPnts = 500) {
  auto gr = make_unique<TGraph>();
  gr->GetXaxis()->SetTitleSize(.05);
  gr->GetYaxis()->SetTitleSize(.05);
  gr->GetXaxis()->SetLabelSize(.05);
  gr->GetYaxis()->SetLabelSize(.05);
  gr->SetLineWidth(3);
  gr->GetXaxis()->SetTitle("E [keV]");
  gr->GetYaxis()->SetTitle("d#Gamma / dE");

  // Assume NH
  const double eMax{18575 + 2};
  const double eMin{18575 - endpointDist};
  cout << "Emin, Emax = " << eMin << ", " << eMax << endl;
  for (int iE{0}; iE < nPnts; iE++) {
    double E{eMin + (eMax - eMin) * double(iE) / double(nPnts - 1)};
    gr->SetPoint(iE, E * 1e-3, SebDecayRate(E, mBeta));
  }
  return gr;
}

// Main
int main(int argc, char *argv[]) {
  cout.precision(10);

  double mBeta{std::stod(argv[1])};  // Effective neutrino mass
  TString outputFile{argv[2]};
  auto fout = make_unique<TFile>(outputFile, "RECREATE");

  // Need to calculate how many times to draw from spectrum for a hypothetical
  // exposure of tritium
  const double chosenEndpoint{18575};
  const double totalSpectrumIntegralMassless{
      SpectrumIntegral(0, chosenEndpoint, 400000, 0)};
  const double spectrumIntegralLasteV{
      SpectrumIntegral(chosenEndpoint - 1, chosenEndpoint, 20000, 0)};
  const double spectrumIntegralLast100eV{
      SpectrumIntegral(chosenEndpoint - 100, chosenEndpoint, 20000, 0)};
  const double lasteVFraction{spectrumIntegralLasteV /
                              totalSpectrumIntegralMassless};
  const double last100eVFraction{spectrumIntegralLast100eV /
                                 totalSpectrumIntegralMassless};
  cout << "Last eV fraction = " << lasteVFraction << endl;
  cout << "Last 100 eV fraction = " << last100eVFraction << endl;

  auto gr1000meV = MakeSpectrumGraph(10, 1, 1000);
  gr1000meV->SetLineColor(kBlack);
  auto gr500meV = MakeSpectrumGraph(10, 0.5, 1000);
  gr500meV->SetLineColor(kRed);
  auto gr100meV = MakeSpectrumGraph(10, 0.1, 1000);
  gr100meV->SetLineColor(kGreen + 1);
  auto gr0meV = MakeSpectrumGraph(10, 0, 1000);
  gr0meV->SetLineColor(kBlue);

  auto grFullSpec = MakeSpectrumGraph(18570, 1, 20000);
  grFullSpec->SetLineColor(kBlack);

  fout->cd();
  grFullSpec->Write("grFullSpec");
  gr1000meV->Write("gr1000meV");
  gr500meV->Write("gr500meV");
  gr100meV->Write("gr100meV");
  gr0meV->Write("gr0meV");

  const double nAtoms{1e20};
  const double tauMean{560975924};
  const double lasteVRate{nAtoms * lasteVFraction / tauMean};
  const double last100eVRate{nAtoms * last100eVFraction / tauMean};

  cout << "Rate (1 eV, 100 eV) = " << lasteVRate << " s^-1, " << last100eVRate
       << " s^-1\n";

  const double bkgRate{1e-6};
  const double deltaEOpt{sqrt(bkgRate / lasteVRate)};
  const double year{31536000};
  const double requiredThrows100eV{std::round(year * last100eVRate)};

  int nDistBins{int(std::round((100 + 5) / deltaEOpt))};

  cout << "Filling " << nDistBins << " bins with " << requiredThrows100eV
       << " throws.\n";

  /////////////// Spectrum testing /////////////////////
  const double nTestAtoms{1e20};
  const double times{1};           // years
  const double specSize{100};      // eV
  const double nuMassLight{0.05};  // eV
  const double nuMassHeavy{1};     // eV
  const double nuMassMassless{0};  // eV

  Spectrum testSpecHeavy(true, nuMassHeavy, times, nTestAtoms, specSize);
  auto specHistHeavy = testSpecHeavy.GetSpectrum();
  Spectrum testSpecLight(true, nuMassLight, times, nTestAtoms, specSize);
  auto specHistLight = testSpecLight.GetSpectrum();
  Spectrum testSpecMassless(true, nuMassMassless, times, nTestAtoms, specSize);
  auto specHistMassless = testSpecMassless.GetSpectrum();

  fout->cd();

  specHistHeavy.Write("specHistHeavy");
  specHistLight.Write("specHistLight");
  specHistMassless.Write("specHistMassless");

  fout->Close();
  return 0;
}