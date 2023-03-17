/*
  FitRangeTest.cpp

  Test the effect of varying the fit range
*/

#include <iostream>
#include <memory>
#include <random>

#include "Spectrum.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "src/include/FundamentalConstants.h"

using namespace spec;

const double chosenEndpoint{18575};  // eV

double dGammadE(double A, double E0, double E, double mBeta) {
  double nuE{E0 - E};
  if (nuE - mBeta < 0) {
    return 0;
  } else {
    double beta{sqrt(1.0 - pow(EMASS_EV / (E + EMASS_EV), 2))};
    double p{sqrt(E * E + 2 * EMASS_EV * E)};
    double eta{ALPHA * 2 / beta};
    double fermiFunction{2 * M_PI * eta / (1 - exp(-2 * M_PI * eta))};
    double rate{G_F * G_F * pow(0.97425, 2) * fermiFunction *
                (1 + 3.0 * pow(-1.2646, 2)) * p * (E + EMASS_EV) /
                (2 * pow(M_PI, 3))};
    return rate * A * nuE * sqrt(nuE * nuE - mBeta * mBeta) / 6.58e-16;
  }
}

double dGammadEFit(double *x, double *par) {
  double E{x[0]};
  double E0{par[0]};
  double mBetaSq{par[1]};
  double A{par[2]};
  double b{par[3]};  // Background term
  return dGammadE(A, E0, E, sqrt(mBetaSq)) + b;
}

int main(int argc, char *argv[]) {
  TString outputFile{argv[1]};
  auto fout = std::make_unique<TFile>(outputFile, "RECREATE");

  // Create a spectrum
  const double nTestAtoms{1e19};
  const double runTime{1};     // years
  const double specSize{200};  // eV
  const double bkgRate{1e-6};  // eV^-1 s^-1
  // Spectrum testSpec(runTime, nTestAtoms, specSize, bkgRate);

  // Generate random neutrino mass
  std::uniform_real_distribution<double> mBetaDist(0.05, 0.8);
  // Seed generator
  long int seed{std::chrono::system_clock::now().time_since_epoch().count()};
  std::mt19937 rng(seed);
  Spectrum testSpec(true, mBetaDist(rng), runTime, nTestAtoms, specSize,
                    chosenEndpoint, bkgRate);

  // Get the filled spectrum
  auto hSpec{testSpec.GetSpectrum()};
  fout->cd();
  hSpec.Write("hSpec");

  const double endFitRange{testSpec.GetQValue() + 10};

  std::cout << "Analytic sigma_{m_beta^2} = " << testSpec.GetSigmaMBetaSq()
            << " eV^2\n"
            << std::endl;

  // First of all, fit the background beyond the endpoint
  const double bkgFRLo{18580};
  const double bkgFRHi{18585};
  auto fBkg = new TF1("fBkg", "pol0", bkgFRLo, bkgFRHi);
  hSpec.Fit(fBkg, "r");
  const double fitBkg{fBkg->GetParameter(0)};
  delete fBkg;

  std::cout << "Fixing endpoint...\n";
  const double endpointFRLo{18375};
  const double endpointFRHi{18475};
  auto fEndpoint =
      new TF1("fEndpoint", dGammadEFit, endpointFRLo, endpointFRHi, 4);
  fEndpoint->SetParNames("E_{0}", "m_{#beta}^2", "A", "Bkg");
  fEndpoint->SetParameter(0, 18575.72);
  fEndpoint->SetParameter(2, 1.375e25);
  fEndpoint->FixParameter(1, 0);
  fEndpoint->FixParameter(3, fitBkg);
  hSpec.Fit(fEndpoint, "r");
  const double fittedEndpoint{fEndpoint->GetParameter(0)};
  std::cout << "Fitted endpoint = " << fittedEndpoint << " eV\n";
  delete fEndpoint;

  std::cout << "\n========================================\n";
  std::cout << "Full fit\n";
  auto fFull =
      new TF1("fFull", dGammadEFit, hSpec.GetBinLowEdge(1), endFitRange, 4);
  fFull->SetParName(0, "E_{0}");
  fFull->SetParName(1, "m_{#beta}^2");
  fFull->SetParName(2, "A");
  fFull->SetParName(3, "Bkg");
  fFull->FixParameter(0, fittedEndpoint);
  fFull->SetParameter(1, pow(0.2, 2));
  fFull->SetParameter(2, 1.375e25);
  fFull->FixParameter(3, fitBkg);
  fFull->SetNpx(2 * hSpec.GetNbinsX());
  hSpec.Fit(fFull, "");
  std::cout << "m_beta^2: True, reco, difference = "
            << pow(testSpec.GetMBeta(), 2) << " eV^2,\t("
            << fFull->GetParameter(1) << " +- " << fFull->GetParError(1)
            << ") eV^2,\t"
            << pow(testSpec.GetMBeta(), 2) - fFull->GetParameter(1)
            << " eV^2\n";
  std::cout << "Reduced chisq = " << fFull->GetChisquare() / fFull->GetNDF()
            << std::endl;

  fout->cd();
  fFull->Write();

  fout->Close();
  return 0;
}