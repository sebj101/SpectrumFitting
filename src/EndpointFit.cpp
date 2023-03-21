/*
  EndpointFit.cpp

  Vary the range over which we fit the tritium endpoint and see what it does to
  our measure values
*/

#include <chrono>
#include <iostream>
#include <memory>
#include <vector>

#include "Spectrum.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "src/include/CommonFuncs.h"
#include "src/include/CommonNumbers.h"

using namespace spec;
using namespace std::chrono;

double dGammadEFit(double *x, double *par) {
  double E{x[0]};
  double E0{par[0]};
  double mBetaSq{par[1]};
  double A{par[2]};
  double b{par[3]};  // Background term
  return A * DifferentialTDecayRate(E0, E, sqrt(mBetaSq)) + b;
}

int main(int argc, char *argv[]) {
  TString outputFile{argv[1]};
  auto fout = std::make_unique<TFile>(outputFile, "RECREATE");

  // Create a spectrum
  const double nTestAtoms{1e19};
  const double runTime{1};     // years
  const double specSize{200};  // eV
  const double bkgRate{1e-6};  // eV^-1 s^-1

  // True endpoint parameters
  const double litQ{18575.72};  // eV
  const double litQUnc{70e-3};  // eV
  const int nSpectra{100};
  TH1D hE0True("hE0True", "True endpoint; E_{0, true} [keV]; N", 40,
               litQ - 3 * litQUnc, litQ + 3 * litQUnc);
  SetHistStyle(hE0True);
  TH1D hE0Reco("hE0Reco", "Reco endpoint; E_{0, reco} [keV]; N", 40,
               litQ - 3 * litQUnc, litQ + 3 * litQUnc);
  SetHistStyle(hE0Reco);
  TH1D hE0Res("hE0Res", "Endpoint; E_{0, reco} - E_{0, true} [eV]; N", 100,
              -litQUnc, litQUnc);
  SetHistStyle(hE0Res);
  for (int n{0}; n < nSpectra; n++) {
    Spectrum testSpec(runTime, nTestAtoms, specSize, bkgRate);
    const double trueE0{testSpec.GetQValue()};
    hE0True.Fill(trueE0);

    // Get the filled spectrum
    auto hSpecNorm{testSpec.GetSpectrumNorm()};
    hSpecNorm.SetTitle(Form("E_{0} = %.3f eV, m_{#beta} = %.3f eV",
                            testSpec.GetQValue(), testSpec.GetMBeta()));
    fout->cd();
    hSpecNorm.Write(Form("hSpecNorm_%d", n));

    // First of all, fit the background beyond the endpoint
    const double bkgFRLo{18580};
    const double bkgFRHi{18585};
    auto fBkgNorm = new TF1("fBkgNorm", "pol0", bkgFRLo, bkgFRHi);
    hSpecNorm.Fit(fBkgNorm, "qr");
    const double fitBkgNorm{fBkgNorm->GetParameter(0)};
    delete fBkgNorm;

    std::cout << "Fixing endpoint...\n";
    const double endpointFRLo{18375};
    const double endpointFRHi{18475};

    auto fEndpointNorm =
        new TF1("fEndpointNorm", dGammadEFit, endpointFRLo, endpointFRHi, 4);
    fEndpointNorm->SetParNames("E_{0}", "m_{#beta}^2", "A", "Bkg");
    fEndpointNorm->SetParameter(0, litQ);
    fEndpointNorm->FixParameter(2, 1.0);
    fEndpointNorm->FixParameter(1, 0);
    fEndpointNorm->FixParameter(3, fitBkgNorm);
    hSpecNorm.Fit(fEndpointNorm, "r");
    const double fittedE0{fEndpointNorm->GetParameter(0)};
    std::cout << "Fitted E0 = " << fittedE0
              << " eV\tReco E0 - True E0 = " << fittedE0 - trueE0 << " eV\n";

    hE0Reco.Fill(fittedE0);
    hE0Res.Fill(fittedE0 - trueE0);

    fout->cd();
    fEndpointNorm->Write(Form("fE0_%d", n));
    delete fEndpointNorm;
  }

  fout->cd();
  hE0True.Write();
  hE0Reco.Write();
  hE0Res.Write();

  fout->Close();
  return 0;
}