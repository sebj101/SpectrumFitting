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
#include "TGraph.h"
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

  // True endpoint parameters
  const double litQ{18575.72};  // eV
  const double litQUnc{70e-3};  // eV

  // Vector of fit ranges
  std::vector<double> fitRanges{200, 100, 50, 20, 10, 5, 2, 1};
  std::vector<TH1D> E0ResVec;
  std::vector<TH1D> E0RecoVec;
  std::vector<TH1D> SigmaE0Vec;
  for (auto &fr : fitRanges) {
    TH1D hE0Res(
        Form("hE0Res_%.0f", fr),
        Form("Fit range = %.0f eV; E_{0, reco} - E_{0, true} [eV]; N", fr), 100,
        -0.005, 0.005);
    SetHistStyle(hE0Res);
    E0ResVec.push_back(hE0Res);

    TH1D hE0Reco(Form("hE0Reco_%.0f", fr),
                 Form("Fit range = %.0f eV; E_{0, reco} [keV]; N", fr), 40,
                 litQ - 3 * litQUnc, litQ + 3 * litQUnc);
    SetHistStyle(hE0Reco);
    E0RecoVec.push_back(hE0Reco);

    TH1D hSigmaE0(Form("hSigmaE0_%.0f", fr),
                  Form("Fit range = %.0f eV; #sigma_{E_{0}} [eV]; N", fr), 500,
                  0, 5e-4);
    SetHistStyle(hSigmaE0);
    SigmaE0Vec.push_back(hSigmaE0);
  }

  // Spectrum parameters
  const double nTestAtoms{1e19};
  const double runTime{1};     // years
  const double specSize{300};  // eV
  const double bkgRate{1e-6};  // eV^-1 s^-1

  const int nSpectra{500};
  TH1D hE0True("hE0True", "True endpoint; E_{0, true} [keV]; N", 40,
               litQ - 3 * litQUnc, litQ + 3 * litQUnc);
  SetHistStyle(hE0True);

  for (int n{0}; n < nSpectra; n++) {
    Spectrum testSpec(runTime, nTestAtoms, specSize, bkgRate);
    const double trueE0{testSpec.GetQValue()};
    hE0True.Fill(trueE0 / 1e3);
    std::cout << "Spectrum " << n << std::endl;

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

    for (unsigned int iR{0}; iR < fitRanges.size(); iR++) {
      const double fr{fitRanges.at(iR)};
      std::cout << "Range = " << fr << " eV\n";
      const double endpointFRLo{hSpecNorm.GetBinLowEdge(1)};  // eV
      const double endpointFRHi{endpointFRLo + fr};           // eV

      auto fEndpointNorm =
          new TF1("fEndpointNorm", dGammadEFit, endpointFRLo, endpointFRHi, 4);
      fEndpointNorm->SetParNames("E_{0}", "m_{#beta}^2", "A", "Bkg");
      fEndpointNorm->SetParameter(0, litQ);
      fEndpointNorm->FixParameter(2, 1.0);
      fEndpointNorm->FixParameter(1, 0);
      fEndpointNorm->FixParameter(3, fitBkgNorm);
      hSpecNorm.Fit(fEndpointNorm, "qr");
      const double fittedE0{fEndpointNorm->GetParameter(0)};
      std::cout << "Fitted E0 = " << fittedE0
                << " eV\tReco E0 - True E0 = " << fittedE0 - trueE0 << " eV\n";
      std::cout << "Reduced chisq = "
                << fEndpointNorm->GetChisquare() / fEndpointNorm->GetNDF()
                << std::endl;
      E0RecoVec.at(iR).Fill(fittedE0 / 1e3);
      E0ResVec.at(iR).Fill(fittedE0 - trueE0);
      SigmaE0Vec.at(iR).Fill(fEndpointNorm->GetParError(0));

      // fout->cd();
      // fEndpointNorm->Write(Form("fE0_%d", n));
      delete fEndpointNorm;
    }
    std::cout << "\n";
  }

  fout->cd();
  hE0True.Write();

  TGraph grResFR;
  SetGraphStyle(grResFR);
  grResFR.SetTitle(
      "Variation with fit range; Fit range [eV]; E_{0} resolution [eV]");
  TGraph grBiasFR;
  SetGraphStyle(grBiasFR);
  grBiasFR.SetTitle(
      "Variation with fit range; Fit range [eV]; E_{0} bias [eV]");

  for (unsigned int iR{0}; iR < fitRanges.size(); iR++) {
    auto fGaus = new TF1("fGaus", "gaus");
    E0ResVec.at(iR).Fit(fGaus);
    grResFR.SetPoint(grResFR.GetN(), fitRanges.at(iR), fGaus->GetParameter(2));
    grBiasFR.SetPoint(grBiasFR.GetN(), fitRanges.at(iR),
                      fGaus->GetParameter(1));

    E0RecoVec.at(iR).Write();
    E0ResVec.at(iR).Write();
    SigmaE0Vec.at(iR).Write();

    delete fGaus;
  }
  fout->cd();
  grResFR.Write("grResFR");
  grBiasFR.Write("grBiasFR");

  fout->Close();
  return 0;
}