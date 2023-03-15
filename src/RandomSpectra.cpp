/*
  RandomSpectra.cpp

  Test the randomly generated endpoint spectra
*/

#include <memory>

#include "Spectrum.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"

using std::make_unique;
using std::unique_ptr;

using namespace spec;

int main(int argc, char *argv[]) {
  TString outputFile{argv[1]};
  auto fout = make_unique<TFile>(outputFile, "RECREATE");

  const int nSpectra{1000};
  const double nTestAtoms{5e19};
  const double runTime{1};     // years
  const double specSize{25};   // eV
  const double bkgRate{1e-6};  // eV^-1 s^-1

  TH1D hEndpoint("hEndpoint", "Endpoint; Endpoint [eV]; N_{spectra}", 100,
                 18575.72 - 0.21, 18575.72 - 0.21);
  TH1D hMBeta("hMBeta", "m_{#beta}; m_{#beta} [eV]; N_{spectra}", 100, 0, 0.3);
  TH1D hNO("hNO", "Ordering; Normally ordered; N_{#spectra}", 2, -0.5, 1.5);
  for (int n{0}; n < nSpectra; n++) {
    Spectrum testSpec(runTime, nTestAtoms, specSize, bkgRate);
    hEndpoint.Fill(testSpec.GetQValue());
    hMBeta.Fill(testSpec.GetMBeta());
    if (testSpec.IsNormallyOrdered())
      hNO.Fill(1);
    else
      hNO.Fill(0);
  }

  fout->cd();
  hEndpoint.Write();
  hMBeta.Write();
  hNO.Write();

  fout->Close();
  return 0;
}