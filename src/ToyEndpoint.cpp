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

// Main
int main(int argc, char *argv[]) {
  cout.precision(10);

  TString outputFile{argv[1]};
  auto fout = make_unique<TFile>(outputFile, "RECREATE");

  /////////////// Spectrum testing /////////////////////
  const double nTestAtoms{1e20};
  const double times{1};           // years
  const double specSize{100};      // eV
  const double nuMassLight{0.05};  // eV
  const double nuMassHeavy{1};     // eV
  const double nuMassMassless{0};  // eV

  Spectrum testSpecHeavy(true, nuMassHeavy, times, nTestAtoms, specSize);
  auto specHistHeavy = testSpecHeavy.GetSpectrum();
  specHistHeavy.SetLineColor(kRed);
  Spectrum testSpecLight(true, nuMassLight, times, nTestAtoms, specSize);
  auto specHistLight = testSpecLight.GetSpectrum();
  specHistLight.SetLineColor(kCyan + 1);
  Spectrum testSpecMassless(true, nuMassMassless, times, nTestAtoms, specSize);
  auto specHistMassless = testSpecMassless.GetSpectrum();
  specHistMassless.SetLineColor(kBlack);

  fout->cd();
  specHistHeavy.Write("specHistHeavy");
  specHistLight.Write("specHistLight");
  specHistMassless.Write("specHistMassless");

  fout->Close();
  return 0;
}