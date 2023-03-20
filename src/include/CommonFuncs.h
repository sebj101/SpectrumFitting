/*
  CommonFuncs.h
*/

#include "src/include/FundamentalConstants.h"

namespace spec {
double DifferentialTDecayRate(double T, double E0, double mBeta) {
  double nuE{E0 - T};
  if (nuE - mBeta < 0) {
    return 0;
  } else {
    double beta{sqrt(1.0 - pow(EMASS_EV / (T + EMASS_EV), 2))};
    double p{sqrt(T * T + 2 * EMASS_EV * T)};
    double eta{ALPHA * 2 / beta};
    double fermiFunction{2 * M_PI * eta / (1 - exp(-2 * M_PI * eta))};
    double rate{G_F * G_F * pow(0.97425, 2) * fermiFunction *
                (1 + 3.0 * pow(-1.2646, 2)) * p * (T + EMASS_EV) /
                (2 * pow(M_PI, 3))};
    rate *= nuE * sqrt(nuE * nuE - mBeta * mBeta);
    rate /= 6.58e-16;
    return rate;
  }
}
}  // namespace spec