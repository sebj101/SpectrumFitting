/*
  NuFitParams.h

  Best fit oscillation parameters from the most recent NuFit 5.2 (2022)
*/

#ifndef NUFIT_PARAMS_H
#define NUFIT_PARAMS_H

#include <math.h>

namespace spec {

// Assuming normal hierarchy
inline constexpr double kNuFitTh12NH{33.41 * M_PI / 180.0};
inline constexpr double kNuFitTh23NH{49.1 * M_PI / 180.0};
inline constexpr double kNuFitTh13NH{8.54 * M_PI / 180.0};
inline constexpr double kNuFitdCPNH{197 * M_PI / 180.0};
inline constexpr double kNuFitDmsq21NH{7.41e-5};
inline constexpr double kNuFitDmsq31NH{2.511e-3};
inline constexpr double kNuFitDmsq32NH{kNuFitDmsq31NH - kNuFitDmsq21NH};
// PMNS matrix elements
inline const double Ue1SqNH{pow(cos(kNuFitTh12NH) * cos(kNuFitTh13NH), 2)};
inline const double Ue2SqNH{pow(sin(kNuFitTh12NH) * cos(kNuFitTh13NH), 2)};
inline const double Ue3SqNH{pow(sin(kNuFitTh13NH), 2)};

// Assuming inverted hierarchy
inline constexpr double kNuFitTh12IH{33.41 * M_PI / 180.0};
inline constexpr double kNuFitTh23IH{49.5 * M_PI / 180.0};
inline constexpr double kNuFitTh13IH{8.57 * M_PI / 180.0};
inline constexpr double kNuFitdCPIH{286 * M_PI / 180.0};
inline constexpr double kNuFitDmsq21IH{7.41e-5};
inline constexpr double kNuFitDmsq32IH{-2.498e-3};
// PMNS matrix elements
inline const double Ue1SqIH{pow(cos(kNuFitTh12IH) * cos(kNuFitTh13IH), 2)};
inline const double Ue2SqIH{pow(sin(kNuFitTh12IH) * cos(kNuFitTh13IH), 2)};
inline const double Ue3SqIH{pow(sin(kNuFitTh13IH), 2)};

}  // namespace spec

#endif