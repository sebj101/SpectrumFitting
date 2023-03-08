/*
  FundamentalConstants.h
*/

#ifndef FUNDAMENTAL_CONSTANTS_H
#define FUNDAMENTAL_CONSTANTS_H

namespace spec {
// Permeability of free space
inline constexpr double EPSILON0{8.8541878128e-12};

// Permeability of free space
inline constexpr double MU0{1.25663706212e-6};

// Elementary charge in coulombs
inline constexpr double QE{1.602176634e-19};

// Speed of light in a vacuum
inline constexpr double CLIGHT{299792458};

// Electron rest mass in kilograms
inline constexpr double EMASS{9.1093837015e-31};

// Electron rest mass in eV
inline constexpr double EMASS_EV{EMASS * CLIGHT * CLIGHT / QE};

// Atomic mass unit/Dalton in kg
inline constexpr double DALTON{1.6605390666e-27};

// Boltzmann constant
inline constexpr double KB{1.380649e-23};  // kg m^2 s^-2 K^-1

// Fine structure constant
inline constexpr double ALPHA{7.2973525698e-3};

// Fermi coupling constant
inline constexpr double G_F{1.1663787e-23};
}  // namespace spec

#endif