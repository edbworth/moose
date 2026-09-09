//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ThermalSolidProperties.h"

/**
 * Tungsten thermal properties as a function of temperature.
 *
 * Correlations from:
 * Milner, J. L., Karkos, P., & Bowers, J. J. (2024).
 * Space Nuclear Propulsion (SNP) Material Property Handbook (No. SNP-HDBK-0008).
 * National Aeronautics and Space Administration (NASA).
 * https://ntrs.nasa.gov/citations/20240004217
 */
class ThermalTungstenProperties : public ThermalSolidProperties
{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"

public:
  static InputParameters validParams();

  ThermalTungstenProperties(const InputParameters & parameters);

  virtual Real k_from_T(const Real & T) const override;
  virtual void k_from_T(const Real & T, Real & k, Real & dk_dT) const override;

  virtual Real cp_from_T(const Real & T) const override;
  virtual void cp_from_T(const Real & T, Real & cp, Real & dcp_dT) const override;

  virtual Real rho_from_T(const Real & T) const override;
  virtual void rho_from_T(const Real & T, Real & rho, Real & drho_dT) const override;

#pragma GCC diagnostic pop

protected:
  // Constants for thermal conductivity k(T)
  // Low temperature (T < 55 K): k = kA0 * (T/1000)^0.874 / (1 + kA1*(T/1000) + ...)
  static constexpr Real _kA0 = 7.348e+05; // [W/m.K]
  static constexpr Real _kA1 = 2.544e+01; // [-]
  static constexpr Real _kA2 = -8.304e+03; // [-]
  static constexpr Real _kA3 = 1.180e+06; // [-]
  static constexpr Real _kPow = 8.740e-01; // [-]

  // High temperature (T >= 55 K): k = (kB0 + kB1*t + ...) / (kC0 + kC1*t + t^2)
  static constexpr Real _kB0 = -3.679; // [W/m.K]
  static constexpr Real _kB1 = 1.181e+02; // [W/m.K]
  static constexpr Real _kB2 = 5.879e+01; // [W/m.K]
  static constexpr Real _kB3 = 2.867; // [W/m.K]
  static constexpr Real _kC0 = -2.052e-02; // [-]
  static constexpr Real _kC1 = 4.741e-01; // [-]

  // Constants for specific heat cp(T)
  // Low temperature (T <= 293 K): cp = cA0 * (T/1000)^cN / (1 + cA1*(T/1000) + ...)
  static constexpr Real _cN = 3.030; // [-]
  static constexpr Real _cA0 = 3.103e+02; // [J/g.K] = 3.103e+05 [J/kg.K]
  static constexpr Real _cA1 = -8.815; // [-]
  static constexpr Real _cA2 = 1.295e+02; // [-]
  static constexpr Real _cA3 = 1.874e+03; // [-]

  // High temperature (T > 293 K): cp = cB0 + cB1*t + cB2*t^2 + cB3*t^3 + cB_2/t^2
  static constexpr Real _cB0 = 1.301e-01; // [J/g.K] = 1.301e+02 [J/kg.K]
  static constexpr Real _cB1 = 2.225e-02; // [J/g.K] = 2.225e+01 [J/kg.K]
  static constexpr Real _cB2 = -7.224e-03; // [J/g.K] = -7.224 [J/kg.K]
  static constexpr Real _cB3 = 3.539e-03; // [J/g.K] = 3.539 [J/kg.K]
  static constexpr Real _cB_2 = -3.061e-04; // [J/g.K] = -3.061e-01 [J/kg.K]

  // Constants for density rho(T)
  // Low temperature (T <= 294 K): rho = rA0 / (1 + (rA1 + rA2*t + ...)/100)^3
  static constexpr Real _rA0 = 19250; // [kg/m^3]
  static constexpr Real _rA1 = -8.529e-02; // [-]
  static constexpr Real _rA2 = -9.915e-02; // [-]
  static constexpr Real _rA3 = 2.257; // [-]
  static constexpr Real _rA4 = -3.157; // [-]

  // High temperature (T > 294 K): rho = rA0 / (1 + (rB0 + rB1*t + ...)/100)^3
  static constexpr Real _rB0 = -1.4e-01; // [-]
  static constexpr Real _rB1 = 4.869e-01; // [-]
  static constexpr Real _rB2 = -3.056e-02; // [-]
  static constexpr Real _rB3 = 2.234e-02; // [-]
};

#pragma GCC diagnostic pop
