//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ThermalTungstenProperties.h"
#include "libmesh/utility.h"
#include <cmath>

registerMooseObject("SolidPropertiesApp", ThermalTungstenProperties);

InputParameters
ThermalTungstenProperties::validParams()
{
  InputParameters params = ThermalSolidProperties::validParams();
  params.addClassDescription("Tungsten thermal properties from NASA SNP Material Property Handbook "
                             "(Milner et al., 2024).");
  return params;
}

ThermalTungstenProperties::ThermalTungstenProperties(const InputParameters & parameters)
  : ThermalSolidProperties(parameters)
{
}

Real
ThermalTungstenProperties::k_from_T(const Real & T) const
{
  if ((T < 1.0) || (T > 3653.0))
    flagInvalidSolution("Tungsten thermal conductivity evaluated outside valid range [1, 3653] K");

  const Real t = T / 1000.0; // Scaled temperature

  if (T < 55.0)
  {
    // Low temperature: k = kA0 * t^kPow / (1 + kA1*t + kA2*t^2 + kA3*t^3)
    const Real numerator = _kA0 * std::pow(t, _kPow);
    const Real denominator = 1.0 + _kA1 * t + _kA2 * Utility::pow<2>(t) + _kA3 * Utility::pow<3>(t);
    return numerator / denominator;
  }
  else
  {
    // High temperature: k = (kB0 + kB1*t + kB2*t^2 + kB3*t^3) / (kC0 + kC1*t + t^2)
    const Real numerator = _kB0 + _kB1 * t + _kB2 * Utility::pow<2>(t) + _kB3 * Utility::pow<3>(t);
    const Real denominator = _kC0 + _kC1 * t + Utility::pow<2>(t);
    return numerator / denominator;
  }
}

void
ThermalTungstenProperties::k_from_T(const Real & T, Real & k, Real & dk_dT) const
{
  if ((T < 1.0) || (T > 3653.0))
    flagInvalidSolution("Tungsten thermal conductivity evaluated outside valid range [1, 3653] K");

  const Real t = T / 1000.0;
  const Real dt_dT = 1.0 / 1000.0;

  if (T < 55.0)
  {
    // Low temperature: k = kA0 * t^kPow / (1 + kA1*t + kA2*t^2 + kA3*t^3)
    const Real numerator = _kA0 * std::pow(t, _kPow);
    const Real dnumerator_dt = _kA0 * _kPow * std::pow(t, _kPow - 1.0);
    const Real denominator = 1.0 + _kA1 * t + _kA2 * Utility::pow<2>(t) + _kA3 * Utility::pow<3>(t);
    const Real ddenominator_dt = _kA1 + 2.0 * _kA2 * t + 3.0 * _kA3 * Utility::pow<2>(t);

    k = numerator / denominator;
    // Quotient rule: d(u/v)/dt = (u'v - uv')/v^2
    const Real dk_dt = (dnumerator_dt * denominator - numerator * ddenominator_dt) /
                       Utility::pow<2>(denominator);
    dk_dT = dk_dt * dt_dT;
  }
  else
  {
    // High temperature: k = (kB0 + kB1*t + kB2*t^2 + kB3*t^3) / (kC0 + kC1*t + t^2)
    const Real numerator = _kB0 + _kB1 * t + _kB2 * Utility::pow<2>(t) + _kB3 * Utility::pow<3>(t);
    const Real dnumerator_dt = _kB1 + 2.0 * _kB2 * t + 3.0 * _kB3 * Utility::pow<2>(t);
    const Real denominator = _kC0 + _kC1 * t + Utility::pow<2>(t);
    const Real ddenominator_dt = _kC1 + 2.0 * t;

    k = numerator / denominator;
    const Real dk_dt = (dnumerator_dt * denominator - numerator * ddenominator_dt) /
                       Utility::pow<2>(denominator);
    dk_dT = dk_dt * dt_dT;
  }
}

Real
ThermalTungstenProperties::cp_from_T(const Real & T) const
{
  if ((T < 11.0) || (T > 3700.0))
    flagInvalidSolution("Tungsten specific heat evaluated outside valid range [11, 3700] K");

  const Real t = T / 1000.0;

  if (T <= 293.0)
  {
    // Low temperature: cp = cA0 * t^cN / (1 + cA1*t + cA2*t^2 + cA3*t^3)
    // Note: cA0 is in J/g.K, convert to J/kg.K by multiplying by 1000
    const Real numerator = _cA0 * 1000.0 * std::pow(t, _cN);
    const Real denominator = 1.0 + _cA1 * t + _cA2 * Utility::pow<2>(t) + _cA3 * Utility::pow<3>(t);
    return numerator / denominator;
  }
  else
  {
    // High temperature: cp = cB0 + cB1*t + cB2*t^2 + cB3*t^3 + cB_2/t^2
    // Note: coefficients are in J/g.K, convert to J/kg.K by multiplying by 1000
    return 1000.0 * (_cB0 + _cB1 * t + _cB2 * Utility::pow<2>(t) + _cB3 * Utility::pow<3>(t) +
                     _cB_2 / Utility::pow<2>(t));
  }
}

void
ThermalTungstenProperties::cp_from_T(const Real & T, Real & cp, Real & dcp_dT) const
{
  if ((T < 11.0) || (T > 3700.0))
    flagInvalidSolution("Tungsten specific heat evaluated outside valid range [11, 3700] K");

  const Real t = T / 1000.0;
  const Real dt_dT = 1.0 / 1000.0;

  if (T <= 293.0)
  {
    // Low temperature: cp = cA0*1000 * t^cN / (1 + cA1*t + cA2*t^2 + cA3*t^3)
    const Real numerator = _cA0 * 1000.0 * std::pow(t, _cN);
    const Real dnumerator_dt = _cA0 * 1000.0 * _cN * std::pow(t, _cN - 1.0);
    const Real denominator = 1.0 + _cA1 * t + _cA2 * Utility::pow<2>(t) + _cA3 * Utility::pow<3>(t);
    const Real ddenominator_dt = _cA1 + 2.0 * _cA2 * t + 3.0 * _cA3 * Utility::pow<2>(t);

    cp = numerator / denominator;
    const Real dcp_dt = (dnumerator_dt * denominator - numerator * ddenominator_dt) /
                        Utility::pow<2>(denominator);
    dcp_dT = dcp_dt * dt_dT;
  }
  else
  {
    // High temperature: cp = 1000*(cB0 + cB1*t + cB2*t^2 + cB3*t^3 + cB_2/t^2)
    cp = 1000.0 * (_cB0 + _cB1 * t + _cB2 * Utility::pow<2>(t) + _cB3 * Utility::pow<3>(t) +
                   _cB_2 / Utility::pow<2>(t));
    const Real dcp_dt = 1000.0 * (_cB1 + 2.0 * _cB2 * t + 3.0 * _cB3 * Utility::pow<2>(t) -
                                  2.0 * _cB_2 / Utility::pow<3>(t));
    dcp_dT = dcp_dt * dt_dT;
  }
}

Real
ThermalTungstenProperties::rho_from_T(const Real & T) const
{
  if ((T < 5.0) || (T > 3600.0))
    flagInvalidSolution("Tungsten density evaluated outside valid range [5, 3600] K");

  const Real t = T / 1000.0;

  if (T <= 294.0)
  {
    // Low temperature: rho = rA0 / (1 + (rA1 + rA2*t + rA3*t^2 + rA4*t^3)/100)^3
    const Real expansion_factor =
        1.0 + (_rA1 + _rA2 * t + _rA3 * Utility::pow<2>(t) + _rA4 * Utility::pow<3>(t)) / 100.0;
    return _rA0 / Utility::pow<3>(expansion_factor);
  }
  else
  {
    // High temperature: rho = rA0 / (1 + (rB0 + rB1*t + rB2*t^2 + rB3*t^3)/100)^3
    const Real expansion_factor =
        1.0 + (_rB0 + _rB1 * t + _rB2 * Utility::pow<2>(t) + _rB3 * Utility::pow<3>(t)) / 100.0;
    return _rA0 / Utility::pow<3>(expansion_factor);
  }
}

void
ThermalTungstenProperties::rho_from_T(const Real & T, Real & rho, Real & drho_dT) const
{
  if ((T < 5.0) || (T > 3600.0))
    flagInvalidSolution("Tungsten density evaluated outside valid range [5, 3600] K");

  const Real t = T / 1000.0;
  const Real dt_dT = 1.0 / 1000.0;

  if (T <= 294.0)
  {
    // Low temperature: rho = rA0 / (1 + (rA1 + rA2*t + rA3*t^2 + rA4*t^3)/100)^3
    const Real polynomial = _rA1 + _rA2 * t + _rA3 * Utility::pow<2>(t) + _rA4 * Utility::pow<3>(t);
    const Real dpolynomial_dt = _rA2 + 2.0 * _rA3 * t + 3.0 * _rA4 * Utility::pow<2>(t);
    const Real expansion_factor = 1.0 + polynomial / 100.0;

    rho = _rA0 / Utility::pow<3>(expansion_factor);
    // d(1/f^3)/dt = -3/f^4 * df/dt
    const Real dexpansion_dt = dpolynomial_dt / 100.0;
    const Real drho_dt = -3.0 * _rA0 / Utility::pow<4>(expansion_factor) * dexpansion_dt;
    drho_dT = drho_dt * dt_dT;
  }
  else
  {
    // High temperature: rho = rA0 / (1 + (rB0 + rB1*t + rB2*t^2 + rB3*t^3)/100)^3
    const Real polynomial = _rB0 + _rB1 * t + _rB2 * Utility::pow<2>(t) + _rB3 * Utility::pow<3>(t);
    const Real dpolynomial_dt = _rB1 + 2.0 * _rB2 * t + 3.0 * _rB3 * Utility::pow<2>(t);
    const Real expansion_factor = 1.0 + polynomial / 100.0;

    rho = _rA0 / Utility::pow<3>(expansion_factor);
    const Real dexpansion_dt = dpolynomial_dt / 100.0;
    const Real drho_dt = -3.0 * _rA0 / Utility::pow<4>(expansion_factor) * dexpansion_dt;
    drho_dT = drho_dt * dt_dT;
  }
}
