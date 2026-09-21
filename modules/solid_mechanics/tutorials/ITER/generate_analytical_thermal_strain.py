#!/usr/bin/env python3
# This file is part of the MOOSE framework
# https://mooseframework.inl.gov
#
# All rights reserved, see COPYRIGHT for full restrictions
# https://github.com/idaholab/moose/blob/master/COPYRIGHT
#
# Licensed under LGPL 2.1, please see LICENSE for details
# https://www.gnu.org/licenses/lgpl-2.1.html

"""
Generate analytical thermal strain table for free thermal expansion verification.

This script computes thermal_strain(T) = ∫[T0 to T] α(T') dT' using the NIST
correlation for OFHC copper with fine numerical quadrature. The output CSV file
is used by MOOSE for analytical displacement functions and error postprocessors.

The analytical solution is INDEPENDENT of MOOSE material property implementations
(e.g., ADComputeInstantaneousThermalExpansionFunctionEigenstrain) - this script
verifies that those implementations are correct.

Usage
-----
Specify parameters explicitly:
    python generate_analytical_thermal_strain.py --T0 4.5 --Tf 300 --time 10 --steps 20

Output to custom location:
    python generate_analytical_thermal_strain.py --T0 4.5 --Tf 300 --time 10 --steps 20 \\
        --output my_analysis/thermal_strain.csv

Use finer quadrature (default 1000 points):
    python generate_analytical_thermal_strain.py --T0 4.5 --Tf 300 --time 10 --steps 20 \\
        --n-quadrature 5000
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid

# NIST thermal expansion coefficient for OFHC Copper
# Source: https://trc.nist.gov/cryogenics/materials/OFHC%20Copper/OFHC_Copper_rev1.htm
NIST_COEFFICIENTS = np.array(
    [
        -17.9081289,
        67.131914,
        -118.809316,
        109.9845997,
        -53.8696089,
        13.30247491,
        -1.30843441,
    ]
)


def alpha_of_temperature(temperature_k: np.ndarray) -> np.ndarray:
    """
    Return instantaneous linear CTE in 1/K for temperature in kelvin.

    Uses NIST correlation for OFHC copper, valid 4-300 K.
    """
    log_temperature = np.log10(temperature_k)
    log_alpha_micro = np.polynomial.polynomial.polyval(
        log_temperature, NIST_COEFFICIENTS
    )
    return 1.0e-6 * np.power(10.0, log_alpha_micro)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate analytical thermal strain table for MOOSE verification"
    )
    parser.add_argument(
        "--T0",
        type=float,
        default=4.5,
        help="Stress-free reference temperature [K] (default: 4.5 K)",
    )
    parser.add_argument(
        "--T-start",
        type=float,
        default=None,
        help="Starting temperature [K] (default: same as T0)",
    )
    parser.add_argument(
        "--Tf",
        type=float,
        default=300.0,
        help="Final temperature [K] (default: 300.0 K)",
    )
    parser.add_argument(
        "--time",
        type=float,
        default=10.0,
        help="Simulation end time [s] (default: 10.0 s)",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=100,
        help="Number of MOOSE time steps for reference (default: 100, but table will be much finer)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("thermal_strain_analytical_vs_time.csv"),
        help="Output CSV file path (default: thermal_strain_analytical_vs_time.csv in current directory)",
    )
    parser.add_argument(
        "--n-quadrature",
        type=int,
        default=10000,
        help="Number of quadrature points for numerical integration (default: 1000)",
    )
    parser.add_argument(
        "--n-output-points",
        type=int,
        default=400,
        help="Number of points in output CSV table for MOOSE PiecewiseLinear (default: 200, must be >> steps for accurate interpolation)",
    )
    return parser.parse_args()


def generate_analytical_thermal_strain_table(
    T0: float,
    T_start: float,
    Tf: float,
    t_end: float,
    n_steps: int,
    output_file: Path,
    n_quadrature: int = 1000,
    n_output_points: int = 200,
) -> None:
    """
    Generate analytical thermal strain vs. time table for MOOSE.

    Computes thermal_strain(T) = ∫[T0 to T] α(T') dT' using NIST correlation
    with fine quadrature, independent of MOOSE material objects.

    Uses two-stage process:
    1. Fine quadrature integration (1000+ points) for accurate ∫α(T)dT
    2. Output to fine table (200+ points) for MOOSE PiecewiseLinear

    The output table is much finer than MOOSE time steps so that linear
    interpolation in PiecewiseLinear accurately represents the nonlinear
    thermal strain function.

    Parameters
    ----------
    T0 : float
        Stress-free reference temperature [K]
    T_start : float
        Starting temperature [K]
    Tf : float
        Final temperature [K]
    t_end : float
        Simulation end time [s]
    n_steps : int
        Number of MOOSE time steps (for reference only)
    output_file : Path
        Output CSV file path
    n_quadrature : int
        Number of quadrature points for numerical integration (default: 1000)
    n_output_points : int
        Number of points in output table (default: 200, must be >> n_steps)

    """
    print("=" * 80)
    print("GENERATING ANALYTICAL THERMAL STRAIN TABLE")
    print("=" * 80)
    print("Simulation parameters:")
    print(f"  Stress-free temperature T₀: {T0:.2f} K")
    print(f"  Starting temperature: {T_start:.2f} K")
    print(f"  Final temperature Tf: {Tf:.2f} K")
    print(f"  Simulation time: {t_end:.2f} s")
    print(f"  Number of MOOSE time steps (for reference only): {n_steps}")
    print(f"  Quadrature points: {n_quadrature}")

    if n_output_points <= n_steps:
        print(
            f"\n  WARNING: Output points ({n_output_points}) should be much larger than"
        )
        print(
            f"           MOOSE steps ({n_steps}) for accurate PiecewiseLinear interpolation."
        )
        print("  Recommendation: n_output_points >= 10 * n_steps")

    # Time discretization for output table (FINE grid for MOOSE PiecewiseLinear)
    time = np.linspace(0, t_end, n_output_points + 1)

    # Temperature history (linear ramp from T_start to Tf)
    temperature = T_start + (Tf - T_start) * time / t_end

    # Fine temperature quadrature for integration
    # Use much finer grid than MOOSE time steps to capture α(T) variation accurately
    T_quad = np.linspace(T0, Tf, n_quadrature)
    alpha_quad = alpha_of_temperature(T_quad)

    # Integrate thermal strain using fine quadrature
    thermal_strain_quad = cumulative_trapezoid(alpha_quad, T_quad, initial=0)

    # Interpolate to MOOSE time step temperatures
    thermal_strain = np.interp(temperature, T_quad, thermal_strain_quad)

    print("\nAnalytical thermal strain range:")
    print(f"  At T = {temperature[0]:.2f} K: {thermal_strain[0]:.6e}")
    print(f"  At T = {temperature[-1]:.2f} K: {thermal_strain[-1]:.6e}")
    print(f"  Max thermal strain: {np.max(thermal_strain):.6e}")

    # Write to CSV for MOOSE PiecewiseLinear function
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as f:
        f.write(
            "# Analytical thermal strain vs. time for free thermal expansion verification\n"
        )
        f.write("# Generated from NIST OFHC Copper correlation (4-300 K)\n")
        f.write(f"# Stress-free temperature T0 = {T0:.2f} K\n")
        f.write(f"# Starting temperature T_start = {T_start:.2f} K\n")
        f.write(f"# Final temperature Tf = {Tf:.2f} K\n")
        f.write(f"# Simulation time = {t_end:.2f} s\n")
        f.write(f"# Number of steps = {n_steps}\n")
        f.write(f"# Integration quadrature points: {n_quadrature}\n")
        f.write("# thermal_strain(T) = ∫[T0 to T] α(T') dT'\n")
        f.write("# time [s], thermal_strain [-]\n")
        for t, strain in zip(time, thermal_strain):
            f.write(f"{t:.10e},{strain:.15e}\n")

    print(f"\nAnalytical thermal strain table written to: {output_file.absolute()}")
    print("=" * 80)
    print("\nNext steps:")
    print("  1. Run MOOSE simulation with input file that references this table")
    print("  2. MOOSE will use it for analytical displacement functions")
    print("  3. MOOSE will compute L2/H1 error norms and write to CSV")
    print("  4. Run verify_free_thermal_expansion.py to check results")


def main():
    """Main entry point."""
    args = parse_args()

    # Use T0 as starting temperature if not specified
    T_start = args.T_start if args.T_start is not None else args.T0

    generate_analytical_thermal_strain_table(
        T0=args.T0,
        T_start=T_start,
        Tf=args.Tf,
        t_end=args.time,
        n_steps=args.steps,
        output_file=args.output,
        n_quadrature=args.n_quadrature,
        n_output_points=args.n_output_points,
    )


if __name__ == "__main__":
    main()
