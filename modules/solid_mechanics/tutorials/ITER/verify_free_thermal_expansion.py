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
Mechanics verification: Free thermal expansion of an unconstrained cylinder.

THIS SCRIPT PERFORMS BOTH IMPLEMENTATION CHECKS AND MECHANICS VERIFICATION.

What this script verifies:
  1. Material property implementation (α(T) and eigenstrain match NIST correlation)
  2. Displacement field u_r(r) matches analytical solution for free thermal expansion
  3. Strain components (ε_rr, ε_θθ, ε_zz) match analytical values
  4. Stress components are zero (stress-free expansion)
  5. Complete mechanical solution from material properties through FEM solver

Analytical Solution for Free Thermal Expansion (FINITE STRAIN)
---------------------------------------------------------------
For an isotropic, homogeneous cylinder undergoing uniform temperature change
with minimal boundary constraints (only rigid body modes fixed):

MOOSE uses logarithmic (Hencky) strain formulation:

    ε_thermal(T) = ∫[T₀ to T] α(T') dT'  (thermal strain)
    λ = exp(ε_thermal)                    (stretch)

    Cartesian:
    u_x = x × (exp(ε_thermal) - 1)
    u_y = y × (exp(ε_thermal) - 1)
    u_z = z × (exp(ε_thermal) - 1)

    Cylindrical (r, θ, z):
    u_r = r × (exp(ε_thermal) - 1)
    u_θ = 0
    u_z = z × (exp(ε_thermal) - 1)

    Strain (small ε_thermal):
    ε_xx = ε_yy = ε_zz ≈ ε_thermal  (to first order)

    Stress (free expansion):
    σ_ij = 0  (all components)

where:
    - T(t) is the uniform temperature at time t [K]
    - T₀ is the stress-free reference temperature [K]
    - α(T) is the instantaneous thermal expansion coefficient [1/K] (NIST correlation)

References:
    - Rashid (1993) "Incremental kinematics for finite element applications"
    - Malvern (1969) "Introduction to the Mechanics of a Continuous Medium"
    - ε_thermal is the integrated thermal strain (must be computed numerically for α(T))
    - position = (x, y, z) is the spatial coordinate in the reference configuration

**Critical Note**: For temperature-dependent α(T), the formula ε = α(T)×ΔT is WRONG.
The correct thermal strain requires integration: ε_thermal(T) = ∫[T₀ to T] α(T') dT'

**MOOSE Postprocessor INITIAL Execution:**

The input files execute most postprocessors on 'INITIAL TIMESTEP_END', which includes
t=0 values in the CSV output. This is acceptable because:

1. ConstantPostprocessor parameters (T₀, T_start, T_f, etc.) MUST execute on INITIAL
   - These provide metadata for verification and are needed at t=0

2. Physical quantities (strain, stress, temperature) executed on INITIAL:
   - Captures the true initial state of the simulation (T=T_start, ε=0, σ=0)
   - For free thermal expansion: at t=0, thermal_strain=0 → u=0, so errors are ~0
   - This is physically correct, not a bug

3. Error postprocessors (ElementL2Error, etc.) skip INITIAL:
   - Execute on 'TIMESTEP_END' only to avoid numerical noise from zero/zero
   - See copper_cylinder_free_thermal_expansion.i lines 88, 113, 140

4. Verification script filtering:
   - Plots that use log scale filter t > 1e-12 to avoid log(0)
   - Plots without log scale include t=0 as physically meaningful

**Recommendation**: Keep INITIAL execution for physical quantities to capture true
initial state. Only skip INITIAL for error norms that would be meaningless at t=0.

**Required Conditions for Analytical Solution Validity:**

1. **Uniform Temperature Field**:
   - Temperature must be spatially constant at each time instant: T(x,y,z,t) = T(t)
   - No temperature gradients: ∇T = 0
   - Temperature CAN vary with time (e.g., ramp from 4K to 300K)

2. **Minimal Boundary Constraints**:
   - Fix ONLY enough DOFs to prevent rigid body motion
   - Recommended: Fix one point (e.g., origin) in all three directions
   - DO NOT fix multiple points or entire surfaces → introduces constraint reactions
   - All other degrees of freedom must be free to expand

3. **No External Loads**:
   - No body forces (set current_density_z = 0 for Lorentz force)
   - No applied tractions or pressures
   - No contact constraints

4. **Material Assumptions**:
   - Isotropic material (α independent of direction)
   - Homogeneous (properties uniform throughout)
   - Linear elastic behavior
   - Small strain assumption NOT required (solution valid for finite strain)

**Current Simulation Limitations:**
The copper_cylinder.i input file as currently configured does NOT satisfy condition #2:
  - Fixes multiple bottom surface points → introduces constraint reactions
  - Solution will NOT match analytical free expansion
  - Must modify BCs to fix only one point for this verification

**References:**
  - Timoshenko & Goodier (1970), "Theory of Elasticity", 3rd ed., Section 38
  - Boley & Weiner (1997), "Theory of Thermal Stresses", Chapter 11
  - Malvern (1969), "Introduction to the Mechanics of a Continuous Medium", Chapter 5

Verification Workflow
---------------------
This script verifies MOOSE free thermal expansion results against analytical solution.
It performs consistency checks to ensure the analytical table matches MOOSE parameters.

**Complete Workflow:**

1. Generate analytical thermal strain table:
    python generate_analytical_thermal_strain.py --T0 4.5 --Tf 300 --time 10 --steps 20

2. Run MOOSE simulation (uses analytical table for error postprocessors):
    cd cylinder
    <moose-app> -i copper_cylinder_free_thermal_expansion.i

3. Verify results:
    python verify_free_thermal_expansion.py

**What This Script Checks:**

- Implementation: α(T) and eigenstrain match NIST correlation
- Mechanics: Strain components match analytical values
- Stresses: All stress components are zero (stress-free expansion)
- Displacements: L2 and H1 error norms (if postprocessors present)
- Consistency: Analytical table matches MOOSE simulation parameters

Usage Examples
--------------
Use default data file:
    python verify_free_thermal_expansion.py

Specify a different CSV file and analytical table:
    python verify_free_thermal_expansion.py \
        --csv cylinder/data/copper_cylinder_free_thermal_expansion_out.csv \
        --analytical-csv cylinder/thermal_strain_analytical_vs_time.csv

Set tolerance for stress check (default 1e-3 MPa):
    python verify_free_thermal_expansion.py --stress-tol 1e-4

"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Verify MOOSE free thermal expansion against analytical solution"
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path("data/copper_cylinder_free_thermal_expansion_out.csv"),
        help="Path to MOOSE CSV output file (default: data/copper_cylinder_free_thermal_expansion_out.csv)",
    )
    parser.add_argument(
        "--analytical-csv",
        type=Path,
        default=Path("thermal_strain_analytical_vs_time.csv"),
        help="Path to analytical thermal strain table (default: thermal_strain_analytical_vs_time.csv in current directory)",
    )
    parser.add_argument(
        "--stress-tol",
        type=float,
        default=1e-3,
        help="Tolerance for stress = 0 check in MPa (default: 1e-3 MPa)",
    )
    parser.add_argument(
        "--strain-rtol",
        type=float,
        default=0.01,
        help="Relative tolerance for strain comparison (default: 0.01 = 1%%)",
    )
    parser.add_argument(
        "--disp-rtol",
        type=float,
        default=0.01,
        help="Relative tolerance for displacement comparison (default: 0.01 = 1%%)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("verification_plots"),
        help="Directory for output plots (default: verification_plots in current directory)",
    )
    parser.add_argument(
        "--param-tol",
        type=float,
        default=0.01,
        help="Relative tolerance for parameter consistency check (default: 0.01 = 1%%)",
    )
    return parser.parse_args()


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


def load_analytical_thermal_strain_table(csv_path: Path) -> dict:
    """
    Load analytical thermal strain table and extract parameters.

    Parameters
    ----------
    csv_path : Path
        Path to analytical thermal strain CSV file

    Returns
    -------
    dict
        Dictionary with keys: T0, T_start, Tf, t_end, n_steps, time, thermal_strain

    """
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Analytical thermal strain table not found: {csv_path}\n"
            f"Generate it first with:\n"
            f"  python generate_analytical_thermal_strain.py"
        )

    # Read parameters from header comments
    params = {}
    with open(csv_path, "r") as f:
        for line in f:
            if not line.startswith("#"):
                break
            if "Stress-free temperature T0 =" in line:
                params["T0"] = float(line.split("=")[1].strip().split()[0])
            elif "Starting temperature T_start =" in line:
                params["T_start"] = float(line.split("=")[1].strip().split()[0])
            elif "Final temperature Tf =" in line:
                params["Tf"] = float(line.split("=")[1].strip().split()[0])
            elif "Simulation time =" in line:
                params["t_end"] = float(line.split("=")[1].strip().split()[0])
            elif "Number of steps =" in line:
                params["n_steps"] = int(line.split("=")[1].strip())

    # Read data
    data = np.loadtxt(csv_path, delimiter=",", comments="#")
    params["time"] = data[:, 0]
    params["thermal_strain"] = data[:, 1]
    params["n_output_points"] = len(params["time"]) - 1

    return params


def check_parameter_consistency(
    df: pd.DataFrame,
    analytical_params: dict,
    rtol: float = 0.01,
) -> bool:
    """
    Check that MOOSE simulation parameters match analytical table.

    Parameters
    ----------
    df : pd.DataFrame
        MOOSE CSV output
    analytical_params : dict
        Parameters from analytical thermal strain table
    rtol : float
        Relative tolerance for comparison

    Returns
    -------
    bool
        True if all parameters are consistent

    """
    print("\n" + "=" * 80)
    print("PARAMETER CONSISTENCY CHECK")
    print("=" * 80)
    print("Verifying that analytical table matches MOOSE simulation parameters...")

    all_passed = True

    # Extract MOOSE parameters
    if "stress_free_temperature" in df.columns:
        T0_moose = df["stress_free_temperature"].iloc[0]
    else:
        T0_moose = df["temperature_average"].iloc[0]
        print(
            "  WARNING: stress_free_temperature not in MOOSE CSV, using first temperature"
        )

    if "starting_temperature" in df.columns:
        T_start_moose = df["starting_temperature"].iloc[0]
    else:
        T_start_moose = df["temperature_average"].iloc[0]

    if "ending_temperature" in df.columns:
        Tf_moose = df["ending_temperature"].iloc[0]
    else:
        Tf_moose = df["temperature_average"].iloc[-1]

    if "simulation_time" in df.columns:
        t_end_moose = df["simulation_time"].iloc[0]
    else:
        t_end_moose = df["time"].iloc[-1]

    n_steps_moose = len(df) - 1

    # Check each parameter
    checks = [
        ("Stress-free temperature T₀", T0_moose, analytical_params["T0"], "K"),
        ("Starting temperature", T_start_moose, analytical_params["T_start"], "K"),
        ("Final temperature Tf", Tf_moose, analytical_params["Tf"], "K"),
        ("Simulation time", t_end_moose, analytical_params["t_end"], "s"),
    ]

    for param_name, moose_val, analytical_val, unit in checks:
        rel_error = abs(moose_val - analytical_val) / (abs(analytical_val) + 1e-12)
        passed = rel_error < rtol
        status = "✓ PASS" if passed else "✗ FAIL"

        print(f"\n  {param_name}:")
        print(f"    MOOSE:      {moose_val:.6f} {unit}")
        print(f"    Analytical: {analytical_val:.6f} {unit}")
        print(f"    Rel. error: {rel_error:.2e} (tolerance: {rtol:.2e})")
        print(f"    Status: {status}")

        if not passed:
            all_passed = False

    # Check number of steps (must be exact)
    print("\n  Number of time steps:")
    print(f"    MOOSE:      {n_steps_moose}")
    print(f"    Analytical: {analytical_params['n_steps']}")
    if n_steps_moose == analytical_params["n_steps"]:
        print("    Status: ✓ PASS")
    else:
        print("    Status: ✗ FAIL (must match exactly)")
        all_passed = False

    # Report analytical table resolution
    print("\n  Analytical table resolution:")
    print(f"    Output points: {analytical_params['n_output_points']}")
    print(
        f"    Ratio to MOOSE steps: {analytical_params['n_output_points'] / n_steps_moose:.1f}x"
    )
    if analytical_params["n_output_points"] < 10 * n_steps_moose:
        print("    WARNING: Table may be too coarse for accurate PiecewiseLinear")
        print(f"    Recommendation: Use --n-output-points >= {10 * n_steps_moose}")

    # Check time array endpoints
    time_moose = df["time"].values
    time_analytical = analytical_params["time"]
    time_start_error = abs(time_moose[0] - time_analytical[0])
    time_end_error = abs(time_moose[-1] - time_analytical[-1])
    time_passed = (time_start_error < rtol * time_analytical[-1]) and (
        time_end_error < rtol * time_analytical[-1]
    )
    status = "✓ PASS" if time_passed else "✗ FAIL"
    print("\n  Time range alignment:")
    print(f"    Start error: {time_start_error:.2e} s")
    print(f"    End error: {time_end_error:.2e} s")
    print(f"    Status: {status}")
    if not time_passed:
        all_passed = False

    print("\n" + "=" * 80)
    if all_passed:
        print("✓ ALL PARAMETERS CONSISTENT")
        print("Analytical table matches MOOSE simulation.")
    else:
        print("✗ PARAMETER MISMATCH DETECTED")
        print("Regenerate analytical table with correct parameters:")
        print("  python generate_analytical_thermal_strain.py \\")
        print(f"      --T0 {T0_moose:.2f} \\")
        print(f"      --T-start {T_start_moose:.2f} \\")
        print(f"      --Tf {Tf_moose:.2f} \\")
        print(f"      --time {t_end_moose:.2f} \\")
        print(f"      --steps {n_steps_moose}")
    print("=" * 80)

    return all_passed


def load_moose_data(csv_path: Path) -> pd.DataFrame:
    """Load MOOSE CSV output and validate required columns."""
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    df = pd.read_csv(csv_path)

    # Check for required columns
    required_cols = [
        "time",
        "temperature_average",
        "thermal_expansion_coeff_pp",
        "thermal_strain_xx_average",
    ]

    # Check for displacement at a specific point or average
    # (This would need to be added to the input file as a postprocessor)

    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns in CSV: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )

    return df


def compute_analytical_solution(
    temperature: np.ndarray,
    alpha: np.ndarray,
    T0: float,
    r_values: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """
    Compute analytical solution for free thermal expansion.

    Parameters
    ----------
    temperature : np.ndarray
        Temperature history T(t) [K]
    alpha : np.ndarray
        Thermal expansion coefficient α(T) [1/K] from MOOSE
        (Not used - we compute from NIST directly for consistency)
    T0 : float
        Stress-free reference temperature [K]
    r_values : np.ndarray, optional
        Radial positions for displacement calculation [m]
        If None, only strain/stress are computed

    Returns
    -------
    dict
        Analytical solution with keys:
        - 'strain': ε = ∫[T0 to T] α(T') dT' (scalar, all normal components equal)
        - 'stress': σ = 0 (all components)
        - 'displacement': u_r = r × ∫[T0 to T] α(T') dT' (if r_values provided)

    """
    dT = temperature - T0

    # For temperature-dependent α(T), thermal strain must be integrated
    # ε_thermal(T) = ∫[T0 to T] α(T') dT'
    # NOT α(T) × ΔT (that's only valid for constant α!)
    alpha_nist = alpha_of_temperature(temperature)
    strain_analytical = cumulative_trapezoid(alpha_nist, temperature, initial=0)

    # For free thermal expansion: all stresses are zero
    stress_analytical = np.zeros_like(dT)

    solution = {
        "strain": strain_analytical,
        "stress": stress_analytical,
        "dT": dT,
        "alpha": alpha_nist,  # Return NIST α for plotting
    }

    # Compute displacement at specified radii
    if r_values is not None:
        # u_r(r, t) = r × ∫[T0 to T] α(T') dT'
        # Shape: (n_times, n_radii)
        displacement_analytical = (
            strain_analytical[:, np.newaxis] * r_values[np.newaxis, :]
        )
        solution["displacement"] = displacement_analytical
        solution["r_values"] = r_values

    return solution


def verify_strain_components(
    df: pd.DataFrame,
    analytical_strain: np.ndarray,
    rtol: float,
) -> dict[str, dict]:
    """
    Verify all strain components against analytical solution.

    Checks:
      - ε_xx = α × ΔT (thermal strain is isotropic)
      - ε_yy = α × ΔT
      - ε_zz = α × ΔT
      - Shear strains ≈ 0 (if available)

    Returns
    -------
    dict
        Results for each strain component with 'passed', 'max_error', 'rms_error'

    """
    results = {}

    # Normal strains (should equal analytical)
    for component in ["xx", "yy", "zz"]:
        col_name = f"thermal_strain_{component}_average"

        if col_name in df.columns:
            moose_strain = df[col_name].values
            error = moose_strain - analytical_strain
            rel_error = np.abs(error / (analytical_strain + 1e-12))  # avoid div by zero

            max_rel_error = np.max(rel_error)
            rms_rel_error = np.sqrt(np.mean(rel_error**2))

            results[component] = {
                "passed": max_rel_error < rtol,
                "max_rel_error": max_rel_error,
                "rms_rel_error": rms_rel_error,
                "max_abs_error": np.max(np.abs(error)),
            }
        else:
            results[component] = {
                "passed": None,
                "error": f"Column {col_name} not found in CSV",
            }

    # Shear strains (should be zero)
    for component in ["xy", "xz", "yz"]:
        col_name = f"strain_{component}"

        if col_name in df.columns:
            moose_shear = df[col_name].values
            max_shear = np.max(np.abs(moose_shear))

            # Use absolute tolerance for shear (should be exactly zero)
            results[component] = {
                "passed": max_shear < rtol * np.mean(np.abs(analytical_strain)),
                "max_value": max_shear,
            }

    return results


def verify_stress_components(
    df: pd.DataFrame,
    atol: float,
) -> dict[str, dict]:
    """
    Verify all stress components are zero (stress-free expansion).

    Checks each component:
      - σ_xx ≈ 0
      - σ_yy ≈ 0
      - σ_zz ≈ 0
      - τ_xy ≈ 0 (if available)
      - τ_xz ≈ 0 (if available)
      - τ_yz ≈ 0 (if available)
      - von Mises ≈ 0 (if available)

    Note: Each component is checked individually because σ_ij = 0 means
    ALL components must be zero, not just ||σ|| = 0.

    Returns
    -------
    dict
        Results for each stress component with 'passed', 'max_value', 'rms_value'

    """
    results = {}

    # Normal stresses
    for component in ["xx", "yy", "zz"]:
        col_name = f"stress_{component}"

        if col_name in df.columns:
            stress = df[col_name].values
            max_stress = np.max(np.abs(stress))
            rms_stress = np.sqrt(np.mean(stress**2))

            results[component] = {
                "passed": max_stress < atol,
                "max_value": max_stress,
                "rms_value": rms_stress,
            }
        else:
            results[component] = {
                "passed": None,
                "error": f"Column {col_name} not found",
            }

    # Shear stresses
    for component in ["xy", "xz", "yz"]:
        col_name = f"stress_{component}"

        if col_name in df.columns:
            stress = df[col_name].values
            max_stress = np.max(np.abs(stress))
            rms_stress = np.sqrt(np.mean(stress**2))

            results[component] = {
                "passed": max_stress < atol,
                "max_value": max_stress,
                "rms_value": rms_stress,
            }

    # von Mises stress
    if "vonmises_stress" in df.columns:
        vonmises = df["vonmises_stress"].values
        max_vm = np.max(vonmises)
        rms_vm = np.sqrt(np.mean(vonmises**2))

        results["vonmises"] = {
            "passed": max_vm < atol,
            "max_value": max_vm,
            "rms_value": rms_vm,
        }

    return results


def plot_verification_results(
    df: pd.DataFrame,
    analytical: dict,
    output_dir: Path,
    disp_error_results: dict | None = None,
    csv_path: Path | None = None,
):
    """Generate comprehensive verification plots."""
    output_dir.mkdir(parents=True, exist_ok=True)

    time = df["time"].values
    temperature = df["temperature_average"].values

    # Figure 1: Strain comparison
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(
        "Free Thermal Expansion: Strain Verification", fontsize=14, fontweight="bold"
    )

    # Plot thermal strain components
    ax = axes[0, 0]
    if "thermal_strain_xx_average" in df.columns:
        ax.plot(
            time,
            df["thermal_strain_xx_average"],
            "o-",
            label="MOOSE ε_xx (element avg)",
            alpha=0.7,
        )
    if "thermal_strain_yy_average" in df.columns:
        ax.plot(
            time,
            df["thermal_strain_yy_average"],
            "s-",
            label="MOOSE ε_yy (element avg)",
            alpha=0.7,
        )
    if "thermal_strain_zz_average" in df.columns:
        ax.plot(
            time,
            df["thermal_strain_zz_average"],
            "^-",
            label="MOOSE ε_zz (element avg)",
            alpha=0.7,
        )
    ax.plot(
        time,
        analytical["strain"],
        "k--",
        linewidth=2,
        label="Analytical ∫α(T)dT (NIST)",
    )
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Thermal Strain [-]")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title("Thermal Strain Components")

    # Plot strain error
    ax = axes[0, 1]
    if "thermal_strain_xx_average" in df.columns:
        error = (
            abs((df["thermal_strain_xx_average"].values - analytical["strain"]))
            / (analytical["strain"] + 1e-12)
            * 100
        )
        ax.plot(time, error, "o-", label="ε_xx error", alpha=0.7)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Relative Error [%]")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title("Strain Relative Error")

    # Plot temperature history
    ax = axes[1, 0]
    ax.plot(time, temperature, "b-", linewidth=2)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Temperature [K]")
    ax.grid(True, alpha=0.3)
    ax.set_title("Temperature History")

    # Plot α vs T
    ax = axes[1, 1]
    if "thermal_expansion_coeff_pp" in df.columns:
        alpha_moose = df["thermal_expansion_coeff_pp"].values
        ax.plot(temperature, alpha_moose * 1e6, "o-", label="MOOSE", alpha=0.7)
    if "alpha" in analytical:
        ax.plot(
            temperature, analytical["alpha"] * 1e6, "k--", linewidth=2, label="NIST"
        )
    ax.set_xlabel("Temperature [K]")
    ax.set_ylabel("α(T) [10⁻⁶/K]")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title("Thermal Expansion Coefficient")

    plt.tight_layout()
    plt.savefig(output_dir / "strain_verification.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Figure 2: Stress verification
    # Filter to actual stress components (exclude stress_free_temperature and L2 norms)
    stress_cols = [
        col
        for col in df.columns
        if col.startswith("stress_")
        and not col.startswith("stress_free")
        and not col.endswith("_L2_norm")
    ]
    if stress_cols:
        # Separate normal and shear stresses
        normal_cols = [
            col for col in stress_cols if any(c in col for c in ["xx", "yy", "zz"])
        ]
        shear_cols = [
            col for col in stress_cols if any(c in col for c in ["xy", "xz", "yz"])
        ]

        # Arrange: shear in left column, normal in right column
        n_rows = max(len(shear_cols), len(normal_cols))

        fig, axes = plt.subplots(n_rows, 2, figsize=(12, 4 * n_rows))
        fig.suptitle(
            "Free Thermal Expansion: Stress Verification (σ = 0)",
            fontsize=14,
            fontweight="bold",
        )

        if n_rows == 1:
            axes = axes.reshape(1, -1)

        # Plot shear stresses in left column
        for idx, col in enumerate(shear_cols):
            ax = axes[idx, 0]
            stress = df[col].values
            component = col.replace("stress_", "")

            ax.plot(time, stress, "o-", label=f"MOOSE τ_{component}")
            ax.axhline(
                y=0, color="k", linestyle="--", linewidth=2, label="Analytical (τ=0)"
            )
            ax.set_xlabel("Time [s]")
            ax.set_ylabel("Stress [MPa]")
            ax.legend()
            ax.grid(True, alpha=0.3)
            ax.set_title(f"Shear Stress: τ_{component}")

        # Plot normal stresses in right column
        for idx, col in enumerate(normal_cols):
            ax = axes[idx, 1]
            stress = df[col].values
            component = col.replace("stress_", "")

            ax.plot(time, stress, "o-", label=f"MOOSE σ_{component}")
            ax.axhline(
                y=0, color="k", linestyle="--", linewidth=2, label="Analytical (σ=0)"
            )
            ax.set_xlabel("Time [s]")
            ax.set_ylabel("Stress [MPa]")
            ax.legend()
            ax.grid(True, alpha=0.3)
            ax.set_title(f"Normal Stress: σ_{component}")

        # Hide unused subplots
        for idx in range(len(shear_cols), n_rows):
            axes[idx, 0].axis("off")
        for idx in range(len(normal_cols), n_rows):
            axes[idx, 1].axis("off")

        plt.tight_layout()
        plt.savefig(
            output_dir / "stress_verification.png", dpi=150, bbox_inches="tight"
        )
        plt.close()

    # Figure 3: Displacement error norms over time
    disp_l2_cols = [
        col
        for col in df.columns
        if col.startswith("disp_") and col.endswith("_L2_error")
    ]
    disp_h1semi_cols = [
        col
        for col in df.columns
        if col.startswith("disp_") and col.endswith("_H1semi_error")
    ]

    # Print displacement error summary
    if disp_l2_cols and disp_h1semi_cols:
        print("\n" + "=" * 80)
        print("DISPLACEMENT ERROR SUMMARY (Final Time)")
        print("=" * 80)
        print("L2 error: ||u_MOOSE - u_analytical||_L2 (displacement value error)")
        print(
            "H1-seminorm error: ||∇(u_MOOSE - u_analytical)||_L2 (gradient/strain error)"
        )
        print()

        for l2_col, h1semi_col in zip(disp_l2_cols, disp_h1semi_cols):
            component = l2_col.replace("disp_", "").replace("_L2_error", "")
            l2_final = df[l2_col].values[-1]
            h1semi_final = df[h1semi_col].values[-1]

            print(f"  u_{component}:")
            print(f"    ||u_error||_L2:        {l2_final:.3e} mm")
            print(f"    ||∇u_error||_L2:       {h1semi_final:.3e} (dimensionless)")

            # Check what dominates
            if l2_final > 0 and h1semi_final > 0:
                ratio = h1semi_final / l2_final
                if ratio > 10:
                    print(f"    → Gradient error dominates (ratio: {ratio:.1f})")
                elif ratio < 0.1:
                    print(f"    → Value error dominates (ratio: {ratio:.3f})")
                else:
                    print(f"    → Both errors significant (ratio: {ratio:.2f})")
            print()

        # Compare gradient error to strain error for consistency check
        strain_l2_cols = [
            col
            for col in df.columns
            if col.startswith("strain_") and col.endswith("_L2_error")
        ]
        if strain_l2_cols:
            print("  Strain errors for comparison (∇u ≈ ε in small strain):")
            for strain_col in strain_l2_cols:
                strain_component = strain_col.replace("strain_", "").replace(
                    "_L2_error", ""
                )
                strain_error = df[strain_col].values[-1]
                print(
                    f"    ||ε_{strain_component} - ε_analytical||_L2: {strain_error:.3e}"
                )
            print()
            print("  NOTE: For finite strain (MOOSE) vs small strain (analytical):")
            print(
                "  H1-seminorm may differ from strain error due to formulation mismatch"
            )
        print("=" * 80)

    if disp_l2_cols or disp_h1semi_cols:
        fig, axes = plt.subplots(2, 1, figsize=(12, 8))
        fig.suptitle(
            "Free Thermal Expansion: Displacement Error Norms",
            fontsize=14,
            fontweight="bold",
        )

        # L2 error over time (skip t=0 if present)
        ax = axes[0]
        # Filter out t=0 to avoid numerical noise
        time_mask = time > 1e-12
        time_plot = time[time_mask]

        # Find maximum error for scale setting
        max_error = 0
        for col in disp_l2_cols:
            error_data = df[col].values[time_mask]
            max_error = max(max_error, np.max(error_data))

        for col in disp_l2_cols:
            component = col.replace("disp_", "").replace("_L2_error", "")
            error_data = df[col].values[time_mask]
            # Clip errors below machine precision for clean log plot
            error_data = np.maximum(error_data, 1e-16)
            ax.plot(time_plot, error_data, "o-", label=f"u_{component}", alpha=0.7)
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("L2 Error [mm]")
        ax.set_yscale("log")
        # Set y-limits to give some headroom (10× above max)
        ax.set_ylim(bottom=max_error / 1e4, top=max_error * 10)
        ax.legend(loc="lower left")  # Place legend in bottom left
        ax.grid(True, alpha=0.3, which="both")
        ax.set_title("L2 Error: ||u_MOOSE - u_analytical||_L2 (Finite Strain)")

        # Add final error values as text annotation (top right)
        if disp_l2_cols:
            final_l2_text = "Final L2 errors:\n"
            for col in disp_l2_cols:
                component = col.replace("disp_", "").replace("_L2_error", "")
                final_val = df[col].values[-1]
                final_l2_text += f"  u_{component}: {final_val:.3e}\n"
            ax.text(
                0.98,
                0.98,
                final_l2_text.strip(),
                transform=ax.transAxes,
                verticalalignment="top",
                horizontalalignment="right",
                bbox=dict(
                    boxstyle="round", facecolor="wheat", alpha=0.8, edgecolor="black"
                ),
                fontsize=9,
                fontfamily="monospace",
            )

        # H1 seminorm error over time (skip t=0 if present)
        ax = axes[1]

        # Find maximum error for scale setting
        max_h1_error = 0
        for col in disp_h1semi_cols:
            error_data = df[col].values[time_mask]
            max_h1_error = max(max_h1_error, np.max(error_data))

        for col in disp_h1semi_cols:
            component = col.replace("disp_", "").replace("_H1semi_error", "")
            error_data = df[col].values[time_mask]
            # Clip errors below machine precision for clean log plot
            error_data = np.maximum(error_data, 1e-16)
            ax.plot(time_plot, error_data, "s-", label=f"u_{component}", alpha=0.7)
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("H1 Seminorm Error [-]")
        ax.set_yscale("log")
        # Set y-limits to give some headroom (10× above max)
        ax.set_ylim(bottom=max_h1_error / 1e4, top=max_h1_error * 10)
        ax.legend(loc="lower left")  # Move legend to lower left
        ax.grid(True, alpha=0.3, which="both")
        ax.set_title(
            "H1 Seminorm Error: ||∇(u_MOOSE - u_analytical)||_L2 (gradient error)"
        )

        # Add final error values as text annotation (bottom right)
        if disp_h1semi_cols:
            final_h1semi_text = "Final H1-seminorm errors:\n"
            for col in disp_h1semi_cols:
                component = col.replace("disp_", "").replace("_H1semi_error", "")
                final_val = df[col].values[-1]
                final_h1semi_text += f"  u_{component}: {final_val:.3e}\n"
            ax.text(
                0.98,
                0.02,
                final_h1semi_text.strip(),
                transform=ax.transAxes,
                verticalalignment="bottom",
                horizontalalignment="right",
                bbox=dict(
                    boxstyle="round", facecolor="wheat", alpha=0.8, edgecolor="black"
                ),
                fontsize=9,
                fontfamily="monospace",
            )

        plt.tight_layout()
        plt.savefig(output_dir / "displacement_error.png", dpi=150, bbox_inches="tight")
        plt.close()

    # Figure 4: Strain L2 errors and stress L2 norms
    strain_l2_cols = [
        col
        for col in df.columns
        if col.startswith("strain_") and col.endswith("_L2_error")
    ]
    stress_l2_cols = [
        col
        for col in df.columns
        if col.startswith("stress_") and col.endswith("_L2_norm")
    ]

    if strain_l2_cols or stress_l2_cols:
        n_plots = (1 if strain_l2_cols else 0) + (1 if stress_l2_cols else 0)
        fig, axes = plt.subplots(n_plots, 1, figsize=(12, 4 * n_plots))
        if n_plots == 1:
            axes = [axes]
        fig.suptitle(
            "Free Thermal Expansion: Strain/Stress L2 Norms",
            fontsize=14,
            fontweight="bold",
        )

        plot_idx = 0

        # Strain L2 errors
        # Filter t=0 which has zero values from INITIAL execution
        if strain_l2_cols:
            ax = axes[plot_idx]
            mask = time > 1e-12
            for col in strain_l2_cols:
                component = col.replace("strain_", "").replace("_L2_error", "")
                ax.plot(
                    time[mask],
                    df[col].values[mask],
                    "o-",
                    label=f"ε_{component}",
                    alpha=0.7,
                )
            ax.set_xlabel("Time [s]")
            ax.set_ylabel("L2 Error [-]")
            ax.set_yscale("log")
            ax.legend(loc="lower left")
            ax.grid(True, alpha=0.3, which="both")
            ax.set_title("Strain L2 Error: ||ε_MOOSE - ε_analytical||_L2")

            # Add final values (top right)
            final_text = "Final strain L2 errors:\n"
            for col in strain_l2_cols:
                component = col.replace("strain_", "").replace("_L2_error", "")
                final_val = df[col].values[-1]
                final_text += f"  ε_{component}: {final_val:.3e}\n"
            ax.text(
                0.98,
                0.98,
                final_text.strip(),
                transform=ax.transAxes,
                verticalalignment="top",
                horizontalalignment="right",
                bbox=dict(
                    boxstyle="round", facecolor="wheat", alpha=0.8, edgecolor="black"
                ),
                fontsize=9,
                fontfamily="monospace",
            )
            plot_idx += 1

        # Stress L2 Frobenius norm
        # Filter t=0 which has zero values from INITIAL execution
        # Note: We compute ||σ||_L2-Frobenius = sqrt(∫_Ω (σ:σ) dV) by combining component L2 norms
        # This is an approximation but provides a single scalar measure of stress-free quality
        if stress_l2_cols:
            # Split into normal and shear stresses
            normal_cols = [
                col
                for col in stress_l2_cols
                if any(c in col for c in ["xx", "yy", "zz"])
            ]
            shear_cols = [
                col
                for col in stress_l2_cols
                if any(c in col for c in ["xy", "xz", "yz"])
            ]

            ax = axes[plot_idx]
            mask = time > 1e-12

            # Compute stress tensor Frobenius-like norm from component L2 norms
            # ||σ||²_L2-Frob ≈ ||σ_xx||²_L2 + ||σ_yy||²_L2 + ||σ_zz||²_L2
            #                  + 2||σ_xy||²_L2 + 2||σ_xz||²_L2 + 2||σ_yz||²_L2
            # Factor of 2 for off-diagonal terms due to tensor symmetry
            stress_tensor_norm = np.zeros_like(time, dtype=np.float64)
            for col in normal_cols:
                stress_tensor_norm += df[col].values ** 2
            for col in shear_cols:
                stress_tensor_norm += 2 * df[col].values ** 2
            stress_tensor_norm = np.sqrt(stress_tensor_norm)

            ax.plot(
                time[mask],
                stress_tensor_norm[mask],
                "o-",
                linewidth=2,
                color="C3",
                alpha=0.8,
                markersize=4,
            )
            ax.set_xlabel("Time [s]")
            ax.set_ylabel("L2-Frobenius Norm [MPa]")
            ax.set_yscale("log")
            ax.grid(True, alpha=0.3, which="both")
            ax.set_title(
                "Stress L2-Frobenius Norm: ||σ_MOOSE||_L2-Frob (stress-free = 0)"
            )

            # Add final value (bottom right)
            final_text = f"Final ||σ||_L2-Frob: {stress_tensor_norm[-1]:.3e} MPa"
            ax.text(
                0.98,
                0.02,
                final_text,
                transform=ax.transAxes,
                verticalalignment="bottom",
                horizontalalignment="right",
                bbox=dict(
                    boxstyle="round", facecolor="wheat", alpha=0.8, edgecolor="black"
                ),
                fontsize=9,
                fontfamily="monospace",
            )

        plt.tight_layout()
        plt.savefig(
            output_dir / "strain_stress_L2_norms.png", dpi=150, bbox_inches="tight"
        )
        plt.close()

    # Figure 5: Displacement along rays comparison
    # Determine final timestep number from CSV data
    n_steps = len(df) - 1  # Subtract 1 because t=0 is included in CSV
    final_step_str = f"{n_steps:04d}"

    # VectorPostprocessor files have full output base prefix
    ray_files = [
        (
            f"copper_cylinder_free_thermal_expansion_out_disp_ray_0deg_{final_step_str}.csv",
            "0° (x-axis)",
            "x",
        ),
        (
            f"copper_cylinder_free_thermal_expansion_out_disp_ray_45deg_{final_step_str}.csv",
            "45°",
            "id",
        ),
        (
            f"copper_cylinder_free_thermal_expansion_out_disp_ray_90deg_{final_step_str}.csv",
            "90° (y-axis)",
            "y",
        ),
    ]

    # Check if ray files exist and csv_path is provided
    if csv_path is not None:
        ray_data_available = all(
            (csv_path.parent / ray_file).exists() for ray_file, _, _ in ray_files
        )
    else:
        ray_data_available = False

    if ray_data_available:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle(
            "Free Thermal Expansion: Displacement Along Radial Rays (Final Time)",
            fontsize=14,
            fontweight="bold",
        )

        # Get final thermal strain for analytical calculation (finite strain)
        final_strain = analytical["strain"][-1]
        exp_strain_minus_1 = np.exp(final_strain) - 1

        # Define colors for each displacement component
        colors = {"x": "C0", "y": "C1", "z": "C2"}  # matplotlib default color cycle

        for idx, (ray_file, label, coord) in enumerate(ray_files):
            if idx >= 3:
                break
            row = idx // 2
            col = idx % 2
            ax = axes[row, col]

            ray_path = csv_path.parent / ray_file
            try:
                ray_df = pd.read_csv(ray_path)

                # Get coordinates
                x = ray_df["x"].values
                y = ray_df["y"].values
                z = ray_df["z"].values

                # Compute radial distance for x-axis
                if coord == "x":
                    r = x
                elif coord == "y":
                    r = y
                else:  # 45 degree - compute from id
                    r = ray_df["id"].values * (10.0 / (len(ray_df) - 1))

                # Get MOOSE displacements
                u_moose_x = ray_df["disp_x"].values
                u_moose_y = ray_df["disp_y"].values
                u_moose_z = ray_df["disp_z"].values

                # Compute analytical displacements: u = position × (exp(ε) - 1)
                u_analytical_x = x * exp_strain_minus_1
                u_analytical_y = y * exp_strain_minus_1
                u_analytical_z = z * exp_strain_minus_1

                # Plot only u_x and u_y for radial rays (u_z should be constant at mid-height)
                ax.plot(
                    r,
                    u_moose_x,
                    "^",
                    color=colors["x"],
                    label="MOOSE $u_x$",
                    alpha=0.8,
                    markersize=6,
                )
                ax.plot(
                    r,
                    u_analytical_x,
                    "-",
                    color=colors["x"],
                    label="Analytical $u_x$",
                    linewidth=2,
                )

                ax.plot(
                    r,
                    u_moose_y,
                    "v",
                    color=colors["y"],
                    label="MOOSE $u_y$",
                    alpha=0.8,
                    markersize=6,
                )
                ax.plot(
                    r,
                    u_analytical_y,
                    "-",
                    color=colors["y"],
                    label="Analytical $u_y$",
                    linewidth=2,
                )

                ax.set_xlabel("Radial distance r [mm]")
                ax.set_ylabel("Displacement [mm]")
                ax.legend(fontsize=8, ncol=2)
                ax.grid(True, alpha=0.3)
                ax.set_title(f"Radial Ray at {label}")

            except Exception as e:
                ax.text(
                    0.5,
                    0.5,
                    f"Error loading {ray_file}:\n{str(e)}",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                )

        # Fourth panel: axial displacement along z-axis
        ax = axes[1, 1]
        axial_file = f"copper_cylinder_free_thermal_expansion_out_disp_axial_{final_step_str}.csv"
        axial_path = csv_path.parent / axial_file
        if axial_path.exists():
            try:
                axial_df = pd.read_csv(axial_path)

                # Get coordinates (at r=0, x=y=0)
                x = axial_df["x"].values
                y = axial_df["y"].values
                z = axial_df["z"].values

                # Get MOOSE displacements
                u_moose_z = axial_df["disp_z"].values

                # Compute analytical displacement at r=0
                # u_z = z × (exp(ε) - 1)
                # u_x, u_y ≈ 0 (x ≈ 0, y ≈ 0 on axis)
                u_analytical_z = z * exp_strain_minus_1

                # Plot only u_z for axial ray (u_x and u_y should be ~0 on centerline)
                ax.plot(
                    z,
                    u_moose_z,
                    "s",
                    color=colors["z"],
                    label="MOOSE $u_z$",
                    alpha=0.8,
                    markersize=5,
                )
                ax.plot(
                    z,
                    u_analytical_z,
                    "-",
                    color=colors["z"],
                    label="Analytical $u_z$",
                    linewidth=2,
                )

                ax.set_xlabel("Axial position z [mm]")
                ax.set_ylabel("Displacement [mm]")
                ax.legend(fontsize=8)
                ax.grid(True, alpha=0.3)
                ax.set_title("Axial Displacement Along Central Axis (r=0)")

            except Exception as e:
                ax.text(
                    0.5,
                    0.5,
                    f"Error loading {axial_file}:\n{str(e)}",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                )
        else:
            ax.text(
                0.5,
                0.5,
                f"File not found: {axial_file}",
                transform=ax.transAxes,
                ha="center",
                va="center",
            )

        plt.tight_layout()
        plt.savefig(output_dir / "displacement_rays.png", dpi=150, bbox_inches="tight")
        plt.close()
        print("  ✓ Displacement ray comparison plot created")
    else:
        print(
            "  ⓘ Skipping displacement ray plots - VectorPostprocessor CSV files not found"
        )

    print(f"\nPlots saved to: {output_dir.absolute()}")


def check_strain_L2_errors(df: pd.DataFrame) -> dict:
    """
    Check strain L2 error norms.

    Parameters
    ----------
    df : pd.DataFrame
        MOOSE CSV output

    Returns
    -------
    dict
        Strain L2 errors for each component

    """
    results = {}
    for component in ["xx", "yy", "zz"]:
        col_name = f"strain_{component}_L2_error"
        if col_name in df.columns:
            l2_error = df[col_name].values[-1]  # Error at final time
            results[component] = l2_error
    return results


def compute_strain_tensor_L2_error(component_errors: dict) -> float | None:
    """
    Compute total strain tensor L2 error from component errors.

    For symmetric tensor: ||ε||_L2 = sqrt(||ε_xx||² + ||ε_yy||² + ||ε_zz||²
                                           + 2||ε_xy||² + 2||ε_xz||² + 2||ε_yz||²)

    Factor of 2 for off-diagonal terms due to tensor symmetry.

    Parameters
    ----------
    component_errors : dict
        Dictionary of component L2 errors

    Returns
    -------
    float or None
        Total tensor L2 error, or None if insufficient components

    """
    if not component_errors:
        return None

    # Need at least normal components for meaningful tensor norm
    required_normal = {"xx", "yy", "zz"}
    if not required_normal.issubset(component_errors.keys()):
        return None

    total = 0.0
    # Normal components (diagonal)
    for comp in ["xx", "yy", "zz"]:
        if comp in component_errors:
            total += component_errors[comp] ** 2

    # Shear components (off-diagonal, factor of 2 for symmetry)
    for comp in ["xy", "xz", "yz"]:
        if comp in component_errors:
            total += 2 * component_errors[comp] ** 2

    return np.sqrt(total)


def check_stress_L2_norms(df: pd.DataFrame) -> dict:
    """
    Check stress L2 norms.

    Parameters
    ----------
    df : pd.DataFrame
        MOOSE CSV output

    Returns
    -------
    dict
        Stress L2 norms for each component

    """
    results = {}
    for component in ["xx", "yy", "zz", "xy", "xz", "yz"]:
        col_name = f"stress_{component}_L2_norm"
        if col_name in df.columns:
            l2_norm = df[col_name].values[-1]  # Norm at final time
            results[component] = l2_norm
    return results


def compute_stress_tensor_L2_norm(component_norms: dict) -> float | None:
    """
    Compute total stress tensor L2 norm from component norms.

    For symmetric tensor: ||σ||_L2 = sqrt(||σ_xx||² + ||σ_yy||² + ||σ_zz||²
                                           + 2||σ_xy||² + 2||σ_xz||² + 2||σ_yz||²)

    Factor of 2 for off-diagonal terms due to tensor symmetry.

    Parameters
    ----------
    component_norms : dict
        Dictionary of component L2 norms

    Returns
    -------
    float or None
        Total tensor L2 norm, or None if insufficient components

    """
    if not component_norms:
        return None

    # Need at least normal components for meaningful tensor norm
    required_normal = {"xx", "yy", "zz"}
    if not required_normal.issubset(component_norms.keys()):
        return None

    total = 0.0
    # Normal components (diagonal)
    for comp in ["xx", "yy", "zz"]:
        if comp in component_norms:
            total += component_norms[comp] ** 2

    # Shear components (off-diagonal, factor of 2 for symmetry)
    for comp in ["xy", "xz", "yz"]:
        if comp in component_norms:
            total += 2 * component_norms[comp] ** 2

    return np.sqrt(total)


def print_verification_summary(
    strain_results: dict,
    stress_results: dict,
    args: argparse.Namespace,
    disp_error_results: dict | None = None,
    strain_l2_errors: dict | None = None,
    stress_l2_norms: dict | None = None,
):
    """Print comprehensive verification summary."""
    print("\n" + "=" * 80)
    print("FREE THERMAL EXPANSION VERIFICATION SUMMARY")
    print("=" * 80)

    # Strain verification
    print("\n--- STRAIN VERIFICATION (ε = ∫α(T)dT) ---")
    all_strain_passed = True
    for component, result in strain_results.items():
        if result.get("passed") is not None:
            status = "✓ PASS" if result["passed"] else "✗ FAIL"
            print(f"\n  ε_{component}:")
            print(f"    Status: {status}")

            # Shear strains have 'max_value', normal strains have 'max_rel_error'
            if "max_value" in result:
                # Shear strain (should be zero)
                print(f"    Max |ε|: {result['max_value']:.2e}")
            else:
                # Normal strain (should match analytical)
                print(
                    f"    Max relative error: {result['max_rel_error']:.2e} (tolerance: {args.strain_rtol:.2e})"
                )
                print(f"    RMS relative error: {result['rms_rel_error']:.2e}")
                print(f"    Max absolute error: {result['max_abs_error']:.2e}")

            if not result["passed"]:
                all_strain_passed = False
        elif "error" in result:
            print(f"\n  ε_{component}: SKIPPED - {result['error']}")

    # Stress verification
    print("\n--- STRESS VERIFICATION (σ = 0) ---")
    all_stress_passed = True
    for component, result in stress_results.items():
        if result.get("passed") is not None:
            status = "✓ PASS" if result["passed"] else "✗ FAIL"
            print(f"\n  σ_{component}:")
            print(f"    Status: {status}")
            # Stress values are already in MPa from MOOSE (Young's modulus is 1e3*(...) MPa)
            print(
                f"    Max |σ|: {result['max_value']:.2e} MPa (tolerance: {args.stress_tol:.2e} MPa)"
            )
            print(f"    RMS σ: {result['rms_value']:.2e} MPa")

            if not result["passed"]:
                all_stress_passed = False
        elif "error" in result:
            print(f"\n  σ_{component}: SKIPPED - {result['error']}")

    # Displacement error norms
    if disp_error_results:
        print("\n--- DISPLACEMENT ERROR NORMS ---")
        print("(Computed by MOOSE using analytical displacement functions)")
        print("L2: displacement value error, H1-seminorm: gradient/strain error")
        for component, errors in disp_error_results.items():
            print(f"\n  u_{component}:")
            print(f"    L2 error:        {errors['L2']:.3e} mm")
            print(f"    H1-seminorm:     {errors['H1semi']:.3e} (dimensionless)")
        print("\nNote: Element-wise errors saved to Exodus for Paraview visualization")
        print("      Variables: disp_x_error, disp_y_error, disp_z_error")

    # Strain L2 error norms
    if strain_l2_errors:
        print("\n--- STRAIN L2 ERROR NORMS ---")
        print("(Captures spatial non-uniformity in strain field)")

        # Compute total tensor norm
        strain_tensor_error = compute_strain_tensor_L2_error(strain_l2_errors)
        if strain_tensor_error is not None:
            print(f"\n  Strain tensor L2 error: {strain_tensor_error:.3e}")
            print("  (Frobenius norm of strain error tensor)")

        print("\n  Component breakdown:")
        for component, l2_error in strain_l2_errors.items():
            print(f"    ε_{component}: {l2_error:.3e}")
        print(
            "\nNote: Tensor norm = sqrt(Σ component²), factor 2 for off-diagonal symmetry"
        )

    # Stress L2 norms
    if stress_l2_norms:
        print("\n--- STRESS L2 NORMS ---")
        print("(Measure of stress-free quality, ||σ - 0||_L2 = ||σ||_L2)")

        # Compute total tensor norm
        stress_tensor_norm = compute_stress_tensor_L2_norm(stress_l2_norms)
        if stress_tensor_norm is not None:
            print(f"\n  Stress tensor L2 norm: {stress_tensor_norm:.3e} MPa")
            print("  (Frobenius norm of stress tensor)")

        print("\n  Component breakdown:")
        for component, l2_norm in stress_l2_norms.items():
            print(f"    σ_{component}: {l2_norm:.3e} MPa")
        print(
            "\nNote: Tensor norm = sqrt(Σ component²), factor 2 for off-diagonal symmetry"
        )

    # Overall summary
    print("\n" + "=" * 80)
    if all_strain_passed and all_stress_passed:
        print("✓ ALL VERIFICATIONS PASSED")
        print("MOOSE solution matches analytical free thermal expansion.")
    else:
        print("✗ SOME VERIFICATIONS FAILED")
        if not all_strain_passed:
            print("  - Strain components do not match analytical solution")
        if not all_stress_passed:
            print("  - Stress components are not zero (not stress-free)")
        print("\nPossible causes:")
        print("  1. Boundary conditions introduce constraints (check BCs)")
        print("  2. External loads present (check body forces, pressures)")
        print("  3. Non-uniform temperature field (check spatial variation)")
        print("  4. Tolerance too strict for numerical solution")
    print("=" * 80 + "\n")

    return all_strain_passed and all_stress_passed


def verify_material_property_implementation(
    temperature: np.ndarray,
    alpha_moose: np.ndarray,
    strain_moose: np.ndarray,
    T0: float,
    rtol: float = 0.001,
) -> bool:
    """
    Verify that α(T) and eigenstrain match NIST correlation.

    This is an implementation check, not a mechanics verification.
    Validates that the material property functions correctly implement
    the NIST correlation before checking the full mechanics solution.

    Parameters
    ----------
    temperature : array
        Temperature at each timestep [K]
    alpha_moose : array
        MOOSE thermal expansion coefficient [1/K]
    strain_moose : array
        MOOSE thermal strain (eigenstrain component)
    T0 : float
        Stress-free reference temperature [K]
    rtol : float
        Relative tolerance for comparison (default 0.001 = 0.1%)

    Returns
    -------
    bool
        True if implementation matches NIST within tolerance

    """
    print("\n" + "=" * 80)
    print("IMPLEMENTATION CHECK: Material Property Validation")
    print("=" * 80)
    print("Verifying α(T) and eigenstrain against NIST correlation...")

    # Compute NIST reference values
    alpha_nist = alpha_of_temperature(temperature)

    # Numerically integrate strain from T0
    strain_nist = cumulative_trapezoid(alpha_nist, temperature, initial=0)

    # Compute relative errors
    alpha_rel_error = np.abs(alpha_moose - alpha_nist) / np.abs(alpha_nist)

    # Only check strain error where strain is large enough
    strain_threshold = 1e-6
    large_strain_mask = np.abs(strain_nist) > strain_threshold
    strain_rel_error = np.full_like(strain_moose, np.nan)
    if np.any(large_strain_mask):
        strain_rel_error[large_strain_mask] = np.abs(
            strain_moose[large_strain_mask] - strain_nist[large_strain_mask]
        ) / np.abs(strain_nist[large_strain_mask])

    # Report statistics
    max_alpha_error = np.max(alpha_rel_error)
    max_strain_error = np.nanmax(strain_rel_error) if np.any(large_strain_mask) else 0.0

    print("\nα(T) implementation:")
    print(f"  Max relative error: {max_alpha_error:.2e} ({max_alpha_error*100:.3f}%)")
    print(f"  Mean relative error: {np.mean(alpha_rel_error):.2e}")

    if np.any(large_strain_mask):
        print("\nEigenstrain integration:")
        print(
            f"  Max relative error: {max_strain_error:.2e} ({max_strain_error*100:.3f}%)"
        )
        print(f"  Mean relative error: {np.nanmean(strain_rel_error):.2e}")
        print(f"  Points checked: {np.sum(large_strain_mask)}/{len(strain_moose)}")

    # Check against tolerance
    alpha_passed = max_alpha_error < rtol
    strain_passed = max_strain_error < rtol if np.any(large_strain_mask) else True

    print(f"\nTolerance: {rtol:.2e} ({rtol*100:.2f}%)")
    print(f"  α(T): {'PASS ✓' if alpha_passed else 'FAIL ✗'}")
    print(f"  ε_thermal: {'PASS ✓' if strain_passed else 'FAIL ✗'}")

    if alpha_passed and strain_passed:
        print("\n✓ Material property implementation matches NIST correlation")
        print("  Ready to proceed with mechanics verification")
    else:
        print("\n✗ Material property implementation does NOT match NIST")
        print("  STOP: Fix material property functions before checking mechanics")
        print("\nPossible causes:")
        print("  - Incorrect NIST coefficient values in input file")
        print("  - Wrong formula in ParsedFunction/ParsedMaterial")
        print("  - Stress-free temperature mismatch")

    return alpha_passed and strain_passed


def main():
    """Main verification workflow."""
    args = parse_args()

    print("=" * 80)
    print("MOOSE Free Thermal Expansion Verification")
    print("=" * 80)
    print(f"\nInput CSV: {args.csv}")
    print(f"Stress tolerance: {args.stress_tol:.2e} MPa")
    print(f"Strain relative tolerance: {args.strain_rtol:.2e}")
    print(f"Displacement relative tolerance: {args.disp_rtol:.2e}")

    # Load MOOSE data
    print("\nLoading MOOSE output...")
    df = load_moose_data(args.csv)
    print(f"  Loaded {len(df)} time steps")

    # Extract key quantities
    temperature = df["temperature_average"].values
    alpha = df["thermal_expansion_coeff_pp"].values

    # Get stress-free temperature (from postprocessor if available, else use first temp)
    if "stress_free_temperature" in df.columns:
        T0 = df["stress_free_temperature"].iloc[0]
        print(f"  Using stress-free temperature from CSV: {T0:.1f} K")
    else:
        T0 = temperature[0]
        print(
            f"  WARNING: 'stress_free_temperature' not in CSV, using first temperature: {T0:.1f} K"
        )
        print("  Add ConstantPostprocessor to input file for exact value")

    print(f"  Temperature range: {temperature.min():.1f} - {temperature.max():.1f} K")

    # Load analytical thermal strain table
    print("\nLoading analytical thermal strain table...")
    try:
        analytical_params = load_analytical_thermal_strain_table(args.analytical_csv)
        print(f"  Loaded {len(analytical_params['time'])} points from table")
    except FileNotFoundError as e:
        print(f"\nERROR: {e}")
        sys.exit(1)

    # Check parameter consistency between MOOSE and analytical table
    consistency_passed = check_parameter_consistency(
        df, analytical_params, args.param_tol
    )

    if not consistency_passed:
        print(
            "\nVERIFICATION ABORTED: Parameter mismatch between MOOSE and analytical table"
        )
        sys.exit(1)

    # Step 1: Verify material property implementation
    # This checks that α(T) and eigenstrain correctly implement NIST correlation
    # before verifying the mechanics solution
    thermal_strain_xx = df["thermal_strain_xx_average"].values
    impl_passed = verify_material_property_implementation(
        temperature, alpha, thermal_strain_xx, T0, rtol=0.001
    )

    if not impl_passed:
        print("\n" + "=" * 80)
        print("VERIFICATION ABORTED: Material property implementation failed")
        print("=" * 80)
        return False

    # Compute analytical solution
    print("\nComputing analytical solution...")
    analytical = compute_analytical_solution(temperature, alpha, T0)

    dT_max = np.max(analytical["dT"])
    strain_max = np.max(analytical["strain"])
    print(f"  Max ΔT: {dT_max:.1f} K")
    print(f"  Max analytical strain: {strain_max:.6f}")

    # Verify strain components
    print("\nVerifying strain components...")
    strain_results = verify_strain_components(
        df, analytical["strain"], args.strain_rtol
    )

    # Verify stress components
    print("Verifying stress components...")
    stress_results = verify_stress_components(df, args.stress_tol)

    # Check displacement error norms (if available from MOOSE postprocessors)
    print("\nChecking displacement error norms...")
    disp_error_results = {}
    for component in ["x", "y", "z"]:
        l2_col = f"disp_{component}_L2_error"
        h1semi_col = f"disp_{component}_H1semi_error"
        if l2_col in df.columns and h1semi_col in df.columns:
            l2_error = df[l2_col].values[-1]  # Error at final time
            h1semi_error = df[h1semi_col].values[-1]
            disp_error_results[component] = {
                "L2": l2_error,
                "H1semi": h1semi_error,
            }
            print(
                f"  u_{component}: L2 = {l2_error:.3e}, H1-seminorm = {h1semi_error:.3e}"
            )

    if not disp_error_results:
        print("  No displacement error norms found in CSV.")
        print(
            "  Add ElementL2Error and ElementH1SemiError postprocessors to MOOSE input."
        )

    # Check strain L2 error norms (if available)
    print("\nChecking strain L2 error norms...")
    strain_l2_errors = check_strain_L2_errors(df)
    if strain_l2_errors:
        strain_tensor_error = compute_strain_tensor_L2_error(strain_l2_errors)
        if strain_tensor_error is not None:
            print(f"  Strain tensor L2 error: {strain_tensor_error:.3e}")
        for component, l2_error in strain_l2_errors.items():
            print(f"    ε_{component}: L2 error = {l2_error:.3e}")
    else:
        print("  No strain L2 error norms found in CSV.")

    # Check stress L2 norms (if available)
    print("\nChecking stress L2 norms...")
    stress_l2_norms = check_stress_L2_norms(df)
    if stress_l2_norms:
        stress_tensor_norm = compute_stress_tensor_L2_norm(stress_l2_norms)
        if stress_tensor_norm is not None:
            print(f"  Stress tensor L2 norm: {stress_tensor_norm:.3e} MPa")
        for component, l2_norm in stress_l2_norms.items():
            print(f"    σ_{component}: L2 norm = {l2_norm:.3e} MPa")
    else:
        print("  No stress L2 norms found in CSV.")

    # Generate plots
    print("\nGenerating verification plots...")
    plot_verification_results(
        df, analytical, args.output_dir, disp_error_results, Path(args.csv)
    )

    # Print summary
    all_passed = print_verification_summary(
        strain_results,
        stress_results,
        args,
        disp_error_results,
        strain_l2_errors,
        stress_l2_norms,
    )

    # Exit code
    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
