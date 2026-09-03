#!/usr/bin/env python3
"""Verify that MOOSE's thermal expansion implementation matches NIST data.

Reads the MOOSE CSV output containing temperature, the calculated thermal
expansion coefficient, and the thermal strain. Compares against the NIST
function and numerically integrated thermal strain, plotting:
  1. NIST α(T) vs MOOSE α(T) over temperature
  2. Integrated thermal strain: MOOSE vs numerical integration of NIST
  3. Relative error in thermal strain

By default, reads from data/copper_cylinder_out.csv.

Examples
--------
Use default data file:

    python verify_copper_thermal_expansion.py

Specify a different CSV file:

    python verify_copper_thermal_expansion.py \
        --csv path/to/output.csv

NIST source:
https://trc.nist.gov/cryogenics/materials/OFHC%20Copper/OFHC_Copper_rev1.htm
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid

# NIST fit: log10(alpha [10^-6/K]) = sum(c_i * log10(T)^i)
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
NIST_T_MIN_K = 4.0
NIST_T_MAX_K = 300.0


def alpha_of_temperature(temperature_k: np.ndarray) -> np.ndarray:
    """Return instantaneous linear CTE in 1/K for temperature in kelvin."""
    temperature_k = np.asarray(temperature_k, dtype=float)
    if np.any(~np.isfinite(temperature_k)):
        raise ValueError("Temperature contains a non-finite value.")
    if np.any((temperature_k < NIST_T_MIN_K) | (temperature_k > NIST_T_MAX_K)):
        low = float(np.min(temperature_k))
        high = float(np.max(temperature_k))
        raise ValueError(
            f"Temperature range [{low:g}, {high:g}] K is outside the "
            f"NIST fit range [{NIST_T_MIN_K:g}, {NIST_T_MAX_K:g}] K."
        )

    log_temperature = np.log10(temperature_k)
    # np.polynomial.polynomial.polyval expects coefficients in ascending order.
    log_alpha_micro = np.polynomial.polynomial.polyval(
        log_temperature, NIST_COEFFICIENTS
    )
    return 1.0e-6 * np.power(10.0, log_alpha_micro)


def integrate_thermal_strain(
    temperature: np.ndarray, stress_free_temperature: float
) -> np.ndarray:
    """Numerically integrate thermal strain from stress-free temperature.

    Returns epsilon_thermal = integral from T_ref to T of alpha(T') dT'
    """
    alpha = alpha_of_temperature(temperature)
    # cumulative_trapezoid with initial=0 prepends 0 and integrates using trapezoidal rule
    # Result: [0, integral(T0->T1), integral(T0->T2), ...]
    # This matches MOOSE's incremental integration starting from T0 (stress-free temp)
    thermal_strain = cumulative_trapezoid(alpha, temperature, initial=0)
    return thermal_strain


def read_moose_csv(
    path: Path,
    time_column: str,
    temperature_column: str,
    alpha_column: str,
    strain_column: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Read time, volume-averaged temperature, volume-averaged thermal expansion coefficient, and volume-averaged thermal strain from MOOSE CSV."""
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        available = reader.fieldnames or []
        requested = [time_column, temperature_column, alpha_column, strain_column]
        missing = [name for name in requested if name not in available]
        if missing:
            raise ValueError(
                f"Missing CSV column(s): {', '.join(missing)}. "
                f"Available columns: {', '.join(available)}"
            )

        rows = list(reader)

    if not rows:
        raise ValueError(f"CSV file is empty: {path}")

    time = np.array([float(row[time_column]) for row in rows])
    temperature = np.array([float(row[temperature_column]) for row in rows])
    alpha_moose = np.array([float(row[alpha_column]) for row in rows])
    strain_moose = np.array([float(row[strain_column]) for row in rows])

    order = np.argsort(time)
    return time[order], temperature[order], alpha_moose[order], strain_moose[order]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path("data/copper_cylinder_out.csv"),
        help="MOOSE CSV file containing time, temperature, thermal expansion coefficient, and thermal strain (default: data/copper_cylinder_out.csv).",
    )
    parser.add_argument("--time-column", default="time")
    parser.add_argument("--temperature-column", default="temperature_average")
    parser.add_argument(
        "--alpha-column",
        default="thermal_expansion_coeff_pp",
        help="Column containing volume-averaged CTE from material property",
    )
    parser.add_argument("--strain-column", default="thermal_strain_xx_average")
    parser.add_argument("--stress-free-temperature", type=float, default=4.5)
    parser.add_argument("--temperature-min", type=float, default=NIST_T_MIN_K)
    parser.add_argument("--temperature-max", type=float, default=NIST_T_MAX_K)
    parser.add_argument("--points", type=int, default=1000)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("copper_thermal_expansion_verification.png"),
    )
    parser.add_argument("--show", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    if args.points < 2:
        raise ValueError("--points must be at least 2.")

    # Read MOOSE output
    time, temperature_moose, alpha_moose, strain_moose = read_moose_csv(
        args.csv,
        args.time_column,
        args.temperature_column,
        args.alpha_column,
        args.strain_column,
    )

    # Check for temperatures outside NIST valid range
    valid_mask = (temperature_moose >= NIST_T_MIN_K) & (
        temperature_moose <= NIST_T_MAX_K
    )
    n_invalid = np.sum(~valid_mask)
    temp_min = float(np.min(temperature_moose))
    temp_max = float(np.max(temperature_moose))

    # Compute NIST reference only at valid MOOSE temperatures
    alpha_nist = np.full_like(temperature_moose, np.nan)
    alpha_nist[valid_mask] = alpha_of_temperature(temperature_moose[valid_mask])

    # Compute reference thermal strain via numerical integration
    # Reference strain integrated from stress-free temperature (command-line parameter)
    strain_nist = np.full_like(temperature_moose, np.nan)
    strain_nist[valid_mask] = integrate_thermal_strain(
        temperature_moose[valid_mask], args.stress_free_temperature
    )

    # Compute relative errors using magnitudes: |test - true| / |true| * 100
    alpha_error = np.abs(alpha_moose - alpha_nist) / np.abs(alpha_nist) * 100.0

    # For strain error, only compute where reference strain is large enough
    # (near stress-free temperature where strain ~ 0, relative error is meaningless)
    strain_threshold = 1e-6  # below this, thermal strain is essentially zero
    strain_error = np.full_like(strain_moose, np.nan)
    large_strain_mask = valid_mask & (np.abs(strain_nist) > strain_threshold)
    strain_error[large_strain_mask] = (
        np.abs(strain_moose[large_strain_mask] - strain_nist[large_strain_mask])
        / np.abs(strain_nist[large_strain_mask])
        * 100.0
    )

    # Generate smooth NIST curve for reference (dense temperature points)
    temperature_grid = np.geomspace(
        args.temperature_min, args.temperature_max, args.points
    )
    alpha_grid = alpha_of_temperature(temperature_grid)

    # Create figure with 2x2 panels
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 9.0), constrained_layout=True)

    # Top left: NIST fit and MOOSE volume-averaged material property vs temperature
    axes[0, 0].plot(
        temperature_grid,
        alpha_grid * 1.0e6,
        linewidth=2.0,
        label="NIST fit",
        color="C0",
    )
    axes[0, 0].plot(
        temperature_moose[valid_mask],
        alpha_moose[valid_mask] * 1.0e6,
        "o",
        markersize=4,
        label="MOOSE (volume avg)",
        alpha=0.7,
        color="C1",
    )
    axes[0, 0].set_xscale("log")
    axes[0, 0].set_xlabel("Temperature [K] (log scale)")
    axes[0, 0].set_ylabel(r"Instantaneous CTE, $\alpha$ [$10^{-6}$ K$^{-1}$]")
    axes[0, 0].set_title("Volume-averaged CTE vs NIST OFHC copper")
    axes[0, 0].grid(True, which="both", alpha=0.3)
    axes[0, 0].legend(frameon=False)

    # Bottom left: Relative error in CTE
    axes[1, 0].plot(
        temperature_moose[valid_mask],
        alpha_error[valid_mask],
        "o",
        markersize=4,
        color="C3",
    )
    axes[1, 0].set_xscale("log")
    axes[1, 0].set_xlim(axes[0, 0].get_xlim())  # Match top plot x-axis
    axes[1, 0].set_xlabel("Temperature [K] (log scale)")
    axes[1, 0].set_ylabel("Relative error [%]")
    axes[1, 0].set_title("CTE error")
    axes[1, 0].grid(True, which="both", alpha=0.3)

    # Add CTE error statistics
    valid_alpha_errors = alpha_error[valid_mask]
    if len(valid_alpha_errors) > 0:
        max_alpha_error = np.max(valid_alpha_errors)
        mean_alpha_error = np.mean(valid_alpha_errors)
        rms_alpha_error = np.sqrt(np.mean(valid_alpha_errors**2))

        alpha_stats_text = (
            f"Max: {max_alpha_error:.2e}%\n"
            f"Mean: {mean_alpha_error:.2e}%\n"
            f"RMS: {rms_alpha_error:.2e}%"
        )
        axes[1, 0].text(
            0.05,
            0.95,
            alpha_stats_text,
            transform=axes[1, 0].transAxes,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
            fontsize=9,
        )

    # Top right: Thermal strain comparison
    axes[0, 1].plot(
        temperature_moose[valid_mask],
        strain_nist[valid_mask],
        linewidth=2.0,
        label="NIST (integrated)",
        linestyle="--",
        color="C0",
    )
    axes[0, 1].plot(
        temperature_moose[valid_mask],
        strain_moose[valid_mask],
        "o",
        markersize=4,
        label="MOOSE (volume avg)",
        alpha=0.7,
        color="C1",
    )
    axes[0, 1].set_xscale("log")
    axes[0, 1].set_xlabel("Temperature [K] (log scale)")
    axes[0, 1].set_ylabel(r"Thermal strain, $\varepsilon_{xx}$ [-]")
    axes[0, 1].set_title(r"Volume-averaged thermal strain vs integrated NIST")
    axes[0, 1].grid(True, which="both", alpha=0.3)
    axes[0, 1].legend(frameon=False)

    # Bottom right: Relative error in thermal strain
    # Only plot points where strain is large enough for meaningful relative error
    plot_mask = valid_mask & (np.abs(strain_nist) > strain_threshold)
    axes[1, 1].plot(
        temperature_moose[plot_mask],
        strain_error[plot_mask],
        "o",
        markersize=4,
        color="C2",
    )
    axes[1, 1].set_xscale("log")
    axes[1, 1].set_xlim(axes[0, 1].get_xlim())  # Match top plot x-axis
    axes[1, 1].set_xlabel("Temperature [K] (log scale)")
    axes[1, 1].set_ylabel("Relative error [%]")
    axes[1, 1].set_title(r"Thermal strain $\varepsilon_{xx}$ error")
    axes[1, 1].grid(True, which="both", alpha=0.3)

    # Print statistics summary
    if len(valid_alpha_errors) > 0:
        print(f"CTE Verification Statistics ({len(valid_alpha_errors)} valid points):")
        print(f"  Maximum error: {max_alpha_error:.3e} %")
        print(f"  Mean error:    {mean_alpha_error:.3e} %")
        print(f"  RMS error:     {rms_alpha_error:.3e} %")
        print()

    # Compute and print statistics (only for points with large enough strain)
    valid_strain_errors = strain_error[large_strain_mask]
    if len(valid_strain_errors) > 0:
        max_strain_error = np.max(valid_strain_errors)
        mean_strain_error = np.mean(valid_strain_errors)
        rms_strain_error = np.sqrt(np.mean(valid_strain_errors**2))

        print(
            f"Thermal Strain Verification Statistics ({len(valid_strain_errors)} valid points):"
        )
        print(f"  Maximum error: {max_strain_error:.3e} %")
        print(f"  Mean error:    {mean_strain_error:.3e} %")
        print(f"  RMS error:     {rms_strain_error:.3e} %")

        # Add statistics as text annotation
        strain_stats_text = (
            f"Max: {max_strain_error:.2e}%\n"
            f"Mean: {mean_strain_error:.2e}%\n"
            f"RMS: {rms_strain_error:.2e}%"
        )
        axes[1, 1].text(
            0.05,
            0.95,
            strain_stats_text,
            transform=axes[1, 1].transAxes,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
            fontsize=9,
        )

    # Add warning footer if temperatures outside NIST range
    if n_invalid > 0:
        warning_text = (
            f"Warning: {n_invalid} of {len(time)} points outside NIST range [{NIST_T_MIN_K}, {NIST_T_MAX_K}] K. "
            f"Data range: [{temp_min:g}, {temp_max:g}] K."
        )
        fig.text(
            0.5,
            0.01,
            warning_text,
            ha="center",
            fontsize=9,
            color="red",
            style="italic",
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300)
    print(f"Saved {args.output}")
    if args.show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
