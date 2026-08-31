#!/usr/bin/env python3
"""
Compare MOOSE simulation results to analytical Lorentz force solution.
Comprehensive verification including axisymmetry and axial uniformity checks.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Physical constants from the input file
vacuum_permeability = 1.25663706e-6  # N/A^2
current_density_z = 1.0e02  # A/mm^2 (after unit conversion from 1e6 A/m^2)


def analytical_lorentz_radial(r):
    """
    Analytical solution for radial Lorentz force per unit volume.
    F_r = -μ₀ * J_z^2 * r / 2

    Parameters:
    -----------
    r : array-like
        Radial distance from center (mm)

    Returns:
    --------
    F_r : array-like
        Radial Lorentz force per unit volume (N/mm^3), negative indicates compression
    """
    return -vacuum_permeability * current_density_z**2 * r / 2


def process_radial_line(filename, line_name):
    """
    Process one radial line sample.

    Parameters:
    -----------
    filename : str
        CSV file to read
    line_name : str
        Descriptive name for this line

    Returns:
    --------
    dict : Contains r, lorentz_radial_sim, lorentz_radial_analytical, errors
    """
    try:
        data = pd.read_csv(filename)
    except FileNotFoundError:
        print(f"Error: Could not find '{filename}'")
        return None

    # Extract coordinates and compute radius
    x = data["x"].values
    y = data["y"].values
    r = np.sqrt(x**2 + y**2)

    # Extract Cartesian components from MOOSE simulation
    lorentz_x_sim = data["lorentz_x_aux"].values
    lorentz_y_sim = data["lorentz_y_aux"].values

    # Mask for non-zero radius to avoid division by zero
    mask = r > 1e-10

    # Compute radial component from simulation (dot product: F·r̂ = (Fx*x + Fy*y)/r)
    lorentz_radial_sim = np.zeros_like(r)
    lorentz_radial_sim[mask] = (
        lorentz_x_sim[mask] * x[mask] + lorentz_y_sim[mask] * y[mask]
    ) / r[mask]
    print(f"Lorentz force value at r=R: {1e9*lorentz_radial_sim[mask][-1]:.6e} N/m^3")

    # Compute analytical radial solution
    lorentz_radial_analytical = analytical_lorentz_radial(r)

    # Compute errors (only for non-zero radius)
    absolute_error = np.zeros_like(r)
    relative_error = np.zeros_like(r)
    absolute_error[mask] = np.abs(
        lorentz_radial_sim[mask] - lorentz_radial_analytical[mask]
    )
    relative_error[mask] = absolute_error[mask] / np.abs(
        lorentz_radial_analytical[mask]
    )

    return {
        "name": line_name,
        "r": r,
        "mask": mask,
        "lorentz_sim": lorentz_radial_sim,
        "lorentz_analytical": lorentz_radial_analytical,
        "abs_error": absolute_error,
        "rel_error": relative_error,
    }


def process_axial_line(filename):
    """
    Process the axial line sample to verify uniformity along z.

    Parameters:
    -----------
    filename : str
        CSV file to read

    Returns:
    --------
    dict : Contains z, r, lorentz_radial values
    """
    try:
        data = pd.read_csv(filename)
    except FileNotFoundError:
        print(f"Error: Could not find '{filename}'")
        return None

    # Extract coordinates
    x = data["x"].values
    y = data["y"].values
    z = data["z"].values
    r = np.sqrt(x**2 + y**2)

    # Extract Cartesian components
    lorentz_x_sim = data["lorentz_x_aux"].values
    lorentz_y_sim = data["lorentz_y_aux"].values

    # Compute radial component
    mask = r > 1e-10
    lorentz_radial_sim = np.zeros_like(r)
    lorentz_radial_sim[mask] = (
        lorentz_x_sim[mask] * x[mask] + lorentz_y_sim[mask] * y[mask]
    ) / r[mask]

    # Analytical solution at this constant radius
    lorentz_radial_analytical = analytical_lorentz_radial(r)

    return {
        "z": z,
        "r": r,
        "lorentz_sim": lorentz_radial_sim,
        "lorentz_analytical": lorentz_radial_analytical,
    }


# Process radial lines at different angles
print("=" * 70)
print("Lorentz Force Verification: MOOSE vs. Analytical")
print("=" * 70)

radial_lines = []
line_configs = [
    ("data/copper_cylinder_out_line_sample_0deg_0001.csv", "0° (x-axis)"),
    ("data/copper_cylinder_out_line_sample_45deg_0001.csv", "45°"),
    ("data/copper_cylinder_out_line_sample_90deg_0001.csv", "90° (y-axis)"),
]

for filename, name in line_configs:
    result = process_radial_line(filename, name)
    if result is not None:
        radial_lines.append(result)
    else:
        print(f"Warning: Skipping {name} due to missing file")

if not radial_lines:
    print("\nError: No radial line samples found. Please run the simulation first:")
    print("  <your-moose-app> -i copper_cylinder.i")
    exit(1)

# Print statistics for each radial line
print("\n1. RADIAL LINES (Axisymmetry Check)")
print("-" * 70)
for line in radial_lines:
    mask = line["mask"]
    max_abs_error = np.max(line["abs_error"][mask])
    mean_abs_error = np.mean(line["abs_error"][mask])
    max_rel_error = np.max(line["rel_error"][mask]) * 100
    mean_rel_error = np.mean(line["rel_error"][mask]) * 100

    print(f"\n{line['name']}:")
    print(f"  Number of points:    {np.sum(mask)}")
    print(
        f"  Radial range:        {line['r'][mask].min():.4f} to {line['r'][mask].max():.4f} mm"
    )
    print(f"  Max absolute error:  {max_abs_error:.6e} N/mm³")
    print(f"  Mean absolute error: {mean_abs_error:.6e} N/mm³")
    print(f"  Max relative error:  {max_rel_error:.4f} %")
    print(f"  Mean relative error: {mean_rel_error:.4f} %")

# Process axial line
print("\n2. AXIAL LINE (Uniformity Check)")
print("-" * 70)
axial_result = process_axial_line("data/copper_cylinder_out_line_sample_axial_0001.csv")
if axial_result is not None:
    axial_abs_error = np.abs(
        axial_result["lorentz_sim"] - axial_result["lorentz_analytical"]
    )
    axial_variation = np.std(axial_result["lorentz_sim"])
    axial_mean = np.mean(axial_result["lorentz_sim"])

    print(f"\nAxial line at r = {axial_result['r'][0]:.4f} mm:")
    print(f"  Number of points:       {len(axial_result['z'])}")
    print(
        f"  Axial range:            {axial_result['z'].min():.4f} to {axial_result['z'].max():.4f} mm"
    )
    print(f"  Mean Lorentz force:     {axial_mean:.6e} N/mm³")
    print(f"  Std dev along axis:     {axial_variation:.6e} N/mm³")
    print(
        f"  Analytical value:       {axial_result['lorentz_analytical'][0]:.6e} N/mm³"
    )
    print(f"  Max absolute error:     {np.max(axial_abs_error):.6e} N/mm³")
    print(f"  Mean absolute error:    {np.mean(axial_abs_error):.6e} N/mm³")

    # Check if axial variation is small (indicates uniformity)
    relative_variation = (
        axial_variation / np.abs(axial_mean) * 100 if axial_mean != 0 else 0
    )
    print(f"  Relative variation:     {relative_variation:.4f} %")
    if relative_variation < 1.0:
        print("  ✓ Axial uniformity confirmed (variation < 1%)")
    else:
        print("  ⚠ Significant axial variation detected")
else:
    print("\nWarning: Axial line sample not found")

print("\n" + "=" * 70)

# Create comprehensive comparison plots
fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)

# Plot 1: Radial Lorentz force comparison for all angles
ax1 = fig.add_subplot(gs[0, 0])
colors = ["red", "green", "blue"]
for i, line in enumerate(radial_lines):
    mask = line["mask"]
    ax1.plot(
        line["r"][mask],
        line["lorentz_sim"][mask],
        "o",
        color=colors[i],
        markersize=3,
        alpha=0.6,
        label=f"MOOSE {line['name']}",
    )
if radial_lines:
    mask = radial_lines[0]["mask"]
    ax1.plot(
        radial_lines[0]["r"][mask],
        radial_lines[0]["lorentz_analytical"][mask],
        "k-",
        linewidth=2,
        label="Analytical",
    )
ax1.set_xlabel("Radius (mm)", fontsize=11)
ax1.set_ylabel("Radial Lorentz Force (N/mm³)", fontsize=11)
ax1.set_title("Radial Force: Axisymmetry Check", fontsize=12, fontweight="bold")
ax1.legend(fontsize=9, loc="best")
ax1.grid(True, alpha=0.3)

# Plot 2: Absolute error for all radial lines
ax2 = fig.add_subplot(gs[0, 1])
for i, line in enumerate(radial_lines):
    mask = line["mask"]
    ax2.plot(
        line["r"][mask],
        line["abs_error"][mask],
        color=colors[i],
        linewidth=1.5,
        label=line["name"],
    )
ax2.set_xlabel("Radius (mm)", fontsize=11)
ax2.set_ylabel("Absolute Error (N/mm³)", fontsize=11)
ax2.set_title("Absolute Error by Angle", fontsize=12, fontweight="bold")
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)
ax2.axhline(y=0, color="k", linestyle="--", linewidth=1)

# Plot 3: Relative error for all radial lines
ax3 = fig.add_subplot(gs[1, 0])
for i, line in enumerate(radial_lines):
    mask = line["mask"]
    ax3.plot(
        line["r"][mask],
        line["rel_error"][mask] * 100,
        color=colors[i],
        linewidth=1.5,
        label=line["name"],
    )
ax3.set_xlabel("Radius (mm)", fontsize=11)
ax3.set_ylabel("Relative Error (%)", fontsize=11)
ax3.set_title("Relative Error by Angle", fontsize=12, fontweight="bold")
ax3.legend(fontsize=9)
ax3.grid(True, alpha=0.3)
ax3.axhline(y=0, color="k", linestyle="--", linewidth=1)

# Plot 4: Axial uniformity check
ax4 = fig.add_subplot(gs[1, 1])
if axial_result is not None:
    ax4.plot(
        axial_result["z"], axial_result["lorentz_sim"], "b-", linewidth=2, label="MOOSE"
    )
    ax4.axhline(
        y=axial_result["lorentz_analytical"][0],
        color="k",
        linestyle="--",
        linewidth=2,
        label="Analytical",
    )
    ax4.set_xlabel("Axial Position z (mm)", fontsize=11)
    ax4.set_ylabel("Radial Lorentz Force (N/mm³)", fontsize=11)
    ax4.set_title(
        f"Axial Uniformity at r = {axial_result['r'][0]:.2f} mm",
        fontsize=12,
        fontweight="bold",
    )
    ax4.legend(fontsize=9)
    ax4.grid(True, alpha=0.3)

# Plot 5: Combined error statistics
ax5 = fig.add_subplot(gs[2, :])
if radial_lines:
    angles = [line["name"] for line in radial_lines]
    max_errors = [np.max(line["abs_error"][line["mask"]]) for line in radial_lines]
    mean_errors = [np.mean(line["abs_error"][line["mask"]]) for line in radial_lines]

    x = np.arange(len(angles))
    width = 0.35

    ax5.bar(x - width / 2, max_errors, width, label="Max Error", color="indianred")
    ax5.bar(x + width / 2, mean_errors, width, label="Mean Error", color="steelblue")
    ax5.set_xlabel("Radial Line", fontsize=11)
    ax5.set_ylabel("Absolute Error (N/mm³)", fontsize=11)
    ax5.set_title("Error Summary Across All Lines", fontsize=12, fontweight="bold")
    ax5.set_xticks(x)
    ax5.set_xticklabels(angles)
    ax5.legend(fontsize=9)
    ax5.grid(True, alpha=0.3, axis="y")

plt.savefig("lorentz_verification.png", dpi=150, bbox_inches="tight")
print("\nPlot saved as 'lorentz_verification.png'")
plt.show()

# Create separate plot for 0 degree line error metrics
if radial_lines:
    # Find the 0 degree line
    line_0deg = None
    for line in radial_lines:
        if "0°" in line["name"]:
            line_0deg = line
            break

    if line_0deg is not None:
        fig2 = plt.figure(figsize=(15, 5))

        mask = line_0deg["mask"]

        # Plot 1: Comparison of simulated vs analytical
        ax1 = fig2.add_subplot(1, 3, 1)
        ax1.plot(
            line_0deg["r"][mask],
            line_0deg["lorentz_sim"][mask],
            "o",
            color="red",
            markersize=4,
            alpha=0.6,
            label="MOOSE 0°",
        )
        ax1.plot(
            line_0deg["r"][mask],
            line_0deg["lorentz_analytical"][mask],
            "k-",
            linewidth=2,
            label="Analytical",
        )
        ax1.set_xlabel("Radius (mm)", fontsize=12)
        ax1.set_ylabel("Radial Lorentz Force (N/mm³)", fontsize=12)
        ax1.set_title("0° Line: MOOSE vs Analytical", fontsize=13, fontweight="bold")
        ax1.legend(fontsize=10)
        ax1.grid(True, alpha=0.3)

        # Plot 2: Absolute error for 0 degree line
        ax2 = fig2.add_subplot(1, 3, 2)
        ax2.plot(
            line_0deg["r"][mask],
            line_0deg["abs_error"][mask],
            "o-",
            color="red",
            linewidth=2,
            markersize=4,
            label="0° line",
        )
        ax2.set_xlabel("Radius (mm)", fontsize=12)
        ax2.set_ylabel("Absolute Error (N/mm³)", fontsize=12)
        ax2.set_title("0° Line: Absolute Error", fontsize=13, fontweight="bold")
        ax2.legend(fontsize=10)
        ax2.grid(True, alpha=0.3)
        ax2.axhline(y=0, color="k", linestyle="--", linewidth=1)

        # Plot 3: Relative error for 0 degree line
        ax3 = fig2.add_subplot(1, 3, 3)
        ax3.plot(
            line_0deg["r"][mask],
            line_0deg["rel_error"][mask] * 100,
            "o-",
            color="red",
            linewidth=2,
            markersize=4,
            label="0° line",
        )
        ax3.set_xlabel("Radius (mm)", fontsize=12)
        ax3.set_ylabel("Relative Error (%)", fontsize=12)
        ax3.set_title("0° Line: Relative Error", fontsize=13, fontweight="bold")
        ax3.legend(fontsize=10)
        ax3.grid(True, alpha=0.3)
        ax3.axhline(y=0, color="k", linestyle="--", linewidth=1)

        plt.tight_layout()
        plt.savefig("lorentz_0deg_error.png", dpi=150, bbox_inches="tight")
        print("Plot saved as 'lorentz_0deg_error.png'")
        plt.show()
    else:
        print("\nWarning: 0° line data not found for error plot")
