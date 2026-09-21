# This file is part of the MOOSE framework
# https://mooseframework.inl.gov
#
# All rights reserved, see COPYRIGHT for full restrictions
# https://github.com/idaholab/moose/blob/master/COPYRIGHT
#
# Licensed under LGPL 2.1, please see LICENSE for details
# https://www.gnu.org/licenses/lgpl-2.1.html

"""Helper functions for parameter consistency checking
- to be integrated into verify_free_thermal_expansion.py."""


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
            f"Generate it first with generate_analytical_thermal_strain.py"
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

    # Check time array alignment
    time_moose = df["time"].values
    time_analytical = analytical_params["time"]
    if len(time_moose) == len(time_analytical):
        time_error = np.max(np.abs(time_moose - time_analytical))
        time_rel_error = time_error / (np.max(time_analytical) + 1e-12)
        time_passed = time_rel_error < rtol
        status = "✓ PASS" if time_passed else "✗ FAIL"
        print("\n  Time array alignment:")
        print(f"    Max abs error: {time_error:.2e} s")
        print(f"    Max rel error: {time_rel_error:.2e}")
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
