#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from plot_steady_profiles import (
    H2O_MASS_PER_MOLECULE,
    averaging_vtr_paths,
    droplet_center,
    load_averaged_rectilinear_grid,
)


TARGET_Y = 300.0e-6
BOLTZMANN = 1.380649e-23
UNIVERSAL_GAS_CONSTANT = 8.31446261815324
H2O_MOLAR_MASS = 0.01801528
H2O_SPECIFIC_GAS_CONSTANT = UNIVERSAL_GAS_CONSTANT / H2O_MOLAR_MASS
OMEGA = 32.0 * math.pi / (32.0 + 9.0 * math.pi)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def parse_surface_dump_last_frame(path: Path) -> tuple[list[str], list[dict[str, float]]]:
    last_columns: list[str] | None = None
    last_rows: list[dict[str, float]] | None = None

    with path.open("r", encoding="utf-8") as handle:
        lines = iter(handle)
        for line in lines:
            if line.strip() != "ITEM: TIMESTEP":
                continue
            next(lines)
            if next(lines).strip() != "ITEM: NUMBER OF SURFS":
                raise ValueError(f"Unexpected dump format in {path}")
            count = int(next(lines).strip())
            bounds_header = next(lines).strip()
            if not bounds_header.startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"Unexpected bounds header in {path}")
            next(lines)
            next(lines)
            next(lines)
            columns = next(lines).strip().split()[2:]
            rows = []
            for _ in range(count):
                raw = next(lines).split()
                rows.append({column: float(value) for column, value in zip(columns, raw)})
            last_columns = columns
            last_rows = rows

    if last_columns is None or last_rows is None:
        raise ValueError(f"No frames found in {path}")
    return last_columns, last_rows


def nearest_profile_value(path: Path, target_y: float, column_index: int) -> float:
    with path.open("r", encoding="utf-8") as handle:
        next(handle)
        best_value = None
        best_distance = None
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            columns = stripped.split()
            dist_str = columns[0]
            distance = float(dist_str)
            value = float(columns[column_index])
            delta = abs(distance - target_y)
            if best_distance is None or delta < best_distance:
                best_distance = delta
                best_value = value
    if best_value is None:
        raise ValueError(f"No profile rows found in {path}")
    return best_value


def xavg_row_values(case_dir: Path, metadata: dict, target_y: float) -> tuple[float, float, float]:
    vtr_paths = averaging_vtr_paths(case_dir)
    grid = load_averaged_rectilinear_grid(vtr_paths)
    _, y_center = droplet_center(metadata)
    apex_y = y_center + metadata["radius"]
    target_absolute_y = apex_y + target_y
    y_index = int(min(range(len(grid["y_centers"])), key=lambda i: abs(grid["y_centers"][i] - target_absolute_y)))

    press_row = grid["press"][0, y_index, :]
    nrho_row = grid["nrho"][0, y_index, :]
    vy_row = grid["velocity"][0, y_index, :, 1]
    x_centers = grid["x_centers"]
    if len(x_centers) > 1:
        dx = float(x_centers[1] - x_centers[0])
    else:
        dx = float(metadata["simulation_bounds"]["xhi"] - metadata["simulation_bounds"]["xlo"])
    pressure = float(press_row.mean())
    mass_flux_row = -(nrho_row * H2O_MASS_PER_MOLECULE * vy_row)
    local_mass_flux_y = float(mass_flux_row.mean())
    symmetry_multiplier = metadata["droplets"][0]["symmetry_multiplier"]
    integrated_mass_flux_y = float(mass_flux_row.sum() * dx * symmetry_multiplier)
    return pressure, local_mass_flux_y, integrated_mass_flux_y


def summarize_case(
    case_dir: Path,
    study_defaults: dict[str, float],
    target_y: float,
    target_flux_y: float,
    sample_mode: str,
    flux_source: str,
) -> dict[str, object]:
    metadata = load_json(case_dir / "metadata.json")
    y_axis_path = case_dir / "profiles_steady" / "y_axis.dat"
    if not y_axis_path.exists():
        raise FileNotFoundError(f"Missing profile file: {y_axis_path}")

    if sample_mode == "centerline":
        pressure = nearest_profile_value(y_axis_path, target_y, 1)
        local_mass_flux_y = nearest_profile_value(y_axis_path, target_flux_y, 4)
        symmetry_multiplier = metadata["droplets"][0]["symmetry_multiplier"]
        integrated_mass_flux_y = local_mass_flux_y * (
            metadata["simulation_bounds"]["xhi"] - metadata["simulation_bounds"]["xlo"]
        ) * symmetry_multiplier
    elif sample_mode == "xavg":
        pressure, _, _ = xavg_row_values(case_dir, metadata, target_y)
        _, local_mass_flux_y, integrated_mass_flux_y = xavg_row_values(case_dir, metadata, target_flux_y)
    else:
        raise ValueError(f"Unsupported sample mode: {sample_mode}")
    equilibrium_pressure = metadata["top_boundary_number_density"] * BOLTZMANN * metadata["top_boundary_temperature_k"]
    liquid_temperature = metadata["temperature_k"]
    condensation_coefficient = metadata["condensation_coefficient"]
    vapor_number_density = study_defaults["vapor_number_density"]
    saturation_pressure = vapor_number_density * BOLTZMANN * liquid_temperature
    xhi = metadata["simulation_bounds"]["xhi"]
    full_domain_width = 2.0 * metadata["simulation_bounds"]["xhi"]

    total_mass_rate = 0.0
    total_surface_length = 0.0
    for droplet_index, droplet_meta in enumerate(metadata["droplets"], start=1):
        dump_path = case_dir / f"surf_droplet{droplet_index}.dump"
        if not dump_path.exists():
            raise FileNotFoundError(f"Missing surface dump: {dump_path}")
        columns, rows = parse_surface_dump_last_frame(dump_path)
        flow_column = f"f_avg_droplet_{droplet_index}[2]"
        if flow_column not in columns:
            raise KeyError(f"Missing {flow_column} in {dump_path}")
        total_mass_rate += sum(row[flow_column] for row in rows) * droplet_meta["symmetry_multiplier"]
        total_surface_length += droplet_meta["analytic_arc_length_full_m"]

    if flux_source == "surf_dump":
        flux_per_wall_length = total_mass_rate / full_domain_width if full_domain_width else 0.0
        flux_per_surface_length = total_mass_rate / total_surface_length if total_surface_length else 0.0
    elif flux_source == "gas_integral":
        flux_per_wall_length = integrated_mass_flux_y / full_domain_width if full_domain_width else 0.0
        flux_per_surface_length = integrated_mass_flux_y / total_surface_length if total_surface_length else 0.0
    else:
        raise ValueError(f"Unsupported flux source: {flux_source}")
    kinetic_prefactor = OMEGA * condensation_coefficient / (
        condensation_coefficient + (1.0 - condensation_coefficient) * OMEGA
    )
    reference_flux = kinetic_prefactor * (pressure - saturation_pressure) / math.sqrt(
        2.0 * math.pi * H2O_SPECIFIC_GAS_CONSTANT * liquid_temperature
    )
    normalized_flux_per_wall_length = flux_per_wall_length / reference_flux if reference_flux else 0.0
    normalized_flux_per_surface_length = flux_per_surface_length / reference_flux if reference_flux else 0.0
    normalized_local_mass_flux_y = local_mass_flux_y / reference_flux if reference_flux else 0.0

    return {
        "case_name": metadata["case_name"],
        "equilibrium_pressure_Pa": equilibrium_pressure,
        "xhi_m": xhi,
        "pressure_at_y_300um_Pa": pressure,
        "local_mass_flux_y_at_target": local_mass_flux_y,
        "reference_flux_model": reference_flux,
        "domain_width_m": full_domain_width,
        "mass_flux_per_wall_length": flux_per_wall_length,
        "mass_flux_per_surface_length": flux_per_surface_length,
        "normalized_flux_per_wall_length": normalized_flux_per_wall_length,
        "normalized_flux_per_surface_length": normalized_flux_per_surface_length,
        "normalized_local_mass_flux_y_at_target": normalized_local_mass_flux_y,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize study pressure and condensation mass flux metrics.")
    parser.add_argument("study_dir", type=Path, help="Study directory under cases/")
    parser.add_argument("--output", default="condensation_flux_summary.dat", help="Output filename within the study directory")
    parser.add_argument("--target-y", type=float, default=TARGET_Y, help="Sample pressure from y_axis.dat at this distance from the droplet surface [m]")
    parser.add_argument("--target-flux-y", type=float, default=TARGET_Y, help="Sample gas-phase mass flux from y_axis.dat at this distance from the droplet surface [m]")
    parser.add_argument("--sample-mode", choices=("centerline", "xavg"), default="centerline", help="How to sample gas-phase pressure and mass flux from the averaged field")
    parser.add_argument("--flux-source", choices=("surf_dump", "gas_integral"), default="surf_dump", help="Source used to build wall/surface normalized flux columns")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    study_dir = args.study_dir.resolve()
    case_list = study_dir / "case_list.txt"
    if not case_list.exists():
        raise FileNotFoundError(f"Case list not found: {case_list}")
    parameters_path = study_dir / "parameters.json"
    if not parameters_path.exists():
        raise FileNotFoundError(f"Study parameters not found: {parameters_path}")
    parameters = load_json(parameters_path)
    study_defaults = parameters["defaults"]
    if "vapor_number_density" not in study_defaults:
        raise KeyError(f"Missing defaults.vapor_number_density in {parameters_path}")

    rows: list[dict[str, object]] = []
    with case_list.open("r", encoding="utf-8") as handle:
        for line in handle:
            case_relpath = line.strip()
            if not case_relpath:
                continue
            case_dir = study_dir / Path(case_relpath).name
            rows.append(
                summarize_case(
                    case_dir,
                    study_defaults,
                    args.target_y,
                    args.target_flux_y,
                    args.sample_mode,
                    args.flux_source,
                )
            )

    output_path = study_dir / args.output
    fieldnames = [
        "case_name",
        "equilibrium_pressure_Pa",
        "xhi_m",
        "pressure_at_y_300um_Pa",
        "local_mass_flux_y_at_target",
        "reference_flux_model",
        "domain_width_m",
        "mass_flux_per_wall_length",
        "mass_flux_per_surface_length",
        "normalized_flux_per_wall_length",
        "normalized_flux_per_surface_length",
        "normalized_local_mass_flux_y_at_target",
    ]
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(" ".join(fieldnames) + "\n")
        for row in rows:
            handle.write(" ".join(str(row[field]) for field in fieldnames) + "\n")

    print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
