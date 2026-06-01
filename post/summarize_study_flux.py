#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_steady_profiles import (
    H2O_MASS_PER_MOLECULE,
    build_cell_table,
    droplet_center,
    load_averaged_grid_rows,
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


def read_profile_table(path: Path) -> dict[str, np.ndarray]:
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().split()
        rows = [[float(value) for value in line.split()] for line in handle if line.strip()]
    if not rows:
        raise ValueError(f"No profile rows found in {path}")
    data = np.asarray(rows, dtype=float)
    return {name: data[:, index] for index, name in enumerate(header)}


def plateau_window_from_profile(path: Path, min_y: float, window: int, tolerance: float) -> dict[str, float]:
    profile = read_profile_table(path)
    distance = profile["dist_m"]
    pressure = profile["press_Pa"]
    candidates = np.where(distance >= min_y)[0]
    if candidates.size == 0:
        candidates = np.arange(len(distance))
    start_index = int(candidates[0])
    best: tuple[float, int, float, float] | None = None
    for index in range(start_index, len(pressure) - window + 1):
        p_window = pressure[index:index + window]
        d_window = distance[index:index + window]
        span = float(np.max(p_window) - np.min(p_window))
        mean_pressure = float(np.mean(p_window))
        mean_distance = float(np.mean(d_window))
        relative_span = span / abs(mean_pressure) if mean_pressure else span
        if best is None or relative_span < best[0]:
            best = (relative_span, index, mean_pressure, mean_distance)
    if best is None:
        raise ValueError(f"Could not find pressure plateau candidates in {path}")
    if best[0] > tolerance:
        print(
            f"Warning: best pressure plateau span {best[0]:.6g} exceeds tolerance {tolerance:.6g} in {path}",
            file=sys.stderr,
        )
    return {
        "relative_span": best[0],
        "start_index": float(best[1]),
        "pressure": best[2],
        "distance": best[3],
    }


def plateau_pressure_from_profile(path: Path, min_y: float, window: int, tolerance: float) -> float:
    return plateau_window_from_profile(path, min_y, window, tolerance)["pressure"]


def one_sided_mass_flux_to_surface(number_density: float, temperature: float, normal_velocity: float) -> float:
    """Mass flux from a drifting Maxwellian toward a surface.

    normal_velocity is positive away from the surface. Negative values increase
    the incoming flux toward the surface.
    """
    if temperature <= 0.0 or number_density <= 0.0:
        return 0.0
    cth = math.sqrt(BOLTZMANN * temperature / H2O_MASS_PER_MOLECULE)
    s = normal_velocity / (math.sqrt(2.0) * cth)
    density = number_density * H2O_MASS_PER_MOLECULE
    return density * (
        cth / math.sqrt(2.0 * math.pi) * math.exp(-(s * s))
        - 0.5 * normal_velocity * math.erfc(s)
    )


def drift_corrected_reference_flux(
    number_density: float,
    gas_temperature: float,
    normal_velocity: float,
    saturation_number_density: float,
    liquid_temperature: float,
    condensation_coefficient: float,
) -> float:
    kinetic_prefactor = OMEGA * condensation_coefficient / (
        condensation_coefficient + (1.0 - condensation_coefficient) * OMEGA
    )
    incoming_flux = one_sided_mass_flux_to_surface(number_density, gas_temperature, normal_velocity)
    saturation_flux = one_sided_mass_flux_to_surface(saturation_number_density, liquid_temperature, 0.0)
    return kinetic_prefactor * (incoming_flux - saturation_flux)


def xavg_row_values(case_dir: Path, metadata: dict, target_y: float) -> tuple[float, float, float]:
    state = xavg_row_state(case_dir, metadata, target_y)
    return state["pressure"], state["mass_flux_y"], state["integrated_mass_flux_y"]


def xavg_row_state(case_dir: Path, metadata: dict, target_y: float) -> dict[str, float]:
    rows, _ = load_averaged_grid_rows(case_dir)
    table = build_cell_table(rows, metadata)
    _, y_center = droplet_center(metadata)
    apex_y = y_center + metadata["radius"]
    target_absolute_y = apex_y + target_y
    unique_y = np.unique(table["y"])
    yline = float(unique_y[np.argmin(np.abs(unique_y - target_absolute_y))])
    row_mask = np.isclose(table["y"], yline)
    row_cells = np.sort(table[row_mask], order="x")
    row_width = float(row_cells["dx"].sum())
    if row_width == 0.0:
        raise ValueError(f"Row width is zero while sampling {case_dir}")
    pressure = float(np.sum(row_cells["press"] * row_cells["dx"]) / row_width)
    number_density = float(np.sum(row_cells["nrho"] * row_cells["dx"]) / row_width)
    temperature = float(np.sum(row_cells["temp"] * row_cells["dx"]) / row_width)
    velocity_y = float(np.sum(row_cells["v"] * row_cells["dx"]) / row_width)
    mass_flux_row = -(row_cells["nrho"] * H2O_MASS_PER_MOLECULE * row_cells["v"])
    local_mass_flux_y = float(np.sum(mass_flux_row * row_cells["dx"]) / row_width)
    symmetry_multiplier = metadata["droplets"][0]["symmetry_multiplier"]
    integrated_mass_flux_y = float(np.sum(mass_flux_row * row_cells["dx"]) * symmetry_multiplier)
    return {
        "sample_y_from_apex_m": target_y,
        "sample_y_absolute_m": target_absolute_y,
        "sample_y_grid_m": yline,
        "pressure": pressure,
        "number_density": number_density,
        "temperature": temperature,
        "velocity_y": velocity_y,
        "mass_flux_y": local_mass_flux_y,
        "integrated_mass_flux_y": integrated_mass_flux_y,
    }


def summarize_case(
    case_dir: Path,
    study_defaults: dict[str, float],
    target_y: float,
    target_flux_y: float,
    sample_mode: str,
    flux_source: str,
    pressure_mode: str,
    plateau_min_y: float,
    plateau_window: int,
    plateau_tolerance: float,
    include_drift_reference: bool,
) -> dict[str, object]:
    metadata = load_json(case_dir / "metadata.json")
    y_axis_path = case_dir / "profiles_steady" / "y_axis.dat"
    if not y_axis_path.exists():
        raise FileNotFoundError(f"Missing profile file: {y_axis_path}")

    if pressure_mode == "plateau":
        plateau = plateau_window_from_profile(y_axis_path, plateau_min_y, plateau_window, plateau_tolerance)
        pressure = plateau["pressure"]
    elif pressure_mode == "target_y":
        plateau = None
        pressure = None
    else:
        raise ValueError(f"Unsupported pressure mode: {pressure_mode}")

    if sample_mode == "centerline":
        if pressure is None:
            pressure = nearest_profile_value(y_axis_path, target_y, 1)
        local_mass_flux_y = nearest_profile_value(y_axis_path, target_flux_y, 4)
        symmetry_multiplier = metadata["droplets"][0]["symmetry_multiplier"]
        integrated_mass_flux_y = local_mass_flux_y * (
            metadata["simulation_bounds"]["xhi"] - metadata["simulation_bounds"]["xlo"]
        ) * symmetry_multiplier
    elif sample_mode == "xavg":
        if pressure is None:
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
    drift_state = {
        "pressure": 0.0,
        "number_density": 0.0,
        "temperature": 0.0,
        "velocity_y": 0.0,
        "mass_flux_y": 0.0,
    }
    drift_reference_flux = 0.0
    drift_normalized_flux_per_wall_length = 0.0
    drift_normalized_flux_per_surface_length = 0.0
    drift_normalized_local_mass_flux_y = 0.0
    if include_drift_reference:
        if pressure_mode == "plateau":
            drift_sample_y = plateau["distance"] if plateau is not None else target_y
        else:
            drift_sample_y = target_y
        drift_state = xavg_row_state(case_dir, metadata, drift_sample_y)
        drift_reference_flux = drift_corrected_reference_flux(
            drift_state["number_density"],
            drift_state["temperature"],
            drift_state["velocity_y"],
            vapor_number_density,
            liquid_temperature,
            condensation_coefficient,
        )
        drift_normalized_flux_per_wall_length = (
            flux_per_wall_length / drift_reference_flux if drift_reference_flux else 0.0
        )
        drift_normalized_flux_per_surface_length = (
            flux_per_surface_length / drift_reference_flux if drift_reference_flux else 0.0
        )
        drift_normalized_local_mass_flux_y = (
            local_mass_flux_y / drift_reference_flux if drift_reference_flux else 0.0
        )

    return {
        "case_name": metadata["case_name"],
        "equilibrium_pressure_Pa": equilibrium_pressure,
        "xhi_m": xhi,
        "reference_pressure_Pa": pressure,
        "local_mass_flux_y_at_target": local_mass_flux_y,
        "reference_flux_model": reference_flux,
        "domain_width_m": full_domain_width,
        "mass_flux_per_wall_length": flux_per_wall_length,
        "mass_flux_per_surface_length": flux_per_surface_length,
        "normalized_flux_per_wall_length": normalized_flux_per_wall_length,
        "normalized_flux_per_surface_length": normalized_flux_per_surface_length,
        "normalized_local_mass_flux_y_at_target": normalized_local_mass_flux_y,
        "drift_reference_pressure_Pa": drift_state["pressure"],
        "drift_reference_number_density_m3": drift_state["number_density"],
        "drift_reference_temperature_K": drift_state["temperature"],
        "drift_reference_velocity_y_m_s": drift_state["velocity_y"],
        "drift_reference_mass_flux_y_kg_m2_s": drift_state["mass_flux_y"],
        "drift_reference_flux_model": drift_reference_flux,
        "drift_normalized_flux_per_wall_length": drift_normalized_flux_per_wall_length,
        "drift_normalized_flux_per_surface_length": drift_normalized_flux_per_surface_length,
        "drift_normalized_local_mass_flux_y_at_target": drift_normalized_local_mass_flux_y,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize study pressure and condensation mass flux metrics.")
    parser.add_argument("study_dir", type=Path, help="Study directory under cases/")
    parser.add_argument("--output", default="condensation_flux_summary.dat", help="Output filename within the study directory")
    parser.add_argument("--target-y", type=float, default=TARGET_Y, help="Sample pressure from y_axis.dat at this distance from the droplet surface [m]")
    parser.add_argument("--target-flux-y", type=float, default=TARGET_Y, help="Sample gas-phase mass flux from y_axis.dat at this distance from the droplet surface [m]")
    parser.add_argument("--sample-mode", choices=("centerline", "xavg"), default="centerline", help="How to sample gas-phase pressure and mass flux from the averaged field")
    parser.add_argument("--flux-source", choices=("surf_dump", "gas_integral"), default="surf_dump", help="Source used to build wall/surface normalized flux columns")
    parser.add_argument("--pressure-mode", choices=("target_y", "plateau"), default="target_y", help="Pressure used in the reference flux model")
    parser.add_argument("--plateau-min-y", type=float, default=50.0e-6, help="Ignore profile points below this distance when detecting pressure plateau [m]")
    parser.add_argument("--plateau-window", type=int, default=20, help="Number of consecutive y_axis points used for plateau detection")
    parser.add_argument("--plateau-tolerance", type=float, default=2.0e-3, help="Relative pressure span accepted as a plateau")
    parser.add_argument("--include-drift-reference", action="store_true", help="Append normalization using the plateau drifting-Maxwellian gas state from grid_steady.dump")
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
                    args.pressure_mode,
                    args.plateau_min_y,
                    args.plateau_window,
                    args.plateau_tolerance,
                    args.include_drift_reference,
                )
            )

    output_path = study_dir / args.output
    fieldnames = [
        "case_name",
        "equilibrium_pressure_Pa",
        "xhi_m",
        "reference_pressure_Pa",
        "local_mass_flux_y_at_target",
        "reference_flux_model",
        "domain_width_m",
        "mass_flux_per_wall_length",
        "mass_flux_per_surface_length",
        "normalized_flux_per_wall_length",
        "normalized_flux_per_surface_length",
        "normalized_local_mass_flux_y_at_target",
        "drift_reference_pressure_Pa",
        "drift_reference_number_density_m3",
        "drift_reference_temperature_K",
        "drift_reference_velocity_y_m_s",
        "drift_reference_mass_flux_y_kg_m2_s",
        "drift_reference_flux_model",
        "drift_normalized_flux_per_wall_length",
        "drift_normalized_flux_per_surface_length",
        "drift_normalized_local_mass_flux_y_at_target",
    ]
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(" ".join(fieldnames) + "\n")
        for row in rows:
            handle.write(" ".join(str(row[field]) for field in fieldnames) + "\n")

    print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
