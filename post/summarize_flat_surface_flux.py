#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from plot_steady_profiles import (
    H2O_MASS_PER_MOLECULE,
    build_cell_table,
    load_averaged_grid_rows,
)
from summarize_study_flux import drift_corrected_reference_flux


BOLTZMANN = 1.380649e-23
UNIVERSAL_GAS_CONSTANT = 8.31446261815324
H2O_MOLAR_MASS = 0.01801528
H2O_SPECIFIC_GAS_CONSTANT = UNIVERSAL_GAS_CONSTANT / H2O_MOLAR_MASS
OMEGA = 32.0 * math.pi / (32.0 + 9.0 * math.pi)
DEFAULT_TARGET_Y = 90.0e-6


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def xavg_row_values(case_dir: Path, target_y: float) -> tuple[float, float]:
    metadata = load_json(case_dir / "metadata.json")
    state = xavg_row_state_from_dump(case_dir, metadata, target_y)
    return state["pressure"], state["mass_flux_y"]


def xavg_row_state_from_dump(case_dir: Path, metadata: dict, target_y: float) -> dict[str, float]:
    rows, _ = load_averaged_grid_rows(case_dir)
    table = build_cell_table(rows, metadata)
    unique_y = np.unique(table["y"])
    yline = float(unique_y[np.argmin(np.abs(unique_y - target_y))])
    row_cells = np.sort(table[np.isclose(table["y"], yline)], order="x")
    row_width = float(row_cells["dx"].sum())
    if row_width == 0.0:
        raise ValueError(f"Row width is zero while sampling {case_dir}")
    mass_flux_row = -(row_cells["nrho"] * H2O_MASS_PER_MOLECULE * row_cells["v"])
    return {
        "sample_y_m": target_y,
        "sample_y_grid_m": yline,
        "pressure": float(np.sum(row_cells["press"] * row_cells["dx"]) / row_width),
        "number_density": float(np.sum(row_cells["nrho"] * row_cells["dx"]) / row_width),
        "temperature": float(np.sum(row_cells["temp"] * row_cells["dx"]) / row_width),
        "velocity_y": float(np.sum(row_cells["v"] * row_cells["dx"]) / row_width),
        "mass_flux_y": float(np.sum(mass_flux_row * row_cells["dx"]) / row_width),
    }


def summarize_case(
    case_dir: Path,
    study_defaults: dict[str, float],
    pressure_target_y: float,
    flux_target_y: float,
) -> dict[str, float | str]:
    metadata = load_json(case_dir / "metadata.json")
    pressure, _ = xavg_row_values(case_dir, pressure_target_y)
    _, local_mass_flux_y = xavg_row_values(case_dir, flux_target_y)
    drift_state = xavg_row_state_from_dump(case_dir, metadata, pressure_target_y)
    equilibrium_pressure = metadata["top_boundary_number_density"] * BOLTZMANN * metadata["top_boundary_temperature_k"]
    liquid_temperature = metadata["temperature_k"]
    condensation_coefficient = metadata["condensation_coefficient"]
    saturation_pressure = study_defaults["vapor_number_density"] * BOLTZMANN * liquid_temperature
    omega_prefactor = OMEGA * condensation_coefficient / (
        condensation_coefficient + (1.0 - condensation_coefficient) * OMEGA
    )
    schrage_prefactor = 2.0 * condensation_coefficient / (2.0 - condensation_coefficient)
    reference_flux_omega = omega_prefactor * (pressure - saturation_pressure) / math.sqrt(
        2.0 * math.pi * H2O_SPECIFIC_GAS_CONSTANT * liquid_temperature
    )
    reference_flux_schrage = schrage_prefactor * (pressure - saturation_pressure) / math.sqrt(
        2.0 * math.pi * H2O_SPECIFIC_GAS_CONSTANT * liquid_temperature
    )
    normalized_flux_omega = local_mass_flux_y / reference_flux_omega if reference_flux_omega else 0.0
    normalized_flux_schrage = local_mass_flux_y / reference_flux_schrage if reference_flux_schrage else 0.0
    drift_reference_flux = drift_corrected_reference_flux(
        drift_state["number_density"],
        drift_state["temperature"],
        drift_state["velocity_y"],
        study_defaults["vapor_number_density"],
        liquid_temperature,
        condensation_coefficient,
    )
    drift_normalized_flux = local_mass_flux_y / drift_reference_flux if drift_reference_flux else 0.0
    return {
        "case_name": metadata["case_name"],
        "equilibrium_pressure_Pa": equilibrium_pressure,
        "pressure_at_target_y_Pa": pressure,
        "local_mass_flux_y_kg_m2_s": local_mass_flux_y,
        "reference_flux_omega": reference_flux_omega,
        "reference_flux_schrage": reference_flux_schrage,
        "normalized_local_mass_flux_y_omega": normalized_flux_omega,
        "normalized_local_mass_flux_y_schrage": normalized_flux_schrage,
        "drift_reference_pressure_Pa": drift_state["pressure"],
        "drift_reference_number_density_m3": drift_state["number_density"],
        "drift_reference_temperature_K": drift_state["temperature"],
        "drift_reference_velocity_y_m_s": drift_state["velocity_y"],
        "drift_reference_mass_flux_y_kg_m2_s": drift_state["mass_flux_y"],
        "drift_reference_flux_model": drift_reference_flux,
        "drift_normalized_local_mass_flux_y": drift_normalized_flux,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize normalized gas-phase condensation flux for flat-surface reference cases.")
    parser.add_argument("study_dir", type=Path, help="Study directory under cases/")
    parser.add_argument("--output", default="flat_flux_summary.dat", help="Output filename within the study directory")
    parser.add_argument("--target-y", type=float, default=90.0e-6, help="Height above the flat wall [m] for x-averaged pressure sampling")
    parser.add_argument("--target-flux-y", type=float, default=DEFAULT_TARGET_Y, help="Height above the flat wall [m] for x-averaged mass-flux sampling")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    study_dir = args.study_dir.resolve()
    parameters = load_json(study_dir / "parameters.json")
    study_defaults = parameters["defaults"]
    case_list = study_dir / "case_list.txt"
    if not case_list.exists():
        raise FileNotFoundError(f"Case list not found: {case_list}")

    rows: list[dict[str, float | str]] = []
    with case_list.open("r", encoding="utf-8") as handle:
        for line in handle:
            case_relpath = line.strip()
            if not case_relpath:
                continue
            case_dir = study_dir / Path(case_relpath).name
            rows.append(summarize_case(case_dir, study_defaults, args.target_y, args.target_flux_y))

    output_path = study_dir / args.output
    fieldnames = [
        "case_name",
        "equilibrium_pressure_Pa",
        "pressure_at_target_y_Pa",
        "local_mass_flux_y_kg_m2_s",
        "reference_flux_omega",
        "reference_flux_schrage",
        "normalized_local_mass_flux_y_omega",
        "normalized_local_mass_flux_y_schrage",
        "drift_reference_pressure_Pa",
        "drift_reference_number_density_m3",
        "drift_reference_temperature_K",
        "drift_reference_velocity_y_m_s",
        "drift_reference_mass_flux_y_kg_m2_s",
        "drift_reference_flux_model",
        "drift_normalized_local_mass_flux_y",
    ]
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(" ".join(fieldnames) + "\n")
        for row in rows:
            handle.write(" ".join(str(row[name]) for name in fieldnames) + "\n")

    print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
