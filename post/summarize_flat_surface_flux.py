#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from plot_steady_profiles import (
    H2O_MASS_PER_MOLECULE,
    averaging_vtr_paths,
    load_averaged_rectilinear_grid,
)


BOLTZMANN = 1.380649e-23
UNIVERSAL_GAS_CONSTANT = 8.31446261815324
H2O_MOLAR_MASS = 0.01801528
H2O_SPECIFIC_GAS_CONSTANT = UNIVERSAL_GAS_CONSTANT / H2O_MOLAR_MASS
OMEGA = 32.0 * math.pi / (32.0 + 9.0 * math.pi)
DEFAULT_TARGET_Y = 30.0e-6


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def xavg_row_values(case_dir: Path, target_y: float) -> tuple[float, float]:
    try:
        vtr_paths = averaging_vtr_paths(case_dir)
    except FileNotFoundError as exc:
        grid_dump = case_dir / "grid_steady.dump"
        if grid_dump.exists():
            raise FileNotFoundError(
                f"Missing exported VTK frames in {case_dir / 'vtk_series' / 'grid_steady'}. "
                f"grid_steady.dump exists, so rerun export_paraview_vtk.py for this case."
            ) from exc
        raise FileNotFoundError(
            f"Missing steady outputs for {case_dir}. "
            f"Expected {grid_dump} from the SPARTA run before postprocessing."
        ) from exc
    grid = load_averaged_rectilinear_grid(vtr_paths)
    y_centers = grid["y_centers"]
    y_index = int(min(range(len(y_centers)), key=lambda i: abs(y_centers[i] - target_y)))
    press_row = grid["press"][0, y_index, :]
    nrho_row = grid["nrho"][0, y_index, :]
    vy_row = grid["velocity"][0, y_index, :, 1]
    pressure = float(press_row.mean())
    local_mass_flux_y = float((-(nrho_row * H2O_MASS_PER_MOLECULE * vy_row)).mean())
    return pressure, local_mass_flux_y


def summarize_case(case_dir: Path, study_defaults: dict[str, float], target_y: float) -> dict[str, float | str]:
    metadata = load_json(case_dir / "metadata.json")
    pressure, local_mass_flux_y = xavg_row_values(case_dir, target_y)
    equilibrium_pressure = metadata["top_boundary_number_density"] * BOLTZMANN * metadata["top_boundary_temperature_k"]
    liquid_temperature = metadata["temperature_k"]
    condensation_coefficient = metadata["condensation_coefficient"]
    saturation_pressure = study_defaults["vapor_number_density"] * BOLTZMANN * liquid_temperature
    kinetic_prefactor = OMEGA * condensation_coefficient / (
        condensation_coefficient + (1.0 - condensation_coefficient) * OMEGA
    )
    reference_flux = kinetic_prefactor * (pressure - saturation_pressure) / math.sqrt(
        2.0 * math.pi * H2O_SPECIFIC_GAS_CONSTANT * liquid_temperature
    )
    normalized_flux = local_mass_flux_y / reference_flux if reference_flux else 0.0
    return {
        "case_name": metadata["case_name"],
        "equilibrium_pressure_Pa": equilibrium_pressure,
        "pressure_at_target_y_Pa": pressure,
        "local_mass_flux_y_kg_m2_s": local_mass_flux_y,
        "reference_flux_model": reference_flux,
        "normalized_local_mass_flux_y": normalized_flux,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize normalized gas-phase condensation flux for flat-surface reference cases.")
    parser.add_argument("study_dir", type=Path, help="Study directory under cases/")
    parser.add_argument("--output", default="flat_flux_summary.dat", help="Output filename within the study directory")
    parser.add_argument("--target-y", type=float, default=DEFAULT_TARGET_Y, help="Height above the flat wall [m] for x-averaged sampling")
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
            rows.append(summarize_case(case_dir, study_defaults, args.target_y))

    output_path = study_dir / args.output
    fieldnames = [
        "case_name",
        "equilibrium_pressure_Pa",
        "pressure_at_target_y_Pa",
        "local_mass_flux_y_kg_m2_s",
        "reference_flux_model",
        "normalized_local_mass_flux_y",
    ]
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(" ".join(fieldnames) + "\n")
        for row in rows:
            handle.write(" ".join(str(row[name]) for name in fieldnames) + "\n")

    print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
