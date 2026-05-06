#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path


TARGET_Y = 300.0e-6


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


def pressure_at_target_y(path: Path, target_y: float) -> float:
    with path.open("r", encoding="utf-8") as handle:
        next(handle)
        best_pressure = None
        best_distance = None
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            dist_str, press_str, *_rest = stripped.split()
            distance = float(dist_str)
            pressure = float(press_str)
            delta = abs(distance - target_y)
            if best_distance is None or delta < best_distance:
                best_distance = delta
                best_pressure = pressure
    if best_pressure is None:
        raise ValueError(f"No profile rows found in {path}")
    return best_pressure


def summarize_case(case_dir: Path) -> dict[str, object]:
    metadata = load_json(case_dir / "metadata.json")
    y_axis_path = case_dir / "profiles_steady" / "y_axis.dat"
    if not y_axis_path.exists():
        raise FileNotFoundError(f"Missing profile file: {y_axis_path}")

    pressure = pressure_at_target_y(y_axis_path, TARGET_Y)
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

    flux_per_wall_length = total_mass_rate / full_domain_width if full_domain_width else 0.0
    flux_per_surface_length = total_mass_rate / total_surface_length if total_surface_length else 0.0

    return {
        "case_name": metadata["case_name"],
        "pressure_at_y_300um_Pa": pressure,
        "domain_width_m": full_domain_width,
        "mass_flux_per_wall_length": flux_per_wall_length,
        "mass_flux_per_surface_length": flux_per_surface_length,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize study pressure and condensation mass flux metrics.")
    parser.add_argument("study_dir", type=Path, help="Study directory under cases/")
    parser.add_argument("--output", default="condensation_flux_summary.dat", help="Output filename within the study directory")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    study_dir = args.study_dir.resolve()
    case_list = study_dir / "case_list.txt"
    if not case_list.exists():
        raise FileNotFoundError(f"Case list not found: {case_list}")

    rows: list[dict[str, object]] = []
    with case_list.open("r", encoding="utf-8") as handle:
        for line in handle:
            case_relpath = line.strip()
            if not case_relpath:
                continue
            case_dir = study_dir / Path(case_relpath).name
            rows.append(summarize_case(case_dir))

    output_path = study_dir / args.output
    fieldnames = [
        "case_name",
        "pressure_at_y_300um_Pa",
        "domain_width_m",
        "mass_flux_per_wall_length",
        "mass_flux_per_surface_length",
    ]
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(" ".join(fieldnames) + "\n")
        for row in rows:
            handle.write(" ".join(str(row[field]) for field in fieldnames) + "\n")

    print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
