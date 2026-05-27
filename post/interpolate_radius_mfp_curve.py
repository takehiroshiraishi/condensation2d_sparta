#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path


BOLTZMANN = 1.380649e-23
STUDY_ROOT = Path(__file__).resolve().parents[1]
WATER_VSS_PATH = STUDY_ROOT / "base" / "water.vss"


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def parse_radius(case_name: str) -> float:
    return float(case_name.split("_r_")[1].replace("p", "."))


def load_vss_diameter(path: Path) -> float:
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if parts[0] == "H2O":
                return float(parts[1])
    raise ValueError(f"Could not find H2O VSS diameter in {path}")


def mean_free_path(surface_temperature: float, pressure: float, diameter: float) -> float:
    return BOLTZMANN * surface_temperature / (math.sqrt(2.0) * math.pi * diameter * diameter * pressure)


def linear_interpolate_or_extrapolate(points: list[tuple[float, float]], target_x: float) -> tuple[float, float, float, str]:
    points = sorted(points)
    for x, y in points:
        if math.isclose(x, target_x, rel_tol=0.0, abs_tol=1.0e-12):
            return y, x, x, "exact"
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x0 <= target_x <= x1:
            t = (target_x - x0) / (x1 - x0)
            return y0 + t * (y1 - y0), x0, x1, "interpolated"
    if target_x < points[0][0]:
        x0, y0 = points[0]
        x1, y1 = points[1]
        t = (target_x - x0) / (x1 - x0)
        return y0 + t * (y1 - y0), x0, x1, "extrapolated_low"
    x0, y0 = points[-2]
    x1, y1 = points[-1]
    t = (target_x - x0) / (x1 - x0)
    return y0 + t * (y1 - y0), x0, x1, "extrapolated_high"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Interpolate a normalized condensation metric at target pressures and plot vs radius.")
    parser.add_argument("study_dirs", nargs="+", type=Path, help="Study directories under cases/")
    parser.add_argument("--target-pressure", type=float, action="append", required=True, help="Target pressure [Pa] for interpolation; pass multiple times for multiple curves")
    parser.add_argument("--output", type=Path, required=True, help="Output .dat file")
    parser.add_argument(
        "--summary-column",
        choices=("normalized_flux_per_wall_length", "normalized_flux_per_surface_length"),
        default="normalized_flux_per_surface_length",
        help="Normalized summary column to interpolate",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_path = args.output.resolve()
    diameter = load_vss_diameter(WATER_VSS_PATH)

    rows_out: list[dict[str, object]] = []
    surface_temperature: float | None = None

    for study_dir_arg in args.study_dirs:
        study_dir = study_dir_arg.resolve()
        summary_path = study_dir / "condensation_flux_summary_surf_dump.dat"
        if not summary_path.exists():
            raise FileNotFoundError(f"Missing summary file: {summary_path}")
        parameters = load_json(study_dir / "parameters.json")
        study_temperature = float(parameters["defaults"]["temperature_k"])
        if surface_temperature is None:
            surface_temperature = study_temperature
        elif not math.isclose(surface_temperature, study_temperature, rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError("All studies must use the same surface temperature for a single combined radius/lambda curve.")

        with summary_path.open("r", encoding="utf-8") as handle:
            header = handle.readline().split()
            data = [dict(zip(header, line.split())) for line in handle if line.strip()]

        grouped: dict[float, list[tuple[float, float]]] = defaultdict(list)
        for row in data:
            radius = parse_radius(row["case_name"])
            pressure_name = "reference_pressure_Pa" if "reference_pressure_Pa" in row else "pressure_at_y_300um_Pa"
            pressure = float(row[pressure_name])
            normalized_flux = float(row[args.summary_column])
            grouped[radius].append((pressure, normalized_flux))

        for target_pressure in args.target_pressure:
            for radius in sorted(grouped):
                interpolated, p_low, p_high, mode = linear_interpolate_or_extrapolate(grouped[radius], target_pressure)
                rows_out.append(
                    {
                        "study_name": study_dir.name,
                        "radius_m": radius,
                        "radius_um": radius * 1.0e6,
                        "target_pressure_Pa": target_pressure,
                        "normalized_flux": interpolated,
                        "pressure_bracket_low_Pa": p_low,
                        "pressure_bracket_high_Pa": p_high,
                        "interp_mode": mode,
                    }
                )

    rows_out.sort(key=lambda row: (float(row["target_pressure_Pa"]), float(row["radius_m"])))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "study_name",
        "radius_m",
        "radius_um",
        "target_pressure_Pa",
        "normalized_flux",
        "pressure_bracket_low_Pa",
        "pressure_bracket_high_Pa",
        "interp_mode",
    ]
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(" ".join(fieldnames) + "\n")
        for row in rows_out:
            handle.write(" ".join(str(row[name]) for name in fieldnames) + "\n")

    print(f"Wrote {output_path}")
    print(f"Surface temperature [K]: {surface_temperature}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
