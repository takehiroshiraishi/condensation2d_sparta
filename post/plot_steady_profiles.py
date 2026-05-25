#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

AVOGADRO = 6.02214076e23
H2O_MOLAR_MASS = 0.01801528
H2O_MASS_PER_MOLECULE = H2O_MOLAR_MASS / AVOGADRO


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def parse_grid_dump_frames(path: Path) -> list[dict]:
    frames: list[dict] = []
    with path.open("r", encoding="utf-8") as handle:
        lines = iter(handle)
        for line in lines:
            if line.strip() != "ITEM: TIMESTEP":
                continue
            timestep = int(next(lines).strip())
            if next(lines).strip() != "ITEM: NUMBER OF CELLS":
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
                row = {column: float(value) for column, value in zip(columns, raw)}
                if "zc" not in row:
                    row["zc"] = 0.0
                rows.append(row)
            frames.append({"timestep": timestep, "rows": rows})
    return frames


def select_averaging_frames(frames: list[dict]) -> list[dict]:
    if not frames:
        raise ValueError("No frames available to average.")
    nonzero_frames = [frame for frame in frames if frame["timestep"] != 0]
    if not nonzero_frames:
        return frames[-1:]
    if len(nonzero_frames) == 1:
        return nonzero_frames
    return nonzero_frames[1:]


def map_frame_fields(rows: list[dict]) -> list[dict]:
    if not rows:
        return []
    source_fields = [
        name
        for name in ("f_grid_avg[1]", "f_grid_avg[2]", "f_grid_avg[3]", "f_grid_avg[4]", "f_grid_avg[5]", "f_grid_avg[6]")
        if name in rows[0]
    ]
    if not source_fields:
        source_fields = [
            name
            for name in ("f_centerline_avg[1]", "f_centerline_avg[2]", "f_centerline_avg[3]", "f_centerline_avg[4]", "f_centerline_avg[5]", "f_centerline_avg[6]")
            if name in rows[0]
        ]
    target_fields = ("nrho", "u", "v", "trot", "temp", "press")
    mapped_rows: list[dict] = []
    for row in rows:
        mapped = dict(row)
        for source_name, target_name in zip(source_fields, target_fields):
            mapped[target_name] = row[source_name]
        mapped.setdefault("u", 0.0)
        mapped.setdefault("v", 0.0)
        mapped.setdefault("trot", 0.0)
        mapped.setdefault("temp", 0.0)
        mapped.setdefault("press", 0.0)
        mapped_rows.append(mapped)
    return mapped_rows


def average_frame_rows(frames: list[dict]) -> list[dict]:
    selected = select_averaging_frames(frames)
    reference_rows = map_frame_fields(selected[0]["rows"])
    sum_rows = {int(row["id"]): dict(row) for row in reference_rows}
    for row in sum_rows.values():
        for field in ("nrho", "u", "v", "trot", "temp", "press"):
            row[field] = float(row.get(field, 0.0))

    for frame in selected[1:]:
        for row in map_frame_fields(frame["rows"]):
            summed = sum_rows[int(row["id"])]
            for field in ("nrho", "u", "v", "trot", "temp", "press"):
                summed[field] += float(row.get(field, 0.0))

    count = float(len(selected))
    averaged_rows = []
    for row_id in sorted(sum_rows):
        row = sum_rows[row_id]
        for field in ("nrho", "u", "v", "trot", "temp", "press"):
            row[field] /= count
        averaged_rows.append(row)
    return averaged_rows


def load_averaged_grid_rows(case_dir: Path, dump_name: str = "grid_steady.dump") -> tuple[list[dict], list[int]]:
    dump_path = case_dir / dump_name
    frames = parse_grid_dump_frames(dump_path)
    selected = select_averaging_frames(frames)
    return average_frame_rows(frames), [frame["timestep"] for frame in selected]


def droplet_center(metadata: dict) -> tuple[float, float]:
    x_center = metadata["droplets"][0]["center_x"]
    radius = metadata["radius"]
    theta_deg = metadata["contact_angle_deg"]
    surface_gap = metadata["surface_gap"]
    y_center = surface_gap - radius * math.cos(math.radians(theta_deg))
    return x_center, y_center


def surface_point(x0: float, y0: float, radius: float, direction_x: float, direction_y: float) -> tuple[float, float]:
    norm = math.hypot(direction_x, direction_y)
    if norm == 0.0:
        raise ValueError("Surface direction must be non-zero.")
    return x0 + radius * direction_x / norm, y0 + radius * direction_y / norm


def reconstruct_edges_from_centers(centers: list[float], lower: float, upper: float) -> list[float]:
    if not centers:
        return []
    if len(centers) == 1:
        return [lower, upper]
    boundary = lower
    boundaries = [boundary]
    for center in centers:
        boundary = 2.0 * center - boundary
        boundaries.append(boundary)
    correction = upper - boundaries[-1]
    adjusted = [edge + (correction if index % 2 == 1 else 0.0) for index, edge in enumerate(boundaries)]
    adjusted[0] = lower
    adjusted[-1] = upper
    return adjusted


def build_cell_table(rows: list[dict], metadata: dict) -> np.ndarray:
    bounds = metadata["simulation_bounds"]
    rows_by_y: dict[float, list[dict]] = {}
    rows_by_x: dict[float, list[dict]] = {}
    for row in rows:
        rows_by_y.setdefault(row["yc"], []).append(row)
        rows_by_x.setdefault(row["xc"], []).append(row)

    x_bounds_by_cell: dict[int, tuple[float, float]] = {}
    for y_center, row_group in rows_by_y.items():
        del y_center
        ordered = sorted(row_group, key=lambda item: item["xc"])
        centers = [row["xc"] for row in ordered]
        edges = reconstruct_edges_from_centers(centers, bounds["xlo"], bounds["xhi"])
        for index, row in enumerate(ordered):
            x_bounds_by_cell[int(row["id"])] = (edges[index], edges[index + 1])

    y_bounds_by_cell: dict[int, tuple[float, float]] = {}
    for x_center, row_group in rows_by_x.items():
        del x_center
        ordered = sorted(row_group, key=lambda item: item["yc"])
        centers = [row["yc"] for row in ordered]
        edges = reconstruct_edges_from_centers(centers, bounds["ylo"], bounds["yhi"])
        for index, row in enumerate(ordered):
            y_bounds_by_cell[int(row["id"])] = (edges[index], edges[index + 1])

    dtype = [
        ("id", int),
        ("x", float),
        ("y", float),
        ("xlo", float),
        ("xhi", float),
        ("ylo", float),
        ("yhi", float),
        ("dx", float),
        ("dy", float),
        ("nrho", float),
        ("u", float),
        ("v", float),
        ("trot", float),
        ("temp", float),
        ("press", float),
    ]
    table = np.empty(len(rows), dtype=dtype)
    for index, row in enumerate(sorted(rows, key=lambda item: int(item["id"]))):
        row_id = int(row["id"])
        xlo, xhi = x_bounds_by_cell[row_id]
        ylo, yhi = y_bounds_by_cell[row_id]
        table[index] = (
            row_id,
            row["xc"],
            row["yc"],
            xlo,
            xhi,
            ylo,
            yhi,
            xhi - xlo,
            yhi - ylo,
            row.get("nrho", 0.0),
            row.get("u", 0.0),
            row.get("v", 0.0),
            row.get("trot", 0.0),
            row.get("temp", 0.0),
            row.get("press", 0.0),
        )
    return table


def positive_spacing(values: np.ndarray) -> float:
    unique = np.unique(values)
    if len(unique) < 2:
        return 1.0
    diffs = np.diff(np.sort(unique))
    positive = diffs[diffs > 0.0]
    if positive.size == 0:
        return 1.0
    return float(positive.min())


def select_axis_x(table: np.ndarray, x0: float, y0: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    yline = float(table["y"][np.argmin(np.abs(table["y"] - y0))])
    x_surface = x0
    dx_tol = positive_spacing(table["x"])
    mask = np.isclose(table["y"], yline) & (table["x"] >= x_surface - 0.5 * dx_tol)
    rows = np.sort(table[mask], order="x")
    s = rows["x"] - x_surface
    mass_flux_y = -(rows["nrho"] * H2O_MASS_PER_MOLECULE * rows["v"])
    return s, rows["trot"], rows["temp"], rows["press"], mass_flux_y


def select_axis_y(table: np.ndarray, x0: float, y0: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    xline = float(table["x"][np.argmin(np.abs(table["x"] - x0))])
    y_surface = y0
    dy_tol = positive_spacing(table["y"])
    mask = np.isclose(table["x"], xline) & (table["y"] >= y_surface - 0.5 * dy_tol)
    rows = np.sort(table[mask], order="y")
    s = rows["y"] - y_surface
    mass_flux_y = -(rows["nrho"] * H2O_MASS_PER_MOLECULE * rows["v"])
    return s, rows["trot"], rows["temp"], rows["press"], mass_flux_y


def select_axis_diag45(table: np.ndarray, x0: float, y0: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    tol = 0.5 * max(positive_spacing(table["x"]), positive_spacing(table["y"]))
    perp = np.abs((table["y"] - y0) - (table["x"] - x0)) / math.sqrt(2.0)
    s = ((table["x"] - x0) + (table["y"] - y0)) / math.sqrt(2.0)
    mask = (perp <= tol) & (s >= -0.5 * tol)
    rows = table[mask]
    s = s[mask]
    order = np.argsort(s)
    mass_flux_y = -(rows["nrho"] * H2O_MASS_PER_MOLECULE * rows["v"])
    return s[order], rows["trot"][order], rows["temp"][order], rows["press"][order], mass_flux_y[order]


def write_profile_table(
    path: Path,
    distance: np.ndarray,
    rotational_temperature: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
    mass_flux_y: np.ndarray,
) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("dist_m press_Pa temp_K trot_K mass_flux_y_kg_m2_s\n")
        for s, trot, temp, press, flux in zip(distance, rotational_temperature, temperature, pressure, mass_flux_y):
            handle.write(f"{s} {press} {temp} {trot} {flux}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot steady-state temperature/pressure profiles from the averaged condensation2d VTK frames.")
    parser.add_argument("case_dir", type=Path, help="Case directory, e.g. run/condensation2d/cases/.../v330K")
    parser.add_argument("--output-dir", default="profiles_steady", help="Output directory relative to the case dir")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    case_dir = args.case_dir.resolve()
    output_dir = case_dir / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    legacy_combined_path = output_dir / "steady_profiles.dat"
    if legacy_combined_path.exists():
        legacy_combined_path.unlink()

    metadata = load_json(case_dir / "metadata.json")
    rows, averaged_timesteps = load_averaged_grid_rows(case_dir)
    table = build_cell_table(rows, metadata)

    x_center, y_center = droplet_center(metadata)
    radius = metadata["radius"]
    x0_x, y0_x = surface_point(x_center, y_center, radius, 1.0, 0.0)
    x0_y, y0_y = surface_point(x_center, y_center, radius, 0.0, 1.0)
    x0_diag, y0_diag = surface_point(x_center, y_center, radius, 1.0, 1.0)

    profiles = {
        "x_axis": select_axis_x(table, x0_x, y0_x),
        "y_axis": select_axis_y(table, x0_y, y0_y),
        "diag_45deg": select_axis_diag45(table, x0_diag, y0_diag),
    }

    for name, values in profiles.items():
        distance, rotational_temperature, temperature, pressure, mass_flux_y = values
        write_profile_table(output_dir / f"{name}.dat", distance, rotational_temperature, temperature, pressure, mass_flux_y)

    print(f"Averaged steady frames: {len(averaged_timesteps)}")
    print(f"First averaged timestep: {averaged_timesteps[0]}")
    print(f"Last averaged timestep: {averaged_timesteps[-1]}")
    print(f"Droplet center used: x={x_center:.12g}, y={y_center:.12g}")
    for name in profiles:
        print(f"Wrote profile table to: {output_dir / f'{name}.dat'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
