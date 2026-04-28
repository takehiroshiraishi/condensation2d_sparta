#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def averaging_vtr_paths(case_dir: Path) -> list[Path]:
    vtr_dir = case_dir / "vtk_series" / "grid_steady"
    candidates = sorted(vtr_dir.glob("grid_steady_*.vtr"))
    if not candidates:
        raise FileNotFoundError(f"No grid_steady_*.vtr files found in {vtr_dir}")

    nonzero_candidates = [path for path in candidates if path.stem != "grid_steady_0000000000"]
    if not nonzero_candidates:
        return candidates[-1:]
    if len(nonzero_candidates) == 1:
        return nonzero_candidates
    return nonzero_candidates[1:]


def parse_data_array(piece: ET.Element, name: str) -> np.ndarray:
    for array in piece.iterfind(".//DataArray"):
        if array.attrib.get("Name") == name:
            text = array.text or ""
            return np.fromstring(text, sep=" ")
    raise KeyError(f"Could not find DataArray '{name}'")


def load_rectilinear_grid(vtr_path: Path) -> dict:
    tree = ET.parse(vtr_path)
    root = tree.getroot()
    piece = root.find(".//Piece")
    if piece is None:
        raise ValueError(f"Could not find Piece element in {vtr_path}")

    x_bounds = parse_data_array(piece, "x_coordinates")
    y_bounds = parse_data_array(piece, "y_coordinates")
    z_bounds = parse_data_array(piece, "z_coordinates")

    x_centers = 0.5 * (x_bounds[:-1] + x_bounds[1:])
    y_centers = 0.5 * (y_bounds[:-1] + y_bounds[1:])
    z_centers = 0.5 * (z_bounds[:-1] + z_bounds[1:])

    nx = len(x_centers)
    ny = len(y_centers)
    nz = len(z_centers)

    temp = parse_data_array(piece, "temp").reshape((nz, ny, nx))
    press = parse_data_array(piece, "press").reshape((nz, ny, nx))

    return {
        "x_centers": x_centers,
        "y_centers": y_centers,
        "z_centers": z_centers,
        "temp": temp,
        "press": press,
    }


def load_averaged_rectilinear_grid(vtr_paths: list[Path]) -> dict:
    grids = [load_rectilinear_grid(path) for path in vtr_paths]
    reference = grids[0]

    avg_temp = np.mean([grid["temp"] for grid in grids], axis=0)
    avg_press = np.mean([grid["press"] for grid in grids], axis=0)

    return {
        "x_centers": reference["x_centers"],
        "y_centers": reference["y_centers"],
        "z_centers": reference["z_centers"],
        "temp": avg_temp,
        "press": avg_press,
    }


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


def build_cell_table(grid: dict) -> np.ndarray:
    x_centers = grid["x_centers"]
    y_centers = grid["y_centers"]
    xx, yy = np.meshgrid(x_centers, y_centers, indexing="xy")

    # 2D cases have nz = 1, so take the single z slice.
    temp = grid["temp"][0]
    press = grid["press"][0]

    dtype = [
        ("x", float),
        ("y", float),
        ("temp", float),
        ("press", float),
    ]
    table = np.empty(xx.size, dtype=dtype)
    table["x"] = xx.ravel()
    table["y"] = yy.ravel()
    table["temp"] = temp.ravel()
    table["press"] = press.ravel()
    return table


def select_axis_x(table: np.ndarray, x0: float, y0: float, dx: float, dy: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    yline = table["y"][np.argmin(np.abs(table["y"] - y0))]
    x_surface = x0
    mask = (np.abs(table["y"] - yline) <= 0.25 * dy) & (table["x"] >= x_surface - 0.5 * dx)
    rows = np.sort(table[mask], order="x")
    s = rows["x"] - x_surface
    return s, rows["temp"], rows["press"]


def select_axis_y(table: np.ndarray, x0: float, y0: float, dx: float, dy: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xline = table["x"][np.argmin(np.abs(table["x"] - x0))]
    y_surface = y0
    mask = (np.abs(table["x"] - xline) <= 0.25 * dx) & (table["y"] >= y_surface - 0.5 * dy)
    rows = np.sort(table[mask], order="y")
    s = rows["y"] - y_surface
    return s, rows["temp"], rows["press"]


def select_axis_diag45(table: np.ndarray, x0: float, y0: float, dx: float, dy: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    tol = 0.5 * max(dx, dy)
    perp = np.abs((table["y"] - y0) - (table["x"] - x0)) / math.sqrt(2.0)
    s = ((table["x"] - x0) + (table["y"] - y0)) / math.sqrt(2.0)
    mask = (perp <= tol) & (s >= -0.5 * max(dx, dy))
    rows = table[mask]
    s = s[mask]
    order = np.argsort(s)
    return s[order], rows["temp"][order], rows["press"][order]


def write_profile_table(path: Path, distance: np.ndarray, temperature: np.ndarray, pressure: np.ndarray) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("dist_m press_Pa temp_K\n")
        for s, temp, press in zip(distance, temperature, pressure):
            handle.write(f"{s} {press} {temp}\n")


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
    vtr_paths = averaging_vtr_paths(case_dir)
    grid = load_averaged_rectilinear_grid(vtr_paths)
    table = build_cell_table(grid)

    x_center, y_center = droplet_center(metadata)
    radius = metadata["radius"]
    x0_x, y0_x = surface_point(x_center, y_center, radius, 1.0, 0.0)
    x0_y, y0_y = surface_point(x_center, y_center, radius, 0.0, 1.0)
    x0_diag, y0_diag = surface_point(x_center, y_center, radius, 1.0, 1.0)
    dx = float(grid["x_centers"][1] - grid["x_centers"][0]) if len(grid["x_centers"]) > 1 else 1.0
    dy = float(grid["y_centers"][1] - grid["y_centers"][0]) if len(grid["y_centers"]) > 1 else 1.0

    profiles = {
        "x_axis": select_axis_x(table, x0_x, y0_x, dx, dy),
        "y_axis": select_axis_y(table, x0_y, y0_y, dx, dy),
        "diag_45deg": select_axis_diag45(table, x0_diag, y0_diag, dx, dy),
    }

    for name, values in profiles.items():
        distance, temperature, pressure = values
        write_profile_table(output_dir / f"{name}.dat", distance, temperature, pressure)

    print(f"Averaged steady frames: {len(vtr_paths)}")
    print(f"First averaged frame: {vtr_paths[0]}")
    print(f"Last averaged frame: {vtr_paths[-1]}")
    print(f"Droplet center used: x={x_center:.12g}, y={y_center:.12g}")
    for name in profiles:
        print(f"Wrote profile table to: {output_dir / f'{name}.dat'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
