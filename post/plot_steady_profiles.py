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


def latest_vtr_path(case_dir: Path) -> Path:
    vtr_dir = case_dir / "vtk_series" / "grid_steady"
    candidates = sorted(vtr_dir.glob("grid_steady_*.vtr"))
    if not candidates:
        raise FileNotFoundError(f"No grid_steady_*.vtr files found in {vtr_dir}")
    return candidates[-1]


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


def droplet_center(metadata: dict) -> tuple[float, float]:
    x_center = metadata["droplets"][0]["center_x"]
    radius = metadata["radius"]
    theta_deg = metadata["contact_angle_deg"]
    surface_gap = metadata["surface_gap"]
    y_center = surface_gap - radius * math.cos(math.radians(theta_deg))
    return x_center, y_center


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
    mask = (np.abs(table["y"] - yline) <= 0.25 * dy) & (table["x"] >= x0 - 0.5 * dx)
    rows = np.sort(table[mask], order="x")
    s = rows["x"] - x0
    return s, rows["temp"], rows["press"]


def select_axis_y(table: np.ndarray, x0: float, y0: float, dx: float, dy: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xline = table["x"][np.argmin(np.abs(table["x"] - x0))]
    mask = (np.abs(table["x"] - xline) <= 0.25 * dx) & (table["y"] >= y0 - 0.5 * dy)
    rows = np.sort(table[mask], order="y")
    s = rows["y"] - y0
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


def write_combined_table(path: Path, profiles: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]) -> None:
    npoints = min(len(values[0]) for values in profiles.values())

    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "dist_x_m press_x_Pa temp_x_K "
            "dist_y_m press_y_Pa temp_y_K "
            "dist_45_m press_45_Pa temp_45_K\n"
        )

        for i in range(npoints):
            row: list[float | str] = []
            for key in ("x_axis", "y_axis", "diag_45deg"):
                s, temp, press = profiles[key]
                row.extend([s[i], press[i], temp[i]])
            handle.write(" ".join(str(value) for value in row) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot steady-state temperature/pressure profiles from the latest condensation2d VTK frame.")
    parser.add_argument("case_dir", type=Path, help="Case directory, e.g. run/condensation2d/cases/.../v330K")
    parser.add_argument("--output-dir", default="profiles_steady", help="Output directory relative to the case dir")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    case_dir = args.case_dir.resolve()
    output_dir = case_dir / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata = load_json(case_dir / "metadata.json")
    vtr_path = latest_vtr_path(case_dir)
    grid = load_rectilinear_grid(vtr_path)
    table = build_cell_table(grid)

    x0, y0 = droplet_center(metadata)
    dx = float(grid["x_centers"][1] - grid["x_centers"][0]) if len(grid["x_centers"]) > 1 else 1.0
    dy = float(grid["y_centers"][1] - grid["y_centers"][0]) if len(grid["y_centers"]) > 1 else 1.0

    profiles = {
        "x_axis": select_axis_x(table, x0, y0, dx, dy),
        "y_axis": select_axis_y(table, x0, y0, dx, dy),
        "diag_45deg": select_axis_diag45(table, x0, y0, dx, dy),
    }

    write_combined_table(output_dir / "steady_profiles.dat", profiles)

    print(f"Read steady frame: {vtr_path}")
    print(f"Droplet center used: x={x0:.12g}, y={y0:.12g}")
    print(f"Wrote combined profile table to: {output_dir / 'steady_profiles.dat'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
