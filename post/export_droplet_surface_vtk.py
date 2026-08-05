#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_last_surf_geom(path: Path) -> tuple[list[str], list[dict[str, float]], int]:
    last_columns: list[str] | None = None
    last_rows: list[dict[str, float]] | None = None
    last_timestep = 0

    with path.open("r", encoding="utf-8") as handle:
        lines = iter(handle)
        for line in lines:
            if line.strip() != "ITEM: TIMESTEP":
                continue
            last_timestep = int(next(lines).strip())
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
    for required in ("id", "v1x", "v1y", "v2x", "v2y"):
        if required not in last_columns:
            raise KeyError(f"Missing {required} in {path}")
    return last_columns, last_rows, last_timestep


def point_key(x: float, y: float) -> tuple[int, int]:
    # SPARTA text dumps use limited precision; quantize for stable endpoint reuse.
    return (round(x / 1.0e-15), round(y / 1.0e-15))


def build_polyline(rows: list[dict[str, float]]) -> tuple[list[tuple[float, float, float]], list[tuple[int, int, int]]]:
    points: list[tuple[float, float, float]] = []
    point_ids: dict[tuple[int, int], int] = {}
    lines: list[tuple[int, int, int]] = []

    for row in rows:
        segment_points = []
        for x_name, y_name in (("v1x", "v1y"), ("v2x", "v2y")):
            x = row[x_name]
            y = row[y_name]
            key = point_key(x, y)
            if key not in point_ids:
                point_ids[key] = len(points)
                points.append((x, y, 0.0))
            segment_points.append(point_ids[key])
        lines.append((int(row["id"]), segment_points[0], segment_points[1]))

    return points, lines


def write_legacy_vtk(path: Path, points: list[tuple[float, float, float]], lines: list[tuple[int, int, int]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# vtk DataFile Version 3.0\n")
        handle.write("SPARTA droplet surface geometry\n")
        handle.write("ASCII\n")
        handle.write("DATASET POLYDATA\n")
        handle.write(f"POINTS {len(points)} float\n")
        for x, y, z in points:
            handle.write(f"{x:.12g} {y:.12g} {z:.12g}\n")
        handle.write(f"LINES {len(lines)} {3 * len(lines)}\n")
        for _, p1, p2 in lines:
            handle.write(f"2 {p1} {p2}\n")
        handle.write(f"CELL_DATA {len(lines)}\n")
        handle.write("SCALARS surf_id int 1\n")
        handle.write("LOOKUP_TABLE default\n")
        for surf_id, _, _ in lines:
            handle.write(f"{surf_id}\n")


def write_vtp(path: Path, points: list[tuple[float, float, float]], lines: list[tuple[int, int, int]]) -> None:
    connectivity = []
    offsets = []
    for index, (_, p1, p2) in enumerate(lines, start=1):
        connectivity.extend([p1, p2])
        offsets.append(2 * index)

    with path.open("w", encoding="utf-8") as handle:
        handle.write('<?xml version="1.0"?>\n')
        handle.write('<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">\n')
        handle.write("  <PolyData>\n")
        handle.write(f'    <Piece NumberOfPoints="{len(points)}" NumberOfLines="{len(lines)}">\n')
        handle.write("      <CellData>\n")
        handle.write('        <DataArray type="Int64" Name="surf_id" format="ascii">\n')
        handle.write("          " + " ".join(str(surf_id) for surf_id, _, _ in lines) + "\n")
        handle.write("        </DataArray>\n")
        handle.write("      </CellData>\n")
        handle.write("      <Points>\n")
        handle.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        handle.write("          " + " ".join(f"{x:.12g} {y:.12g} {z:.12g}" for x, y, z in points) + "\n")
        handle.write("        </DataArray>\n")
        handle.write("      </Points>\n")
        handle.write("      <Lines>\n")
        handle.write('        <DataArray type="Int64" Name="connectivity" format="ascii">\n')
        handle.write("          " + " ".join(str(value) for value in connectivity) + "\n")
        handle.write("        </DataArray>\n")
        handle.write('        <DataArray type="Int64" Name="offsets" format="ascii">\n')
        handle.write("          " + " ".join(str(value) for value in offsets) + "\n")
        handle.write("        </DataArray>\n")
        handle.write("      </Lines>\n")
        handle.write("    </Piece>\n")
        handle.write("  </PolyData>\n")
        handle.write("</VTKFile>\n")


def polygon_from_surface(rows: list[dict[str, float]]) -> list[tuple[float, float, float]]:
    unique_points: dict[tuple[int, int], tuple[float, float, float]] = {}
    for row in rows:
        for x_name, y_name in (("v1x", "v1y"), ("v2x", "v2y")):
            x = row[x_name]
            y = row[y_name]
            unique_points[point_key(x, y)] = (x, y, 0.0)

    surface_points = sorted(unique_points.values(), key=lambda point: (point[1], point[0]))
    if not surface_points:
        raise ValueError("No droplet surface points found")

    contact = surface_points[0]
    polygon = [(0.0, 0.0, 0.0)]
    if abs(contact[0]) > 1.0e-15 or abs(contact[1]) > 1.0e-15:
        polygon.append(contact)
    polygon.extend(point for point in surface_points[1:] if abs(point[0]) > 1.0e-15 or abs(point[1]) > 1.0e-15)
    return polygon


def write_liquid_field_vtu(
    path: Path,
    polygon: list[tuple[float, float, float]],
    liquid_temperature: float,
) -> None:
    center_x = sum(point[0] for point in polygon) / len(polygon)
    center_y = sum(point[1] for point in polygon) / len(polygon)
    points = [(center_x, center_y, 0.0), *polygon]
    connectivity = []
    offsets = []
    cell_types = []
    cell_count = len(polygon)
    for index in range(cell_count):
        p1 = index + 1
        p2 = 1 if index == cell_count - 1 else index + 2
        connectivity.extend([0, p1, p2])
        offsets.append(3 * (index + 1))
        cell_types.append(5)  # VTK_TRIANGLE

    with path.open("w", encoding="utf-8") as handle:
        handle.write('<?xml version="1.0"?>\n')
        handle.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        handle.write("  <UnstructuredGrid>\n")
        handle.write(f'    <Piece NumberOfPoints="{len(points)}" NumberOfCells="{cell_count}">\n')
        handle.write("      <CellData>\n")
        for name, value in (
            ("press", 0.0),
            ("p_star", 0.0),
            ("phase", 1.0),
            ("temp", liquid_temperature),
        ):
            handle.write(f'        <DataArray type="Float64" Name="{name}" format="ascii">\n')
            handle.write("          " + " ".join(f"{value:.12g}" for _ in range(cell_count)) + "\n")
            handle.write("        </DataArray>\n")
        handle.write("      </CellData>\n")
        handle.write("      <Points>\n")
        handle.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        handle.write("          " + " ".join(f"{x:.12g} {y:.12g} {z:.12g}" for x, y, z in points) + "\n")
        handle.write("        </DataArray>\n")
        handle.write("      </Points>\n")
        handle.write("      <Cells>\n")
        handle.write('        <DataArray type="Int64" Name="connectivity" format="ascii">\n')
        handle.write("          " + " ".join(str(value) for value in connectivity) + "\n")
        handle.write("        </DataArray>\n")
        handle.write('        <DataArray type="Int64" Name="offsets" format="ascii">\n')
        handle.write("          " + " ".join(str(value) for value in offsets) + "\n")
        handle.write("        </DataArray>\n")
        handle.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
        handle.write("          " + " ".join(str(value) for value in cell_types) + "\n")
        handle.write("        </DataArray>\n")
        handle.write("      </Cells>\n")
        handle.write("    </Piece>\n")
        handle.write("  </UnstructuredGrid>\n")
        handle.write("</VTKFile>\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert SPARTA droplet surface geometry dumps to ParaView PolyData.")
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("--droplet-index", type=int, default=1)
    parser.add_argument("--output-prefix", default=None)
    parser.add_argument("--liquid-output", default=None)
    args = parser.parse_args()

    case_dir = args.case_dir.resolve()
    input_path = case_dir / f"surf_geom_droplet{args.droplet_index}.dump"
    if not input_path.exists():
        raise FileNotFoundError(f"Missing surface geometry dump: {input_path}")

    _, rows, timestep = parse_last_surf_geom(input_path)
    points, lines = build_polyline(rows)
    output_prefix = args.output_prefix or f"droplet_surface{args.droplet_index}"
    legacy_path = case_dir / f"{output_prefix}.vtk"
    vtp_path = case_dir / f"{output_prefix}.vtp"
    write_legacy_vtk(legacy_path, points, lines)
    write_vtp(vtp_path, points, lines)

    metadata_path = case_dir / "metadata.json"
    liquid_temperature = 0.0
    if metadata_path.exists():
        with metadata_path.open("r", encoding="utf-8") as handle:
            liquid_temperature = float(json.load(handle).get("temperature_k", 0.0))
    liquid_output = args.liquid_output or f"droplet_liquid_field{args.droplet_index}.vtu"
    liquid_path = case_dir / liquid_output
    write_liquid_field_vtu(liquid_path, polygon_from_surface(rows), liquid_temperature)

    print(f"Wrote {legacy_path}")
    print(f"Wrote {vtp_path}")
    print(f"Wrote {liquid_path}")
    print(f"Surface frame timestep: {timestep}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
