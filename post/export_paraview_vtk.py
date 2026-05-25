#!/usr/bin/env python3

from __future__ import annotations

import argparse
from itertools import product
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_steady_profiles import build_cell_table, load_json, map_frame_fields, parse_grid_dump_frames


def infer_bounds(centers: list[float], spacing: float) -> list[float]:
    if len(centers) == 1:
        return [centers[0] - 0.5 * spacing, centers[0] + 0.5 * spacing]
    bounds = [centers[0] - 0.5 * (centers[1] - centers[0])]
    for left, right in zip(centers, centers[1:]):
        bounds.append(0.5 * (left + right))
    bounds.append(centers[-1] + 0.5 * (centers[-1] - centers[-2]))
    return bounds


def build_frame(raw_frame: dict, metadata: dict) -> dict:
    rows = map_frame_fields(raw_frame["rows"])
    x_centers = sorted({row["xc"] for row in rows})
    y_centers = sorted({row["yc"] for row in rows})
    z_centers = sorted({row["zc"] for row in rows})
    bounds = metadata["simulation_bounds"]
    nx, ny, nz = metadata["grid_cells"]
    dx = (bounds["xhi"] - bounds["xlo"]) / nx
    dy = (bounds["yhi"] - bounds["ylo"]) / ny
    dz = (bounds["zhi"] - bounds["zlo"]) / max(nz, 1)

    is_rectilinear = len(rows) == len(x_centers) * len(y_centers) * len(z_centers)
    ordered_rows = []
    if is_rectilinear:
        row_map = {(row["xc"], row["yc"], row["zc"]): row for row in rows}
        ordered_rows = [row_map[(x, y, z)] for z, y, x in product(z_centers, y_centers, x_centers)]
    else:
        ordered_rows = sorted(rows, key=lambda row: int(row["id"]))
    for row in ordered_rows:
        row["w"] = 0.0

    frame = {
        "timestep": raw_frame["timestep"],
        "ordered_rows": ordered_rows,
        "is_rectilinear": is_rectilinear,
    }
    if is_rectilinear:
        frame["x_bounds"] = infer_bounds(x_centers, dx)
        frame["y_bounds"] = infer_bounds(y_centers, dy)
        frame["z_bounds"] = infer_bounds(z_centers, dz)
    else:
        table = build_cell_table(ordered_rows, metadata)
        frame["cell_table"] = table
        frame["z_bounds"] = [bounds["zlo"], bounds["zhi"]]
    return frame


def write_legacy_vtk(path: Path, frame: dict) -> None:
    if not frame["is_rectilinear"]:
        raise ValueError("Legacy rectilinear VTK output is not available for refined grids; use --mode all.")
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# vtk DataFile Version 3.0\n")
        handle.write(f"SPARTA averaged grid dump timestep {frame['timestep']}\n")
        handle.write("ASCII\n")
        handle.write("DATASET RECTILINEAR_GRID\n")
        handle.write(f"DIMENSIONS {len(frame['x_bounds'])} {len(frame['y_bounds'])} {len(frame['z_bounds'])}\n")
        for axis_name, values in (("X", frame["x_bounds"]), ("Y", frame["y_bounds"]), ("Z", frame["z_bounds"])):
            handle.write(f"{axis_name}_COORDINATES {len(values)} float\n")
            handle.write(" ".join(f"{value:.12g}" for value in values) + "\n")
        handle.write(f"CELL_DATA {len(frame['ordered_rows'])}\n")
        for field_name in ("nrho", "trot", "temp", "press"):
            handle.write(f"SCALARS {field_name} float 1\nLOOKUP_TABLE default\n")
            for row in frame["ordered_rows"]:
                handle.write(f"{row.get(field_name, 0.0):.12g}\n")
        handle.write("VECTORS velocity float\n")
        for row in frame["ordered_rows"]:
            handle.write(f"{row['u']:.12g} {row['v']:.12g} 0.0\n")


def write_vtr(path: Path, frame: dict) -> None:
    if not frame["is_rectilinear"]:
        raise ValueError("Rectilinear VTR output is not available for refined grids.")
    whole_extent = f"0 {len(frame['x_bounds'])-1} 0 {len(frame['y_bounds'])-1} 0 {len(frame['z_bounds'])-1}"
    with path.open("w", encoding="utf-8") as handle:
        handle.write('<?xml version="1.0"?>\n')
        handle.write('<VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">\n')
        handle.write(f'  <RectilinearGrid WholeExtent="{whole_extent}">\n')
        handle.write(f'    <Piece Extent="{whole_extent}">\n')
        handle.write("      <CellData>\n")
        for field_name in ("nrho", "trot", "temp", "press"):
            handle.write(f'      <DataArray type="Float64" Name="{field_name}" format="ascii">\n')
            handle.write("        " + " ".join(f"{row.get(field_name, 0.0):.12g}" for row in frame["ordered_rows"]) + "\n")
            handle.write("      </DataArray>\n")
        handle.write('      <DataArray type="Float64" Name="velocity" NumberOfComponents="3" format="ascii">\n')
        handle.write("        " + " ".join(f"{row['u']:.12g} {row['v']:.12g} 0.0" for row in frame["ordered_rows"]) + "\n")
        handle.write("      </DataArray>\n")
        handle.write("      </CellData>\n")
        handle.write("      <Coordinates>\n")
        for axis_name, values in (("x", frame["x_bounds"]), ("y", frame["y_bounds"]), ("z", frame["z_bounds"])):
            handle.write(f'      <DataArray type="Float64" Name="{axis_name}_coordinates" format="ascii">\n')
            handle.write("        " + " ".join(f"{value:.12g}" for value in values) + "\n")
            handle.write("      </DataArray>\n")
        handle.write("      </Coordinates>\n")
        handle.write("    </Piece>\n")
        handle.write("  </RectilinearGrid>\n")
        handle.write("</VTKFile>\n")


def write_vtu(path: Path, frame: dict) -> None:
    table = frame["cell_table"]
    points: list[tuple[float, float, float]] = []
    connectivity: list[int] = []
    offsets: list[int] = []
    cell_types: list[int] = []
    point_index = 0
    zlo, zhi = frame["z_bounds"]
    for row in table:
        xlo = float(row["xlo"])
        xhi = float(row["xhi"])
        ylo = float(row["ylo"])
        yhi = float(row["yhi"])
        points.extend(
            [
                (xlo, ylo, zlo),
                (xhi, ylo, zlo),
                (xhi, yhi, zlo),
                (xlo, yhi, zlo),
                (xlo, ylo, zhi),
                (xhi, ylo, zhi),
                (xhi, yhi, zhi),
                (xlo, yhi, zhi),
            ]
        )
        connectivity.extend(range(point_index, point_index + 8))
        point_index += 8
        offsets.append(point_index)
        cell_types.append(12)

    with path.open("w", encoding="utf-8") as handle:
        handle.write('<?xml version="1.0"?>\n')
        handle.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        handle.write("  <UnstructuredGrid>\n")
        handle.write(f'    <Piece NumberOfPoints="{len(points)}" NumberOfCells="{len(table)}">\n')
        handle.write("      <CellData>\n")
        for field_name in ("nrho", "trot", "temp", "press"):
            handle.write(f'      <DataArray type="Float64" Name="{field_name}" format="ascii">\n')
            handle.write("        " + " ".join(f"{float(row[field_name]):.12g}" for row in table) + "\n")
            handle.write("      </DataArray>\n")
        handle.write('      <DataArray type="Float64" Name="velocity" NumberOfComponents="3" format="ascii">\n')
        handle.write("        " + " ".join(f"{float(row['u']):.12g} {float(row['v']):.12g} 0.0" for row in table) + "\n")
        handle.write("      </DataArray>\n")
        handle.write("      </CellData>\n")
        handle.write("      <Points>\n")
        handle.write('      <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        handle.write("        " + " ".join(f"{x:.12g} {y:.12g} {z:.12g}" for x, y, z in points) + "\n")
        handle.write("      </DataArray>\n")
        handle.write("      </Points>\n")
        handle.write("      <Cells>\n")
        handle.write('      <DataArray type="Int64" Name="connectivity" format="ascii">\n')
        handle.write("        " + " ".join(str(value) for value in connectivity) + "\n")
        handle.write("      </DataArray>\n")
        handle.write('      <DataArray type="Int64" Name="offsets" format="ascii">\n')
        handle.write("        " + " ".join(str(value) for value in offsets) + "\n")
        handle.write("      </DataArray>\n")
        handle.write('      <DataArray type="UInt8" Name="types" format="ascii">\n')
        handle.write("        " + " ".join(str(value) for value in cell_types) + "\n")
        handle.write("      </DataArray>\n")
        handle.write("      </Cells>\n")
        handle.write("    </Piece>\n")
        handle.write("  </UnstructuredGrid>\n")
        handle.write("</VTKFile>\n")


def write_pvd(path: Path, entries: list[tuple[int, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write('<?xml version="1.0"?>\n')
        handle.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n  <Collection>\n')
        for timestep, filename in entries:
            handle.write(f'    <DataSet timestep="{timestep}" group="" part="0" file="{filename}"/>\n')
        handle.write("  </Collection>\n</VTKFile>\n")


def export_last(frames: list[dict], metadata: dict, output_path: Path) -> None:
    frame = build_frame(frames[-1], metadata)
    if frame["is_rectilinear"]:
        write_legacy_vtk(output_path, frame)
        return
    fallback_path = output_path.with_suffix(".vtu")
    write_vtu(fallback_path, frame)


def export_all(frames: list[dict], metadata: dict, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    entries = []
    for raw_frame in frames:
        frame = build_frame(raw_frame, metadata)
        if frame["is_rectilinear"]:
            filename = f"{stem}_{frame['timestep']:010d}.vtr"
            write_vtr(output_dir / filename, frame)
        else:
            filename = f"{stem}_{frame['timestep']:010d}.vtu"
            write_vtu(output_dir / filename, frame)
        entries.append((frame["timestep"], filename))
    write_pvd(output_dir / f"{stem}.pvd", entries)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert SPARTA condensation2d grid dumps to ParaView-readable files.")
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("--grid-dump", default="grid_steady.dump")
    parser.add_argument("--slice-dump", default="line_x0.dump")
    parser.add_argument("--grid-output", default="grid_steady_legacy.vtk")
    parser.add_argument("--slice-output", default="line_x0_legacy.vtk")
    parser.add_argument("--mode", choices=("last", "all"), default="last")
    parser.add_argument("--series-dir", default="vtk_series")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    case_dir = args.case_dir.resolve()
    metadata = load_json(case_dir / "metadata.json")
    grid_dump = case_dir / args.grid_dump
    line_dump = case_dir / args.slice_dump

    if args.mode == "last":
        if grid_dump.exists():
            export_last(parse_grid_dump_frames(grid_dump), metadata, case_dir / args.grid_output)
            if (case_dir / args.grid_output).exists():
                print(f"Wrote {case_dir / args.grid_output}")
            else:
                print(f"Wrote {(case_dir / args.grid_output).with_suffix('.vtu')}")
        if line_dump.exists():
            export_last(parse_grid_dump_frames(line_dump), metadata, case_dir / args.slice_output)
            if (case_dir / args.slice_output).exists():
                print(f"Wrote {case_dir / args.slice_output}")
            else:
                print(f"Wrote {(case_dir / args.slice_output).with_suffix('.vtu')}")
        return 0

    series_dir = case_dir / args.series_dir
    if grid_dump.exists():
        export_all(parse_grid_dump_frames(grid_dump), metadata, series_dir / "grid_steady", "grid_steady")
        print(f"Wrote {series_dir / 'grid_steady' / 'grid_steady.pvd'}")
    if line_dump.exists():
        export_all(parse_grid_dump_frames(line_dump), metadata, series_dir / "line_x0", "line_x0")
        print(f"Wrote {series_dir / 'line_x0' / 'line_x0.pvd'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
