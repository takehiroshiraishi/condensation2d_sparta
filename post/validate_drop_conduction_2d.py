#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from validate_drop_conduction_1d import (  # noqa: E402
    DEFAULT_CONDUCTIVITY,
    DEFAULT_LATENT,
    load_json,
    parse_surf_dump,
    selected_flux_rows,
    segment_length,
    surface_edges,
    width_at_y,
)


def build_grid(geom_rows: list[dict[str, float]], grid: float) -> tuple[int, int, float, float, list[int]]:
    yedges = surface_edges(geom_rows)
    ylo = yedges[0]
    yhi = yedges[-1]
    xhi = max(width_at_y(geom_rows, y) for y in yedges)
    nx = math.ceil(xhi / grid)
    ny = math.ceil((yhi - ylo) / grid)
    if nx <= 0 or ny <= 0:
        raise ValueError("Invalid 2D liquid grid")

    mask = [0] * (nx * ny)
    for iy in range(ny):
        y = ylo + (iy + 0.5) * grid
        local_width = width_at_y(geom_rows, y)
        for ix in range(nx):
            x = (ix + 0.5) * grid
            if x <= local_width + 1.0e-14:
                mask[iy * nx + ix] = 1

    if not any(mask):
        raise ValueError("No liquid cells found in 2D grid")
    return nx, ny, ylo, yhi, mask


def find_cell(
    x: float,
    y: float,
    nx: int,
    ny: int,
    ylo: float,
    grid: float,
    mask: list[int],
) -> int:
    ix = min(max(math.floor(x / grid), 0), nx - 1)
    iy = min(max(math.floor((y - ylo) / grid), 0), ny - 1)
    index = iy * nx + ix
    if mask[index]:
        return index

    best = -1
    best_dist2 = float("inf")
    for radius in range(1, nx + ny + 1):
        for jy in range(iy - radius, iy + radius + 1):
            if jy < 0 or jy >= ny:
                continue
            for jx in range(ix - radius, ix + radius + 1):
                if jx < 0 or jx >= nx:
                    continue
                if abs(jx - ix) != radius and abs(jy - iy) != radius:
                    continue
                candidate = jy * nx + jx
                if not mask[candidate]:
                    continue
                cx = (jx + 0.5) * grid
                cy = ylo + (jy + 0.5) * grid
                dist2 = (cx - x) ** 2 + (cy - y) ** 2
                if dist2 < best_dist2:
                    best = candidate
                    best_dist2 = dist2
        if best >= 0:
            return best
    return -1


def deposit_heat(
    geom_rows: list[dict[str, float]],
    flux_by_id: dict[int, float],
    nx: int,
    ny: int,
    ylo: float,
    grid: float,
    mask: list[int],
    latent: float,
) -> tuple[list[float], float]:
    source = [0.0] * (nx * ny)
    total_heat = 0.0
    for geom in geom_rows:
        surf_id = int(geom["id"])
        xmid = 0.5 * (geom["v1x"] + geom["v2x"])
        ymid = 0.5 * (geom["v1y"] + geom["v2y"])
        local_width = width_at_y(geom_rows, ymid)
        xin = xmid
        yin = ymid
        if xmid >= 0.5 * local_width:
            xin = xmid - 0.25 * grid
        else:
            yin = ymid - 0.25 * grid
        xin = max(0.0, xin)
        yin = max(ylo, yin)
        cell = find_cell(xin, yin, nx, ny, ylo, grid, mask)
        if cell < 0:
            continue
        heat = flux_by_id[surf_id] * segment_length(geom) * latent
        source[cell] += heat
        total_heat += heat
    return source, total_heat


def solve_steady_2d(
    nx: int,
    ny: int,
    mask: list[int],
    source: list[float],
    twall: float,
    conductivity: float,
    tolerance: float,
    maxiter: int,
    omega: float,
) -> tuple[list[float], int, float]:
    temperature = [twall] * (nx * ny)
    final_delta = 0.0
    for iteration in range(1, maxiter + 1):
        max_delta = 0.0
        for iy in range(ny):
            for ix in range(nx):
                index = iy * nx + ix
                if not mask[index]:
                    continue

                sum_g = 0.0
                sum_gt = source[index]

                if ix > 0 and mask[iy * nx + ix - 1]:
                    sum_g += conductivity
                    sum_gt += conductivity * temperature[iy * nx + ix - 1]
                if ix < nx - 1 and mask[iy * nx + ix + 1]:
                    sum_g += conductivity
                    sum_gt += conductivity * temperature[iy * nx + ix + 1]

                if iy == 0:
                    sum_g += 2.0 * conductivity
                    sum_gt += 2.0 * conductivity * twall
                elif mask[(iy - 1) * nx + ix]:
                    sum_g += conductivity
                    sum_gt += conductivity * temperature[(iy - 1) * nx + ix]

                if iy < ny - 1 and mask[(iy + 1) * nx + ix]:
                    sum_g += conductivity
                    sum_gt += conductivity * temperature[(iy + 1) * nx + ix]

                if sum_g <= 0.0:
                    continue
                solved = sum_gt / sum_g
                updated = temperature[index] + omega * (solved - temperature[index])
                delta = abs(updated - temperature[index])
                max_delta = max(max_delta, delta)
                temperature[index] = updated

        final_delta = max_delta
        if max_delta < tolerance:
            return temperature, iteration, final_delta
    return temperature, maxiter, final_delta


def write_cell_output(
    path: Path,
    nx: int,
    ny: int,
    ylo: float,
    grid: float,
    mask: list[int],
    source: list[float],
    temperature: list[float],
) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# x_m y_m mask heat_W_m temperature_K delta_T_K\n")
        for iy in range(ny):
            y = ylo + (iy + 0.5) * grid
            for ix in range(nx):
                index = iy * nx + ix
                x = (ix + 0.5) * grid
                handle.write(
                    f"{x:.16e} {y:.16e} {mask[index]} {source[index]:.16e} "
                    f"{temperature[index]:.16e} {temperature[index]-temperature[0]:.16e}\n"
                )
            handle.write("\n")


def write_surface_output(
    path: Path,
    geom_rows: list[dict[str, float]],
    flux_by_id: dict[int, float],
    nx: int,
    ny: int,
    ylo: float,
    grid: float,
    mask: list[int],
    source: list[float],
    temperature: list[float],
    twall: float,
) -> None:
    del source
    rows = []
    for geom in geom_rows:
        xmid = 0.5 * (geom["v1x"] + geom["v2x"])
        ymid = 0.5 * (geom["v1y"] + geom["v2y"])
        local_width = width_at_y(geom_rows, ymid)
        xin = xmid
        yin = ymid
        if xmid >= 0.5 * local_width:
            xin = xmid - 0.25 * grid
        else:
            yin = ymid - 0.25 * grid
        cell = find_cell(max(0.0, xin), max(ylo, yin), nx, ny, ylo, grid, mask)
        temp = temperature[cell] if cell >= 0 else twall
        rows.append((ymid, xmid, int(geom["id"]), flux_by_id[int(geom["id"])], temp))
    rows.sort()

    with path.open("w", encoding="utf-8") as handle:
        handle.write("# y_mid_m x_mid_m surf_id mflux_kg_m2_s surface_temperature_K delta_T_K\n")
        for ymid, xmid, surf_id, flux, temp in rows:
            handle.write(f"{ymid:.16e} {xmid:.16e} {surf_id} {flux:.16e} {temp:.16e} {temp-twall:.16e}\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate the steady 2D liquid-conduction model using fixed SPARTA surface flux dumps."
    )
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("--droplet-index", type=int, default=1)
    parser.add_argument("--average-frames", action="store_true")
    parser.add_argument("--grid", type=float, default=0.25e-6, help="Cartesian liquid grid spacing [m]")
    parser.add_argument("--twall", type=float, default=None)
    parser.add_argument("--latent", type=float, default=DEFAULT_LATENT)
    parser.add_argument("--conductivity", type=float, default=DEFAULT_CONDUCTIVITY)
    parser.add_argument("--tolerance", type=float, default=1.0e-8)
    parser.add_argument("--maxiter", type=int, default=50000)
    parser.add_argument("--omega", type=float, default=1.0, help="SOR relaxation factor; 1.0 is Gauss-Seidel")
    parser.add_argument("--output-prefix", type=Path, default=None)
    args = parser.parse_args()

    case_dir = args.case_dir.resolve()
    metadata = load_json(case_dir / "metadata.json")
    twall = args.twall if args.twall is not None else float(metadata["temperature_k"])
    flux_path = case_dir / f"surf_droplet{args.droplet_index}.dump"
    geom_path = case_dir / f"surf_geom_droplet{args.droplet_index}.dump"

    flux_frames = parse_surf_dump(flux_path)
    geom_frames = parse_surf_dump(geom_path)
    frame_label, flux_rows = selected_flux_rows(flux_frames, args.average_frames)
    geom_rows = geom_frames[-1]["rows"]  # type: ignore[assignment]
    flux_columns = flux_frames[-1]["columns"]
    flux_column = f"f_avg_droplet_{args.droplet_index}[1]"
    integrated_column = f"f_avg_droplet_{args.droplet_index}[2]"
    if flux_column not in flux_columns:
        raise KeyError(f"Missing {flux_column} in {flux_path}")

    geom_by_id = {int(row["id"]): row for row in geom_rows}  # type: ignore[union-attr]
    flux_by_id = {int(row["id"]): row[flux_column] for row in flux_rows}
    integrated_by_id = {}
    if integrated_column in flux_columns:
        integrated_by_id = {int(row["id"]): row[integrated_column] for row in flux_rows}
    if geom_by_id.keys() != flux_by_id.keys():
        raise ValueError("Surface IDs do not match between flux and geometry dumps")

    nx, ny, ylo, _, mask = build_grid(list(geom_by_id.values()), args.grid)
    source, total_heat = deposit_heat(
        list(geom_by_id.values()),
        flux_by_id,
        nx,
        ny,
        ylo,
        args.grid,
        mask,
        args.latent,
    )
    temperature, iterations, final_delta = solve_steady_2d(
        nx,
        ny,
        mask,
        source,
        twall,
        args.conductivity,
        args.tolerance,
        args.maxiter,
        args.omega,
    )

    liquid_temps = [temperature[i] for i, inside in enumerate(mask) if inside]
    liquid_area = sum(mask) * args.grid * args.grid
    total_mass = total_heat / args.latent
    total_integrated_mass = sum(integrated_by_id.values()) if integrated_by_id else 0.0
    mass_error = 0.0
    if integrated_by_id and abs(total_integrated_mass) > 0.0:
        mass_error = abs(total_mass - total_integrated_mass) / abs(total_integrated_mass)

    # Bottom wall heat removal from fixed-T half-cell conductance.
    wall_heat = 0.0
    for ix in range(nx):
        index = ix
        if mask[index]:
            wall_heat += 2.0 * args.conductivity * (temperature[index] - twall)
    heat_balance_error = abs(wall_heat - total_heat) / max(abs(total_heat), 1.0e-300)

    output_prefix = args.output_prefix
    if output_prefix is None:
        output_prefix = case_dir / "profiles_steady" / "validate_drop_conduction_2d"
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    cell_output = output_prefix.with_suffix(".cells.dat")
    surface_output = output_prefix.with_suffix(".surface.dat")
    write_cell_output(cell_output, nx, ny, ylo, args.grid, mask, source, temperature)
    write_surface_output(
        surface_output,
        list(geom_by_id.values()),
        flux_by_id,
        nx,
        ny,
        ylo,
        args.grid,
        mask,
        source,
        temperature,
        twall,
    )

    print(f"case: {case_dir}")
    print(f"flux frame: {frame_label}")
    print(f"grid: {nx} x {ny}, dx = {args.grid:.8e} m")
    print(f"liquid cells: {sum(mask)}")
    print(f"liquid area from Cartesian cells [m2/m]: {liquid_area:.8e}")
    print(f"total mass rate from mflux*arc_length [kg/m/s]: {total_mass:.8e}")
    if integrated_by_id:
        print(f"total mass rate from integrated dump column [kg/m/s]: {total_integrated_mass:.8e}")
        print(f"relative mass-rate mismatch: {mass_error:.8e}")
    print(f"total latent heat input [W/m]: {total_heat:.8e}")
    print(f"wall heat removal [W/m]: {wall_heat:.8e}")
    print(f"relative steady heat-balance mismatch: {heat_balance_error:.8e}")
    print(f"iterations: {iterations}")
    print(f"final max update [K]: {final_delta:.8e}")
    print(f"2D liquid Tmin/Tmax [K]: {min(liquid_temps):.6f} {max(liquid_temps):.6f}")
    print(f"2D max rise [K]: {max(liquid_temps)-twall:.6f}")
    print(f"wrote: {cell_output}")
    print(f"wrote: {surface_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
