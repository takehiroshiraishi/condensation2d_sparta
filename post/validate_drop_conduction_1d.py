#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


DEFAULT_LATENT = 2.43e6
DEFAULT_CONDUCTIVITY = 0.6
DEFAULT_RHO = 997.0
DEFAULT_CP = 4180.0


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def parse_surf_dump(path: Path) -> list[dict[str, object]]:
    frames: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8") as handle:
        lines = iter(handle)
        for line in lines:
            if line.strip() != "ITEM: TIMESTEP":
                continue
            timestep = int(next(lines).strip())
            if next(lines).strip() != "ITEM: NUMBER OF SURFS":
                raise ValueError(f"Unexpected dump format in {path}")
            count = int(next(lines).strip())
            bounds_header = next(lines).strip()
            if not bounds_header.startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"Unexpected bounds header in {path}")
            bounds = [tuple(float(v) for v in next(lines).split()) for _ in range(3)]
            columns = next(lines).strip().split()[2:]
            rows = []
            for _ in range(count):
                raw = next(lines).split()
                rows.append({column: float(value) for column, value in zip(columns, raw)})
            frames.append({"timestep": timestep, "bounds": bounds, "columns": columns, "rows": rows})
    if not frames:
        raise ValueError(f"No frames found in {path}")
    return frames


def selected_flux_rows(frames: list[dict[str, object]], average_frames: bool) -> tuple[str, list[dict[str, float]]]:
    data_frames = [frame for frame in frames if int(frame["timestep"]) != 0]
    if not data_frames:
        data_frames = frames[-1:]

    if not average_frames:
        frame = data_frames[-1]
        return str(frame["timestep"]), frame["rows"]  # type: ignore[return-value]

    columns = data_frames[-1]["columns"]
    rows_by_id: dict[int, dict[str, float]] = {}
    counts: dict[int, int] = {}
    for frame in data_frames:
        for row in frame["rows"]:  # type: ignore[union-attr]
            surf_id = int(row["id"])
            if surf_id not in rows_by_id:
                rows_by_id[surf_id] = {column: 0.0 for column in columns}  # type: ignore[arg-type]
                counts[surf_id] = 0
            for column in columns:  # type: ignore[union-attr]
                rows_by_id[surf_id][column] += row[column]
            counts[surf_id] += 1

    rows = []
    for surf_id, row in rows_by_id.items():
        averaged = {column: value / counts[surf_id] for column, value in row.items()}
        averaged["id"] = float(surf_id)
        rows.append(averaged)
    rows.sort(key=lambda row: int(row["id"]))
    return f"average_nonzero_frames_{len(data_frames)}", rows


def surface_edges(geom_rows: list[dict[str, float]]) -> list[float]:
    values = []
    for row in geom_rows:
        values.append(row["v1y"])
        values.append(row["v2y"])
    values.sort()
    unique = []
    for value in values:
        if not unique or abs(value - unique[-1]) > 1.0e-14:
            unique.append(value)
    if len(unique) < 2:
        raise ValueError("Could not determine surface-based y bins")
    return unique


def x_at_y(row: dict[str, float], y: float) -> float:
    y1 = row["v1y"]
    y2 = row["v2y"]
    x1 = row["v1x"]
    x2 = row["v2x"]
    if abs(y2 - y1) < 1.0e-300:
        return 0.5 * (x1 + x2)
    fraction = (y - y1) / (y2 - y1)
    return x1 + fraction * (x2 - x1)


def width_at_y(geom_rows: list[dict[str, float]], y: float) -> float:
    candidates = []
    for row in geom_rows:
        ylo = min(row["v1y"], row["v2y"])
        yhi = max(row["v1y"], row["v2y"])
        if ylo - 1.0e-14 <= y <= yhi + 1.0e-14:
            candidates.append(max(0.0, x_at_y(row, y)))
    if not candidates:
        return 0.0
    return max(candidates)


def direct_slice_area(geom_rows: list[dict[str, float]], ylo: float, yhi: float) -> float:
    area = 0.0
    for row in geom_rows:
        sylo = min(row["v1y"], row["v2y"])
        syhi = max(row["v1y"], row["v2y"])
        alo = max(ylo, sylo)
        ahi = min(yhi, syhi)
        if ahi <= alo:
            continue
        xlo = x_at_y(row, alo)
        xhi = x_at_y(row, ahi)
        area += 0.5 * (xlo + xhi) * (ahi - alo)
    return area


def segment_length(row: dict[str, float]) -> float:
    return math.hypot(row["v2x"] - row["v1x"], row["v2y"] - row["v1y"])


def find_bin(yedges: list[float], y: float) -> int:
    if y <= yedges[0]:
        return 0
    if y >= yedges[-1]:
        return len(yedges) - 2
    for index in range(len(yedges) - 1):
        if yedges[index] <= y < yedges[index + 1]:
            return index
    return len(yedges) - 2


def solve_steady_temperature(
    yedges: list[float],
    widths: list[float],
    heat: list[float],
    twall: float,
    conductivity: float,
) -> list[float]:
    nbin = len(heat)
    lower = [0.0] * nbin
    diag = [0.0] * nbin
    upper = [0.0] * nbin
    rhs = [0.0] * nbin

    for i in range(nbin):
        center = 0.5 * (yedges[i] + yedges[i + 1])
        if i == 0:
            down_distance = center - yedges[i]
        else:
            down_center = 0.5 * (yedges[i - 1] + yedges[i])
            down_distance = center - down_center
        gdown = conductivity * widths[i] / down_distance if down_distance > 0.0 else 0.0

        gup = 0.0
        if i < nbin - 1:
            up_center = 0.5 * (yedges[i + 1] + yedges[i + 2])
            up_distance = up_center - center
            gup = conductivity * widths[i + 1] / up_distance if up_distance > 0.0 else 0.0

        lower[i] = 0.0 if i == 0 else -gdown
        upper[i] = 0.0 if i == nbin - 1 else -gup
        diag[i] = gdown + gup
        rhs[i] = heat[i]
        if i == 0:
            rhs[i] += gdown * twall

    cp = [0.0] * nbin
    dp = [0.0] * nbin
    cp[0] = upper[0] / diag[0]
    dp[0] = rhs[0] / diag[0]
    for i in range(1, nbin):
        denom = diag[i] - lower[i] * cp[i - 1]
        cp[i] = upper[i] / denom
        dp[i] = (rhs[i] - lower[i] * dp[i - 1]) / denom
    for i in range(nbin - 2, -1, -1):
        dp[i] -= cp[i] * dp[i + 1]
    return dp


def solve_transient_temperature(
    yedges: list[float],
    widths: list[float],
    volumes: list[float],
    heat: list[float],
    previous_temperature: float,
    twall: float,
    conductivity: float,
    liquid_rho: float,
    liquid_cp: float,
    dtcond: float,
) -> list[float]:
    nbin = len(heat)
    lower = [0.0] * nbin
    diag = [0.0] * nbin
    upper = [0.0] * nbin
    rhs = [0.0] * nbin

    for i in range(nbin):
        capdt = liquid_rho * liquid_cp * volumes[i] / dtcond
        center = 0.5 * (yedges[i] + yedges[i + 1])
        if i == 0:
            down_distance = center - yedges[i]
        else:
            down_center = 0.5 * (yedges[i - 1] + yedges[i])
            down_distance = center - down_center
        gdown = conductivity * widths[i] / down_distance if down_distance > 0.0 else 0.0

        gup = 0.0
        if i < nbin - 1:
            up_center = 0.5 * (yedges[i + 1] + yedges[i + 2])
            up_distance = up_center - center
            gup = conductivity * widths[i + 1] / up_distance if up_distance > 0.0 else 0.0

        lower[i] = 0.0 if i == 0 else -gdown
        upper[i] = 0.0 if i == nbin - 1 else -gup
        diag[i] = capdt + gdown + gup
        rhs[i] = capdt * previous_temperature + heat[i]
        if i == 0:
            rhs[i] += gdown * twall

    cp_work = [0.0] * nbin
    dp = [0.0] * nbin
    cp_work[0] = upper[0] / diag[0]
    dp[0] = rhs[0] / diag[0]
    for i in range(1, nbin):
        denom = diag[i] - lower[i] * cp_work[i - 1]
        cp_work[i] = upper[i] / denom
        dp[i] = (rhs[i] - lower[i] * dp[i - 1]) / denom
    for i in range(nbin - 2, -1, -1):
        dp[i] -= cp_work[i] * dp[i + 1]
    return dp


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate the 1D liquid-conduction reduction using fixed SPARTA surface flux dumps."
    )
    parser.add_argument("case_dir", type=Path, help="Non-conduction case directory containing surf dumps")
    parser.add_argument("--droplet-index", type=int, default=1)
    parser.add_argument("--output", type=Path, default=None, help="Output .dat path")
    parser.add_argument("--average-frames", action="store_true", help="Average all nonzero dump frames instead of using the latest")
    parser.add_argument("--twall", type=float, default=None, help="Wall/liquid base temperature [K]. Defaults to metadata temperature_k")
    parser.add_argument("--latent", type=float, default=DEFAULT_LATENT, help="Latent heat [J/kg]")
    parser.add_argument("--conductivity", type=float, default=DEFAULT_CONDUCTIVITY, help="Liquid conductivity [W/m/K]")
    parser.add_argument("--rho", type=float, default=DEFAULT_RHO, help="Liquid density [kg/m3]")
    parser.add_argument("--cp", type=float, default=DEFAULT_CP, help="Liquid heat capacity [J/kg/K]")
    parser.add_argument("--dtcond", type=float, default=None, help="Conduction update interval for one transient update [s]")
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

    yedges = surface_edges(list(geom_by_id.values()))
    widths = [width_at_y(list(geom_by_id.values()), y) for y in yedges]
    nbin = len(yedges) - 1
    volumes_trap = [0.5 * (widths[i] + widths[i + 1]) * (yedges[i + 1] - yedges[i]) for i in range(nbin)]
    volumes_direct = [direct_slice_area(list(geom_by_id.values()), yedges[i], yedges[i + 1]) for i in range(nbin)]
    heat = [0.0] * nbin
    mass_rate = [0.0] * nbin
    integrated_mass_rate = [0.0] * nbin

    for surf_id, geom in geom_by_id.items():
        ymid = 0.5 * (geom["v1y"] + geom["v2y"])
        ibin = find_bin(yedges, ymid)
        length = segment_length(geom)
        segment_mass_rate = flux_by_id[surf_id] * length
        mass_rate[ibin] += segment_mass_rate
        heat[ibin] += segment_mass_rate * args.latent
        if integrated_by_id:
            integrated_mass_rate[ibin] += integrated_by_id[surf_id]

    steady = solve_steady_temperature(yedges, widths, heat, twall, args.conductivity)
    transient = None
    if args.dtcond is not None:
        transient = solve_transient_temperature(
            yedges,
            widths,
            volumes_trap,
            heat,
            twall,
            twall,
            args.conductivity,
            args.rho,
            args.cp,
            args.dtcond,
        )

    output = args.output
    if output is None:
        output = case_dir / "profiles_steady" / "validate_drop_conduction_1d.dat"
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", encoding="utf-8") as handle:
        handle.write(
            "# y_mid_m ylo_m yhi_m width_lo_m width_hi_m volume_trap_m2 volume_direct_m2 "
            "mass_rate_kg_m_s integrated_mass_rate_kg_m_s heat_W_m heat_per_volume_W_m3 "
            "steady_temperature_K steady_delta_T_K transient_one_update_K transient_one_update_delta_T_K\n"
        )
        for i in range(nbin):
            ymid = 0.5 * (yedges[i] + yedges[i + 1])
            heat_per_volume = heat[i] / volumes_trap[i] if volumes_trap[i] else 0.0
            t_transient = transient[i] if transient is not None else float("nan")
            handle.write(
                f"{ymid:.16e} {yedges[i]:.16e} {yedges[i+1]:.16e} "
                f"{widths[i]:.16e} {widths[i+1]:.16e} "
                f"{volumes_trap[i]:.16e} {volumes_direct[i]:.16e} "
                f"{mass_rate[i]:.16e} {integrated_mass_rate[i]:.16e} "
                f"{heat[i]:.16e} {heat_per_volume:.16e} "
                f"{steady[i]:.16e} {steady[i]-twall:.16e} "
                f"{t_transient:.16e} {t_transient-twall:.16e}\n"
            )

    total_mass = sum(mass_rate)
    total_integrated_mass = sum(integrated_mass_rate)
    total_heat = sum(heat)
    volume_error = max(
        abs(a - b) / max(abs(b), 1.0e-300) for a, b in zip(volumes_trap, volumes_direct)
    )
    mass_error = 0.0
    if integrated_by_id and abs(total_integrated_mass) > 0.0:
        mass_error = abs(total_mass - total_integrated_mass) / abs(total_integrated_mass)

    print(f"case: {case_dir}")
    print(f"flux frame: {frame_label}")
    print(f"bins: {nbin}")
    print(f"total mass rate from mflux*arc_length [kg/m/s]: {total_mass:.8e}")
    if integrated_by_id:
        print(f"total mass rate from integrated dump column [kg/m/s]: {total_integrated_mass:.8e}")
        print(f"relative mass-rate mismatch: {mass_error:.8e}")
    print(f"total latent heat input [W/m]: {total_heat:.8e}")
    print(f"liquid cross-section area from 1D slices [m2/m]: {sum(volumes_trap):.8e}")
    print(f"max relative slice-volume trapezoid/direct mismatch: {volume_error:.8e}")
    print(f"steady Tmin/Tmax [K]: {min(steady):.6f} {max(steady):.6f}")
    print(f"steady max rise [K]: {max(steady) - twall:.6f}")
    if transient is not None:
        print(f"one-update transient Tmin/Tmax [K]: {min(transient):.6f} {max(transient):.6f}")
        print(f"one-update transient max rise [K]: {max(transient) - twall:.6f}")
    print(f"wrote: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
