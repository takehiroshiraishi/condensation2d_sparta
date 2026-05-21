#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path


def parse_last_surf_frame(path: Path) -> tuple[list[str], list[dict[str, float]]]:
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
                row = {column: float(value) for column, value in zip(columns, raw)}
                rows.append(row)
            last_columns = columns
            last_rows = rows

    if last_columns is None or last_rows is None:
        raise ValueError(f"No frames found in {path}")
    return last_columns, last_rows


def main() -> int:
    parser = argparse.ArgumentParser(description="Write local condensation flux profile along the droplet surface.")
    parser.add_argument("case_dir", type=Path, help="Case directory")
    parser.add_argument("--droplet-index", type=int, default=1, help="Droplet index to process")
    args = parser.parse_args()

    case_dir = args.case_dir.resolve()
    droplet_index = args.droplet_index
    flux_path = case_dir / f"surf_droplet{droplet_index}.dump"
    geom_path = case_dir / f"surf_geom_droplet{droplet_index}.dump"
    if not flux_path.exists():
        raise FileNotFoundError(f"Missing surface flux dump: {flux_path}")
    if not geom_path.exists():
        raise FileNotFoundError(f"Missing surface geometry dump: {geom_path}")

    flux_columns, flux_rows = parse_last_surf_frame(flux_path)
    geom_columns, geom_rows = parse_last_surf_frame(geom_path)
    flux_column = f"f_avg_droplet_{droplet_index}[1]"
    if flux_column not in flux_columns:
        raise KeyError(f"Missing {flux_column} in {flux_path}")

    flux_by_id = {int(row["id"]): row[flux_column] for row in flux_rows}
    geom_by_id = {int(row["id"]): row for row in geom_rows}
    if flux_by_id.keys() != geom_by_id.keys():
        raise ValueError(f"Surface IDs do not match between {flux_path} and {geom_path}")

    segments: list[dict[str, float]] = []
    for surf_id, geom in geom_by_id.items():
        x_mid = 0.5 * (geom["v1x"] + geom["v2x"])
        y_mid = 0.5 * (geom["v1y"] + geom["v2y"])
        seg_len = ((geom["v2x"] - geom["v1x"]) ** 2 + (geom["v2y"] - geom["v1y"]) ** 2) ** 0.5
        segments.append(
            {
                "id": surf_id,
                "x_mid": x_mid,
                "y_mid": y_mid,
                "seg_len": seg_len,
                "flux": flux_by_id[surf_id],
            }
        )

    # Half-droplet surfaces run from apex (smallest x) to contact line (largest x).
    segments.sort(key=lambda row: (row["x_mid"], -row["y_mid"]))

    distance = 0.0
    output_rows = []
    for index, row in enumerate(segments):
        if index == 0:
            distance = 0.5 * row["seg_len"]
        else:
            distance += 0.5 * segments[index - 1]["seg_len"] + 0.5 * row["seg_len"]
        output_rows.append(
            {
                "arc_dist_m": distance,
                "local_mass_flux": row["flux"],
                "x_mid": row["x_mid"],
                "y_mid": row["y_mid"],
            }
        )

    output_dir = case_dir / "profiles_steady"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "local_flux.dat"
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("arc_dist_m local_mass_flux x_mid y_mid\n")
        for row in output_rows:
            handle.write(
                f"{row['arc_dist_m']} {row['local_mass_flux']} {row['x_mid']} {row['y_mid']}\n"
            )

    print(f"Wrote local flux profile to: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
