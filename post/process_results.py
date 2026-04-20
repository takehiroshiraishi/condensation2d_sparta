#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path


STUDY_ROOT = Path(__file__).resolve().parents[1]
CASES_DIR = STUDY_ROOT / "cases"
RESULTS_DIR = STUDY_ROOT / "results"
PLOTS_DIR = RESULTS_DIR / "plots"
MANIFEST_PATH = CASES_DIR / "case_manifest.csv"


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_manifest_case_dirs() -> list[Path]:
    if not MANIFEST_PATH.exists():
        return sorted(metadata_path.parent for metadata_path in CASES_DIR.glob("*/metadata.json"))

    case_dirs: list[Path] = []
    with MANIFEST_PATH.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            case_relpath = row.get("case_relpath", "")
            if not case_relpath:
                continue
            case_dirs.append(STUDY_ROOT / case_relpath)
    return case_dirs


def parse_dump_frames(path: Path) -> list[dict]:
    frames: list[dict] = []
    with path.open("r", encoding="utf-8") as handle:
        lines = iter(handle)
        for line in lines:
            if not line:
                continue
            if line.strip() != "ITEM: TIMESTEP":
                continue
            timestep = int(next(lines).strip())
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
            frames.append({"timestep": timestep, "columns": columns, "rows": rows})
    return frames


def analyze_case(case_dir: Path) -> dict | None:
    metadata_path = case_dir / "metadata.json"
    if not metadata_path.exists():
        return None

    metadata = load_json(metadata_path)
    droplet_count = metadata["droplet_count"]
    latent_heat = metadata["latent_heat_j_per_kg"]

    top_boundary_velocity = metadata.get("top_boundary_velocity", metadata.get("inlet_velocity"))

    row: dict[str, object] = {
        "case_name": metadata["case_name"],
        "geometry_mode": metadata["geometry_mode"],
        "top_boundary_velocity": top_boundary_velocity,
        "spacing": metadata["spacing"],
        "droplet_count": droplet_count,
        "status": "complete",
    }

    total_area = 0.0
    total_mass_rate = 0.0
    total_mass_rate_from_flux = 0.0

    for droplet_index, droplet_meta in enumerate(metadata["droplets"], start=1):
        dump_path = case_dir / f"surf_droplet{droplet_index}.dump"
        if not dump_path.exists():
            row["status"] = "missing_surface_dump"
            return row

        frames = parse_dump_frames(dump_path)
        if not frames:
            row["status"] = "empty_surface_dump"
            return row

        last_frame = frames[-1]
        area = sum(entry["area"] for entry in last_frame["rows"])
        area *= droplet_meta["symmetry_multiplier"]
        flux_column = f"f_avg_droplet_{droplet_index}[1]"
        flow_column = f"f_avg_droplet_{droplet_index}[2]"

        area_weighted_flux = sum(entry["area"] * entry[flux_column] for entry in last_frame["rows"])
        mass_rate_from_flux = area_weighted_flux * droplet_meta["symmetry_multiplier"]
        mass_rate = sum(entry[flow_column] for entry in last_frame["rows"])
        mass_rate *= droplet_meta["symmetry_multiplier"]

        heat_rate = mass_rate * latent_heat
        heat_flux = heat_rate / area if area else 0.0

        row[f"droplet_{droplet_index}_area"] = area
        row[f"droplet_{droplet_index}_mass_rate"] = mass_rate
        row[f"droplet_{droplet_index}_mass_rate_from_flux"] = mass_rate_from_flux
        row[f"droplet_{droplet_index}_heat_rate"] = heat_rate
        row[f"droplet_{droplet_index}_heat_flux"] = heat_flux

        total_area += area
        total_mass_rate += mass_rate
        total_mass_rate_from_flux += mass_rate_from_flux

    row["total_area"] = total_area
    row["total_mass_rate"] = total_mass_rate
    row["total_mass_rate_from_flux"] = total_mass_rate_from_flux
    row["mean_mass_rate_per_droplet"] = total_mass_rate / droplet_count if droplet_count else 0.0
    row["mean_heat_rate_per_droplet"] = row["mean_mass_rate_per_droplet"] * latent_heat if droplet_count else 0.0
    row["mean_heat_flux_per_droplet"] = (
        row["mean_heat_rate_per_droplet"] / (total_area / droplet_count) if droplet_count and total_area else 0.0
    )
    row["positive_condensation_mass_rate"] = max(total_mass_rate, 0.0)
    row["positive_condensation_heat_flux"] = max(row["mean_heat_flux_per_droplet"], 0.0)
    return row


def attach_single_reference_comparisons(rows: list[dict]) -> None:
    single_reference = {}
    for row in rows:
        if row.get("status") != "complete":
            continue
        if row["geometry_mode"] not in {"single_full", "single_quarter", "single_open_quarter", "single_open_half"}:
            continue
        key = (row["top_boundary_velocity"],)
        single_reference[key] = row

    for row in rows:
        if row.get("status") != "complete":
            row["normalized_heat_flux_vs_single"] = ""
            row["normalized_mass_rate_vs_single"] = ""
            continue

        reference = single_reference.get((row["top_boundary_velocity"],))
        if reference is None:
            row["normalized_heat_flux_vs_single"] = ""
            row["normalized_mass_rate_vs_single"] = ""
            continue

        ref_heat_flux = reference["mean_heat_flux_per_droplet"]
        ref_mass_rate = reference["mean_mass_rate_per_droplet"]
        row["normalized_heat_flux_vs_single"] = row["mean_heat_flux_per_droplet"] / ref_heat_flux if ref_heat_flux else ""
        row["normalized_mass_rate_vs_single"] = row["mean_mass_rate_per_droplet"] / ref_mass_rate if ref_mass_rate else ""


def write_summary(rows: list[dict]) -> Path:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    summary_path = RESULTS_DIR / "summary.csv"
    fieldnames = sorted({key for row in rows for key in row.keys()})
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return summary_path


def make_plots(rows: list[dict]) -> None:
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        print("matplotlib not available; CSV summary was written but plots were skipped.")
        return

    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    pair_rows = [
        row
        for row in rows
        if row.get("status") == "complete" and row["geometry_mode"] in {"pair_x", "array_quarter"}
    ]
    by_velocity = defaultdict(list)
    for row in pair_rows:
        if row.get("normalized_heat_flux_vs_single") == "":
            continue
        by_velocity[row["top_boundary_velocity"]].append(row)

    if by_velocity:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for top_boundary_velocity, group in sorted(by_velocity.items()):
            group = sorted(group, key=lambda item: item["spacing"])
            ax.plot(
                [item["spacing"] for item in group],
                [item["normalized_heat_flux_vs_single"] for item in group],
                marker="o",
                label=f"vtop={top_boundary_velocity:g}",
            )
        ax.axhline(1.0, color="black", linestyle="--", linewidth=1.0)
        ax.set_xlabel("Droplet spacing")
        ax.set_ylabel("Neighbored metric / single metric")
        ax.set_title("Normalized Condensation Metric vs Spacing")
        ax.legend()
        fig.tight_layout()
        fig.savefig(PLOTS_DIR / "normalized_heat_flux_vs_spacing.png", dpi=150)
        plt.close(fig)

    velocity_groups = defaultdict(list)
    for row in rows:
        if row.get("status") != "complete":
            continue
        if row["geometry_mode"] == "single_full":
            label = "single_full"
        elif row["geometry_mode"] == "single_open_quarter":
            label = "single_open_quarter"
        elif row["geometry_mode"] == "single_quarter":
            label = "single_quarter"
        elif row["geometry_mode"] == "single_open_half":
            label = "single_open_half"
        elif row["geometry_mode"] == "pair_x":
            label = f"pair_x_dx={row['spacing']:g}"
        elif row["geometry_mode"] == "array_quarter":
            label = f"array_quarter_dx={row['spacing']:g}"
        elif row["geometry_mode"] == "array_half":
            label = f"array_half_dx={row['spacing']:g}"
        else:
            label = row["geometry_mode"]
        velocity_groups[label].append(row)

    if velocity_groups:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for label, group in sorted(velocity_groups.items()):
            group = sorted(group, key=lambda item: item["top_boundary_velocity"])
            ax.plot(
                [item["top_boundary_velocity"] for item in group],
                [item["mean_heat_flux_per_droplet"] for item in group],
                marker="o",
                label=label,
            )
        ax.set_xlabel("Top boundary velocity")
        ax.set_ylabel("Heat-flux-like metric per droplet")
        ax.set_title("Condensation Metric vs Top Boundary Velocity")
        ax.legend()
        fig.tight_layout()
        fig.savefig(PLOTS_DIR / "heat_flux_vs_top_boundary_velocity.png", dpi=150)
        plt.close(fig)


def main() -> int:
    rows = []
    for case_dir in load_manifest_case_dirs():
        analyzed = analyze_case(case_dir)
        if analyzed is not None:
            rows.append(analyzed)

    attach_single_reference_comparisons(rows)
    summary_path = write_summary(rows)
    make_plots(rows)
    print(f"Wrote {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
