#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
from pathlib import Path


STUDY_ROOT = Path(__file__).resolve().parents[1]
BASE_DIR = STUDY_ROOT / "base"
CASES_DIR = STUDY_ROOT / "cases"
TEMPLATE_PATH = BASE_DIR / "in.flat_surface.template"
DEFAULT_CONFIG_NAME = "parameters.json"
ASSET_FILES = ("water.species", "water.vss")


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def dump_json(path: Path, data: dict) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
        handle.write("\n")


def study_name_from_config_path(config_path: Path) -> str:
    if config_path.name != DEFAULT_CONFIG_NAME:
        raise ValueError(f"Expected config file named {DEFAULT_CONFIG_NAME}: {config_path}")
    if config_path.parent == CASES_DIR or config_path.parent.name == "_templates":
        raise ValueError(f"Config file must live under cases/<study_name>/: {config_path}")
    return config_path.parent.name


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def format_float(value: float) -> str:
    return f"{value:.12g}"


def slug_float(value: float) -> str:
    sign = "m" if value < 0.0 else ""
    text = f"{abs(value):.6g}"
    return sign + text.replace(".", "p").replace("+", "")


def compute_grid_count(length: float, cell_size: float, axis_name: str) -> int:
    raw_count = length / cell_size
    rounded_count = round(raw_count)
    require(
        rounded_count > 0 and math.isclose(raw_count, rounded_count, rel_tol=0.0, abs_tol=1.0e-9),
        f"{axis_name} length {format_float(length)} is not an integer multiple of cell_size {format_float(cell_size)}.",
    )
    return int(rounded_count)


def boundary_state_pairs(config: dict) -> list[tuple[float, float]]:
    cases_cfg = config.get("cases", {})
    sweep = config["sweep"]
    top_temperatures = sweep["top_boundary_temperature_k"]
    top_number_densities = sweep["top_boundary_number_density"]
    mode = cases_cfg.get("boundary_state_mode", "paired")

    if mode == "cartesian":
        return [
            (top_temperature, top_number_density)
            for top_temperature in top_temperatures
            for top_number_density in top_number_densities
        ]

    if mode != "paired":
        raise ValueError(f"Unsupported boundary_state_mode: {mode}")

    if len(top_temperatures) == len(top_number_densities):
        return list(zip(top_temperatures, top_number_densities))
    if len(top_temperatures) == 1:
        return [(top_temperatures[0], top_number_density) for top_number_density in top_number_densities]
    if len(top_number_densities) == 1:
        return [(top_temperature, top_number_densities[0]) for top_temperature in top_temperatures]

    raise ValueError(
        "For boundary_state_mode='paired', top_boundary_temperature_k and "
        "top_boundary_number_density must have the same length, or one list must have length 1."
    )


def build_geometry(defaults: dict) -> dict:
    box_length_x = defaults["box_length_x"]
    box_height = defaults["box_height"]
    cell_size = defaults["cell_size"]
    xlo = -0.5 * box_length_x
    xhi = 0.5 * box_length_x
    ylo = 0.0
    yhi = box_height
    zlo = -0.5 * cell_size
    zhi = 0.5 * cell_size
    grid_nx = compute_grid_count(xhi - xlo, cell_size, "x")
    grid_ny = compute_grid_count(yhi - ylo, cell_size, "y")
    return {
        "xlo": xlo,
        "xhi": xhi,
        "ylo": ylo,
        "yhi": yhi,
        "zlo": zlo,
        "zhi": zhi,
        "grid_nx": grid_nx,
        "grid_ny": grid_ny,
        "cell_size": cell_size,
    }


def render_pre_run_section(defaults: dict) -> str:
    start_sampling_step = defaults["start_sampling_step"]
    if start_sampling_step <= 0:
        return ""
    return f"# Equilibration run before any steady-state averaging begins.\nrun                 {start_sampling_step}"


def render_centerline_section(defaults: dict, geometry: dict) -> str:
    dump_every = defaults["sampling_steps"]
    dx = (geometry["xhi"] - geometry["xlo"]) / geometry["grid_nx"]
    band_cells = defaults["centerline_band_cells"]
    half_thickness = 0.500001 * band_cells * dx
    line_xlo = -half_thickness
    line_xhi = half_thickness
    return "\n".join(
        [
            "# Time-averaged centerline fields in a thin band around x = 0.",
            "# f_centerline_avg columns are: nrho u v trot temp press.",
            f"region              xline block {format_float(line_xlo)} {format_float(line_xhi)} INF INF INF INF",
            "group               centerline grid region xline center",
            "compute             centerline_flow grid centerline water nrho u v trot",
            "compute             centerline_thermo thermal/grid centerline water temp press",
            f"fix                 centerline_avg ave/grid centerline 1 {dump_every} {dump_every} &",
            "                    c_centerline_flow[*] c_centerline_thermo[*] ave one",
            f"dump                centerline_dump grid centerline {dump_every} line_x0.dump id xc yc f_centerline_avg[*]",
        ]
    )


def render_full_grid_section(defaults: dict) -> str:
    dump_every = defaults["sampling_steps"]
    return "\n".join(
        [
            "# Time-averaged whole-domain grid fields for steady-state inspection.",
            "# f_grid_avg columns are: nrho u v trot temp press.",
            "compute             grid_flow grid all water nrho u v trot",
            "compute             grid_thermo thermal/grid all water temp press",
            f"fix                 grid_avg ave/grid all 1 {dump_every} {dump_every} &",
            "                    c_grid_flow[*] c_grid_thermo[*] ave one",
            f"dump                grid_dump grid all {dump_every} grid_steady.dump id xc yc f_grid_avg[*]",
        ]
    )


def render_diagnostics_section(defaults: dict, geometry: dict) -> str:
    return "\n".join([render_centerline_section(defaults, geometry), "", render_full_grid_section(defaults)])


def render_case_input(template_text: str, defaults: dict, geometry: dict, top_temperature: float, top_number_density: float) -> str:
    production_run_steps = defaults["run_steps"] - defaults["start_sampling_step"]
    placeholders = {
        "__SEED__": str(defaults["seed"]),
        "__GRIDCUT__": format_float(defaults["gridcut"]),
        "__XLO__": format_float(geometry["xlo"]),
        "__XHI__": format_float(geometry["xhi"]),
        "__YLO__": format_float(geometry["ylo"]),
        "__YHI__": format_float(geometry["yhi"]),
        "__ZLO__": format_float(geometry["zlo"]),
        "__ZHI__": format_float(geometry["zhi"]),
        "__TWALL__": format_float(defaults["temperature_k"]),
        "__COEFF__": format_float(defaults["condensation_coefficient"]),
        "__VAPOR_NRHO__": format_float(defaults["vapor_number_density"]),
        "__TOP_TEMP__": format_float(top_temperature),
        "__TOP_NRHO__": format_float(top_number_density),
        "__GRID_NX__": str(geometry["grid_nx"]),
        "__GRID_NY__": str(geometry["grid_ny"]),
        "__FNUM__": format_float(defaults["fnum"]),
        "__TIMESTEP__": format_float(defaults["time_step"]),
        "__PRE_RUN_SECTION__": render_pre_run_section(defaults),
        "__DIAGNOSTICS_SECTION__": render_diagnostics_section(defaults, geometry),
        "__STATS_EVERY__": str(defaults["stats_every"]),
        "__PRODUCTION_RUN_STEPS__": str(production_run_steps),
    }
    rendered = template_text
    for key, value in placeholders.items():
        rendered = rendered.replace(key, value)
    return rendered


def render_paraview_grid_description(geometry: dict) -> str:
    return "\n".join(
        [
            "# Case-local grid description for tools/paraview/grid2paraview.py",
            "dimension           2",
            f"create_box          {format_float(geometry['xlo'])} {format_float(geometry['xhi'])} {format_float(geometry['ylo'])} {format_float(geometry['yhi'])} {format_float(geometry['zlo'])} {format_float(geometry['zhi'])}",
            f"create_grid         {geometry['grid_nx']} {geometry['grid_ny']} 1",
        ]
    )


def render_run_single_script() -> str:
    return "\n".join(
        [
            "#!/bin/bash",
            "#SBATCH --nodes 1",
            "#SBATCH --job-name flatcond2d",
            "#SBATCH --ntasks-per-node 32",
            "#SBATCH --cpus-per-task 1",
            "#SBATCH --time 7-00:00:00",
            "",
            "set -euo pipefail",
            "",
            'case_dir="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"',
            'cases_dir="$(cd "$case_dir/../.." && pwd)"',
            'sparta_bin="${SPARTA_BIN:-$cases_dir/spa_mpi}"',
            "",
            "module purge",
            "module load gcc",
            "module load openmpi",
            "",
            'if [[ ! -x "$sparta_bin" ]]; then',
            '  echo "Could not find spa_mpi at $sparta_bin." >&2',
            '  echo "Place spa_mpi at condensation2d/cases/spa_mpi or set SPARTA_BIN." >&2',
            "  exit 1",
            "fi",
            "",
            'cd "$case_dir"',
            'srun "$sparta_bin" < in.condensation > log.txt',
        ]
    )


def render_study_run_script() -> str:
    return "\n".join(
        [
            "#!/usr/bin/env bash",
            "",
            "set -euo pipefail",
            "",
            'study_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
            'case_list="$study_dir/case_list.txt"',
            "",
            'while IFS= read -r case_relpath; do',
            '  [[ -z "$case_relpath" ]] && continue',
            '  case_dir="$study_dir/$(basename "$case_relpath")"',
            '  echo "Submitting $case_dir/run_single.sh"',
            '  (cd "$case_dir" && sbatch ./run_single.sh)',
            'done < "$case_list"',
        ]
    )


def render_study_profiles_script() -> str:
    return "\n".join(
        [
            "#!/usr/bin/env bash",
            "",
            "set -euo pipefail",
            "",
            'study_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
            'study_root="$(cd "$study_dir/../.." && pwd)"',
            'case_list="$study_dir/case_list.txt"',
            "",
            'while IFS= read -r case_relpath; do',
            '  [[ -z "$case_relpath" ]] && continue',
            '  case_dir="$study_dir/$(basename "$case_relpath")"',
            '  echo "Converting steady outputs to VTK for $case_dir"',
            '  python3 "$study_root/post/export_paraview_vtk.py" --mode all "$case_dir"',
            'done < "$case_list"',
            "",
            'echo "Flat-surface reference cases only write grid and centerline outputs."',
        ]
    )


def generate_cases(config_path: Path, force: bool) -> list[dict]:
    config = load_json(config_path)
    defaults = config["defaults"]
    require(defaults["start_sampling_step"] >= 0, "start_sampling_step must be non-negative.")
    require(defaults["sampling_steps"] > 0, "sampling_steps must be positive.")
    require(defaults["run_steps"] > defaults["start_sampling_step"], "run_steps must exceed start_sampling_step.")
    require(
        (defaults["run_steps"] - defaults["start_sampling_step"]) % defaults["sampling_steps"] == 0,
        "run_steps - start_sampling_step must be divisible by sampling_steps.",
    )
    template_text = TEMPLATE_PATH.read_text(encoding="utf-8")
    study_name = study_name_from_config_path(config_path)
    study_dir = CASES_DIR / study_name
    study_dir.mkdir(parents=True, exist_ok=True)
    geometry = build_geometry(defaults)
    manifest_rows: list[dict] = []

    for top_temperature, top_number_density in boundary_state_pairs(config):
        case_name = f"flat_ttop_{slug_float(top_temperature)}_ntop_{slug_float(top_number_density)}"
        case_dir = study_dir / case_name
        if case_dir.exists():
            if not force:
                raise FileExistsError(f"{case_dir} already exists. Re-run with --force to replace it.")
            shutil.rmtree(case_dir)
        case_dir.mkdir(parents=True, exist_ok=True)

        rendered_input = render_case_input(template_text, defaults, geometry, top_temperature, top_number_density)
        (case_dir / "in.condensation").write_text(rendered_input + "\n", encoding="utf-8")
        (case_dir / "pv_grid.txt").write_text(render_paraview_grid_description(geometry) + "\n", encoding="utf-8")
        run_single_path = case_dir / "run_single.sh"
        run_single_path.write_text(render_run_single_script() + "\n", encoding="utf-8")
        run_single_path.chmod(0o755)
        for asset_name in ASSET_FILES:
            shutil.copy2(BASE_DIR / asset_name, case_dir / asset_name)

        dump_json(
            case_dir / "metadata.json",
            {
                "case_name": case_name,
                "study_name": study_name,
                "geometry_mode": "flat_surface",
                "top_boundary_temperature_k": top_temperature,
                "top_boundary_number_density": top_number_density,
                "temperature_k": defaults["temperature_k"],
                "condensation_coefficient": defaults["condensation_coefficient"],
                "vapor_number_density": defaults["vapor_number_density"],
                "cell_size": defaults["cell_size"],
                "box_length_x": defaults["box_length_x"],
                "box_height": defaults["box_height"],
                "grid_cells": [geometry["grid_nx"], geometry["grid_ny"], 1],
                "simulation_bounds": {
                    "xlo": geometry["xlo"],
                    "xhi": geometry["xhi"],
                    "ylo": geometry["ylo"],
                    "yhi": geometry["yhi"],
                    "zlo": geometry["zlo"],
                    "zhi": geometry["zhi"],
                },
                "surface_model": "evaprefpart_ylo_full_wall",
            },
        )

        manifest_rows.append(
            {
                "case_name": case_name,
                "case_relpath": f"cases/{study_name}/{case_name}",
                "geometry_mode": "flat_surface",
                "top_boundary_temperature_k": top_temperature,
                "top_boundary_number_density": top_number_density,
            }
        )

    manifest_path = study_dir / "case_manifest.csv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "case_name",
                "case_relpath",
                "geometry_mode",
                "top_boundary_temperature_k",
                "top_boundary_number_density",
            ],
        )
        writer.writeheader()
        for row in manifest_rows:
            writer.writerow(row)

    case_list_path = study_dir / "case_list.txt"
    with case_list_path.open("w", encoding="utf-8") as handle:
        for row in manifest_rows:
            handle.write(f"{row['case_relpath']}\n")

    run_script_path = study_dir / "run.sh"
    run_script_path.write_text(render_study_run_script() + "\n", encoding="utf-8")
    run_script_path.chmod(0o755)
    profiles_script_path = study_dir / "plot_profiles.sh"
    profiles_script_path.write_text(render_study_profiles_script() + "\n", encoding="utf-8")
    profiles_script_path.chmod(0o755)
    return manifest_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate planar flat-surface condensation reference cases.")
    parser.add_argument("--config", type=Path, required=True, help="Path to the JSON sweep configuration.")
    parser.add_argument("--force", action="store_true", help="Replace existing generated case directories.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config_path = args.config.resolve()
    study_dir = CASES_DIR / study_name_from_config_path(config_path)
    rows = generate_cases(config_path, args.force)
    print(f"Generated {len(rows)} flat-surface case directories under {study_dir}")
    print(f"Parameters: {config_path}")
    print(f"Study manifest: {study_dir / 'case_manifest.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
