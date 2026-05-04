# Condensation2D Agent Guide

This directory is a self-contained 2D SPARTA study platform for vapor condensation over a circular-cap droplet cross-section. Treat it as a generator-driven workflow: most case directories under `cases/` are generated outputs, not primary source files.

## Scope

- Active work in this subtree should prefer `run/condensation2d` conventions over older `run/condensation` patterns.
- The authoritative sources for study generation are:
  - `cases/init.sh`
  - `cases/_templates/`
  - `scripts/generate_cases.py`
  - `base/in.condensation.template`

## Directory Roles

- `base/`: shared SPARTA template and common species/VSS data.
- `cases/`: study bootstrap plus generated study directories.
- `cases/_templates/`: tracked templates copied into new studies by `cases/init.sh`.
- `cases/<study_name>/`: generated or user-edited study area. This is where `parameters.json`, `generate.sh`, `run.sh`, `plot_profiles.sh`, and all case directories live.
- `scripts/`: generation and run helpers.
- `post/`: post-processing and VTK/profile export helpers.
- `results/`: derived summaries/plots.

## Source Of Truth

- Do not treat `cases/<study_name>/<case_name>/in.condensation` as the primary place to make broad workflow changes.
- If a change should affect future generated cases, edit:
  - `scripts/generate_cases.py` for generated sections, naming, scripts, manifests, diagnostics, and workflow helpers.
  - `base/in.condensation.template` for template-level SPARTA lines surrounding placeholders.
  - `cases/_templates/parameters.json` for new-study defaults.
  - `cases/_templates/*.gp` or other templates for study bootstrap helper files.
- Only edit a generated case file directly when the user explicitly wants a one-off/manual experiment for that case.

## Study Workflow

1. Create a study with `./cases/init.sh <study_name>`.
2. Edit `cases/<study_name>/parameters.json`.
3. Generate cases with `cases/<study_name>/generate.sh` or `python3 scripts/generate_cases.py --config cases/<study_name>/parameters.json --force`.
4. Submit all cases from the study directory with `./run.sh`.
5. Build VTK/profile outputs with `./plot_profiles.sh`.

## Naming And Config Rules

- Study name is derived from the config path `cases/<study_name>/parameters.json`. It is no longer stored in `parameters.json`.
- Case names are generated from geometry mode plus swept parameters. If a parameter creates distinct geometries, ensure the name changes too to avoid collisions.
- `box_height` may be swept from `parameters.json`; if multiple values are present, case names include an `hbox_...` suffix.

## Boundary Condition Context

- The current default top boundary in generated cases is `emit/face` plus `vanish`, not the older `outletv` setup.
- There are archived/manual experiments in some study directories using `subsonic` pressure boundaries or older outlet variants. Do not generalize from a single case file without checking whether the generator already encodes that behavior.
- If changing the top boundary model for future studies, implement it in the generator/template and preserve the prior mode as an explicit option when the user asks to keep the current setup archived.

## Rotational Temperature / Diagnostics

- Generated diagnostics currently include `trot` in grid and line dumps.
- Post-processing in `post/export_paraview_vtk.py` and `post/plot_steady_profiles.py` assumes specific dump column ordering. If SPARTA dump columns change, update both scripts together.
- `plot_steady_profiles.py` uses VTK frames under `vtk_series/grid_steady/`, not the legacy single-frame `.vtk` files.

## HPC / Run Assumptions

- Generated `run_single.sh` scripts assume `spa_mpi` is available at `cases/spa_mpi` unless `SPARTA_BIN` is set.
- On compute nodes, generated job scripts rely on `SLURM_SUBMIT_DIR` to resolve the actual case directory.
- `run.sh` submits each case’s `run_single.sh`; it is not itself an `sbatch` script.

## Editing Guidance

- Prefer updating generators/templates over mass-editing generated case directories.
- Preserve manual/archive studies in `cases/` unless the user explicitly asks to regenerate or delete them.
- When changing profile output formats, also update any gnuplot templates or helper scripts under `cases/_templates/`.
- Check `README.md`, but do not assume it is perfectly current; prefer the code paths in `scripts/`, `post/`, and `cases/init.sh` when they disagree.
