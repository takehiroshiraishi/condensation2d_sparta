# Condensation2D Study Platform

This directory is the 2D analogue of `run/condensation`. It keeps the same study-platform structure and boundary-condition logic, but the droplet is represented by a 2D circular-cap cross-section ("kamaboko") instead of the 3D spherical-cap surface mesh.

## Layout

- `base/`: reusable template and shared species/VSS files
- `cases/`: generated self-contained case directories plus study manifests
- `scripts/`: case generation and run helpers
- `post/`: post-processing and ParaView export helpers
- `results/`: summary CSVs and plots
- `test/`: preserved copied reference area from the original condensation folder

## Implemented Geometry Modes

- `single_open_half`: half-domain single-droplet reference with symmetry at `x=0` and an open outer boundary at `x=xhi`
- `array_half`: half-domain representative neighboring-droplet cell with symmetry at both lateral faces

In both modes:

- `y=0` is the wall
- `y=yhi` is the `outletv` top boundary
- the droplet is a 2D circular arc generated per case

## Generate Cases

```bash
python3 scripts/generate_cases.py --config cases/neighbor_condensation_2d/parameters.json --force
```

If there is exactly one `cases/*/parameters.json` file, `--config` is optional.

Generated cases live under:

- `cases/<study_name>/<case_name>/`

Each case receives:

- `in.condensation`
- `metadata.json`
- `pv_grid.txt`
- `run_single.sh`
- `sdata.droplet`
- `water.species`
- `water.vss`

The generator also refreshes:

- `cases/case_manifest.csv`
- `cases/case_list.txt`

## Run One Case

```bash
scripts/run_case.sh cases/neighbor_condensation_2d/single_open_half_vtop_m1 ../../src/spa_serial
```

MPI example:

```bash
SPARTA_LAUNCH="mpirun -np 4" scripts/run_case.sh cases/neighbor_condensation_2d/array_half_dx_8e-06_vtop_m1 ../../src/spa_mpi
```

Per-case Slurm launcher:

```bash
cd cases/neighbor_condensation_2d/single_open_half_vtop_m1
sbatch run_single.sh
```

`run_single.sh` looks for the executable in this order:

- `SPARTA_BIN` if you export it before `sbatch`
- `./spa_mpi` inside the case directory
- one shared `spa_mpi` at the study root: `run/condensation2d/spa_mpi`

For HPC use, the cleanest setup is usually one shared `spa_mpi` at `run/condensation2d/spa_mpi` or setting `SPARTA_BIN=/path/to/spa_mpi`, rather than creating a separate symlink in every case directory.

## Run a Sweep

```bash
scripts/run_sweep.sh cases/case_list.txt ../../src/spa_serial
```

The Slurm array template is:

- `scripts/slurm_array_template.sh`

## Retrieve Results From HPC

To pull case results back from HPC without deleting local-only files such as
ParaView exports, use:

```bash
scripts/sync_from_hpc.sh user@hpc:/path/to/condensation2d
```

This syncs:

- remote `cases/neighbor_condensation_2d/`
- into local `cases/neighbor_condensation_2d/`

It intentionally does not use `rsync --delete`, so local files created only on
your machine are preserved.

If your study name differs, pass it explicitly:

```bash
scripts/sync_from_hpc.sh user@hpc:/path/to/condensation2d my_study_name
```

## Post-Process Results

```bash
python3 post/process_results.py
```

Generated cases write:

- `surf_droplet*.dump`
- `surf_geom_droplet*.dump`
- `grid_steady.dump`
- `line_x0.dump`

The 2D condensation metrics are interpreted per unit out-of-plane depth.

## ParaView Export

Last frame only:

```bash
python3 post/export_paraview_vtk.py cases/neighbor_condensation_2d/single_open_half_vtop_m1
paraview cases/neighbor_condensation_2d/single_open_half_vtop_m1/grid_steady_legacy.vtk
```

Time series:

```bash
python3 post/export_paraview_vtk.py --mode all cases/neighbor_condensation_2d/single_open_half_vtop_m1
paraview cases/neighbor_condensation_2d/single_open_half_vtop_m1/vtk_series/grid_steady/grid_steady.pvd
```

The centerline export is:

- `line_x0.dump`
- `vtk_series/line_x0/line_x0.pvd`

## Which Files To Edit

- Edit `cases/<study_name>/parameters.json` to change the sweep, `cell_size`, droplet radius, contact angle, run length, and output cadence.
- Edit `base/in.condensation.template` to change the SPARTA calculation itself.
- Edit `scripts/generate_cases.py` to change geometry generation, case naming, symmetry logic, or diagnostics.
- Edit `post/process_results.py` to change derived metrics, normalization, CSV columns, or plots.
- Edit `post/export_paraview_vtk.py` if you want different visualization export behavior.

## Notes and Assumptions

- The droplet surface is a generated 2D circular arc, not a closed circle and not a 3D projection.
- The boundary-condition logic matches the current 3D study setup:
  - droplet emits via `emit/surf`
  - droplet uses `evapref`
  - top boundary uses `outletv`
- The grid is controlled by a single isotropic `cell_size`. The generator requires each domain length to be an integer multiple of `cell_size`.
- Surface quantities from `dump surf` are line-based in 2D, so the reported condensation metrics are per unit out-of-plane depth.
- `single_open_half` reconstructs full-droplet metrics with a symmetry multiplier of 2.
- `array_half` also uses a symmetry multiplier of 2 for the representative neighboring-droplet cell.
