# Repository Guidelines

## Project Structure & Module Organization
This repository is a SPARTA DSMC codebase. Core C++ sources live in `src/`; machine-specific build settings are in `src/MAKE/`. Shared physics data and sample surfaces live in `data/`. Reference problems are under `examples/` and `bench/`, while custom run folders such as `run/condensation/test/` hold input decks like `in.condensation` and `in.evapChannel` plus generated outputs (`log.sparta`, `image_*.ppm`, `video.gif`).

## Build, Test, and Development Commands
Build from `src/` with the provided make targets:

```bash
cd src
make serial        # builds spa_serial
make mpi           # builds spa_mpi
make clean-all     # removes Obj_* build directories
```

Run a case from its working directory so relative paths resolve correctly:

```bash
../../src/spa_mpi < in.condensation
../../src/spa_serial < in.evapChannel
```

Use `examples/README` as the reference pattern for quick smoke tests and visualization.

## Coding Style & Naming Conventions
Match the existing C++ style in `src/`: two-space indentation, opening braces on the next line for functions, and no unnecessary alignment changes. File names use lowercase snake_case, typically in paired forms such as `compute_grid.cpp` and `compute_grid.h`. New input scripts should keep the `in.*` naming pattern and group commands in a stable order: variables, geometry, physics, output, then `run`.

## Testing Guidelines
There is no separate unit-test framework in this checkout. Validate changes by rebuilding the affected target and running at least one representative example or run case. Compare new `log.sparta` output against the archived `examples/*/log.*` or `bench/log.*` files for statistically similar behavior, not bitwise identity. Keep generated artifacts out of commits unless they are intentional fixtures.

## Commit & Pull Request Guidelines
This workspace does not include `.git` history, so use short imperative commit subjects under 72 characters, for example `Add condensation test input` or `Fix grid balance regression`. PRs should describe the physical/modeling change, list the exact run commands used for verification, and attach screenshots or GIFs only when dump/image behavior changes.

## Configuration Tips
Input decks often reference shared assets through relative paths such as `../../data/air.vss`; preserve that layout when adding new cases. Avoid hard-coding local viewer commands unless they are optional post-processing steps.
