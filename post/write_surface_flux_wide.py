#!/usr/bin/env python3

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path


def format_um(xhi_m: float) -> str:
    value_um = xhi_m * 1.0e6
    if abs(value_um - round(value_um)) < 1.0e-9:
        return str(int(round(value_um)))
    return f"{value_um:g}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write a wide surface-flux table from condensation_flux_summary.dat."
    )
    parser.add_argument("study_dir", type=Path, help="Study directory under cases/")
    parser.add_argument(
        "--summary",
        default="condensation_flux_summary.dat",
        help="Input summary filename within the study directory",
    )
    parser.add_argument(
        "--output",
        default="condensation_flux_surface_wide.dat",
        help="Output filename within the study directory",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    study_dir = args.study_dir.resolve()
    summary_path = study_dir / args.summary
    if not summary_path.exists():
        raise FileNotFoundError(f"Summary file not found: {summary_path}")

    with summary_path.open("r", encoding="utf-8") as handle:
        header = handle.readline().split()
        rows = []
        for line in handle:
            parts = line.split()
            if not parts:
                continue
            row = dict(zip(header, parts))
            rows.append(row)

    grouped: dict[float, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[float(row["xhi_m"])].append(row)

    xhis = sorted(grouped)
    if not xhis:
        raise ValueError(f"No rows found in {summary_path}")

    for xhi in xhis:
        grouped[xhi].sort(key=lambda row: float(row["equilibrium_pressure_Pa"]))

    row_count = len(grouped[xhis[0]])
    if any(len(grouped[xhi]) != row_count for xhi in xhis):
        raise ValueError("Each xhi group must have the same number of rows to build a wide table")

    fieldnames: list[str] = []
    for xhi in xhis:
        suffix = format_um(xhi)
        fieldnames.extend(
            [
                f"pressure_at_y_300um_xhi_{suffix}um_Pa",
                f"mass_flux_per_surface_length_xhi_{suffix}um",
                f"normalized_flux_per_surface_length_xhi_{suffix}um",
            ]
        )

    output_path = study_dir / args.output
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(" ".join(fieldnames) + "\n")
        for index in range(row_count):
            values: list[str] = []
            for xhi in xhis:
                row = grouped[xhi][index]
                values.extend(
                    [
                        row["pressure_at_y_300um_Pa"],
                        row["mass_flux_per_surface_length"],
                        row["normalized_flux_per_surface_length"],
                    ]
                )
            handle.write(" ".join(values) + "\n")

    print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
