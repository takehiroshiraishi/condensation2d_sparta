#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 CASE_DIR [SPARTA_BIN]" >&2
  exit 1
fi

case_dir="$1"
sparta_bin="${2:-${SPARTA_BIN:-../../src/spa_serial}}"
sparta_launch="${SPARTA_LAUNCH:-}"

if [[ ! -d "$case_dir" ]]; then
  echo "Case directory not found: $case_dir" >&2
  exit 1
fi

if [[ ! -f "$case_dir/in.condensation" ]]; then
  echo "Missing in.condensation in $case_dir" >&2
  exit 1
fi

pushd "$case_dir" >/dev/null

if [[ -n "$sparta_launch" ]]; then
  read -r -a launch_parts <<< "$sparta_launch"
  "${launch_parts[@]}" "$sparta_bin" < in.condensation > log.sparta 2>&1
else
  "$sparta_bin" < in.condensation > log.sparta 2>&1
fi

popd >/dev/null
