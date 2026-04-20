#!/usr/bin/env bash

set -euo pipefail

study_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
case_list="${1:-$study_root/cases/case_list.txt}"
sparta_bin="${2:-${SPARTA_BIN:-../../src/spa_serial}}"

if [[ ! -f "$case_list" ]]; then
  echo "Case list not found: $case_list" >&2
  exit 1
fi

while IFS= read -r case_relpath; do
  [[ -z "$case_relpath" ]] && continue
  echo "Running $case_relpath"
  "$study_root/scripts/run_case.sh" "$study_root/$case_relpath" "$sparta_bin"
done < "$case_list"
