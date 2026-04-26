#!/usr/bin/env bash

set -euo pipefail

study_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
case_list="${1:-}"
sparta_bin="${2:-${SPARTA_BIN:-../../src/spa_serial}}"

if [[ -z "$case_list" ]]; then
  mapfile -t matches < <(find "$study_root/cases/studies" -mindepth 2 -maxdepth 2 -name case_list.txt | sort)
  if [[ ${#matches[@]} -eq 1 ]]; then
    case_list="${matches[0]}"
  else
    echo "Usage: $0 CASE_LIST [SPARTA_BIN]" >&2
    echo "Pass a study-specific case list under cases/studies/<study_name>/case_list.txt." >&2
    exit 1
  fi
fi

if [[ ! -f "$case_list" ]]; then
  echo "Case list not found: $case_list" >&2
  exit 1
fi

while IFS= read -r case_relpath; do
  [[ -z "$case_relpath" ]] && continue
  echo "Running $case_relpath"
  "$study_root/scripts/run_case.sh" "$study_root/$case_relpath" "$sparta_bin"
done < "$case_list"
