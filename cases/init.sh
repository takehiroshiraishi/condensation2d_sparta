#!/usr/bin/env bash

set -euo pipefail

cases_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
study_root="$(cd "$cases_dir/.." && pwd)"
template_path="$cases_dir/_templates/parameters.json"

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 STUDY_NAME" >&2
  exit 1
fi

study_name="$1"
study_dir="$cases_dir/$study_name"
parameters_path="$study_dir/parameters.json"

if [[ ! "$study_name" =~ ^[A-Za-z0-9._-]+$ ]]; then
  echo "Invalid study name: $study_name" >&2
  echo "Allowed characters: letters, digits, dot, underscore, hyphen." >&2
  exit 1
fi

if [[ ! -f "$template_path" ]]; then
  echo "Template not found: $template_path" >&2
  exit 1
fi

if [[ -e "$study_dir" ]]; then
  echo "Study directory already exists: $study_dir" >&2
  exit 1
fi

mkdir -p "$study_dir"
sed "s/replace_with_study_name/$study_name/g" "$template_path" > "$parameters_path"

echo "Initialized $study_dir"
echo "Edit $parameters_path, then run:"
echo "  (cd \"$study_root\" && python3 scripts/generate_cases.py --config cases/$study_name/parameters.json --force)"
