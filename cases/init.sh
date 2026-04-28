#!/usr/bin/env bash

set -euo pipefail

cases_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
study_root="$(cd "$cases_dir/.." && pwd)"
parameters_template_path="$cases_dir/_templates/parameters.json"
plot_template_path="$cases_dir/_templates/plot_y_axis.gp"

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 STUDY_NAME" >&2
  exit 1
fi

study_name="$1"
study_dir="$cases_dir/$study_name"
parameters_path="$study_dir/parameters.json"
generate_script_path="$study_dir/generate.sh"
plot_script_path="$study_dir/plot_y_axis.gp"

if [[ ! "$study_name" =~ ^[A-Za-z0-9._-]+$ ]]; then
  echo "Invalid study name: $study_name" >&2
  echo "Allowed characters: letters, digits, dot, underscore, hyphen." >&2
  exit 1
fi

if [[ ! -f "$parameters_template_path" ]]; then
  echo "Template not found: $parameters_template_path" >&2
  exit 1
fi

if [[ ! -f "$plot_template_path" ]]; then
  echo "Template not found: $plot_template_path" >&2
  exit 1
fi

if [[ -e "$study_dir" ]]; then
  echo "Study directory already exists: $study_dir" >&2
  exit 1
fi

mkdir -p "$study_dir"
cp "$parameters_template_path" "$parameters_path"
python3 - "$parameters_path" "$study_name" <<'EOF'
import json
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
study_name = sys.argv[2]

with path.open("r", encoding="utf-8") as handle:
    data = json.load(handle)

data["study_name"] = study_name

with path.open("w", encoding="utf-8") as handle:
    json.dump(data, handle, indent=2, sort_keys=False)
    handle.write("\n")
EOF
cp "$plot_template_path" "$plot_script_path"
cat > "$generate_script_path" <<EOF
#!/usr/bin/env bash

set -euo pipefail

study_dir="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"
study_root="\$(cd "\$study_dir/../.." && pwd)"

cd "\$study_root"
python3 scripts/generate_cases.py --config "cases/$study_name/parameters.json" --force
EOF
chmod +x "$generate_script_path"

echo "Initialized $study_dir"
echo "Edit $parameters_path, then run:"
echo "  $generate_script_path"
