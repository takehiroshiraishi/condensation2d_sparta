#!/usr/bin/env bash
# Fill in the SBATCH placeholders before submitting.

#SBATCH --job-name=condensation
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16

set -euo pipefail

study_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
case_list="${CASELIST:-$study_root/cases/case_list.txt}"
sparta_bin="${SPARTA_BIN:-$study_root/../../src/spa_mpi}"

case_relpath="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$case_list")"
if [[ -z "$case_relpath" ]]; then
  echo "No case found for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi

export SPARTA_LAUNCH="${SPARTA_LAUNCH:-srun}"
"$study_root/scripts/run_case.sh" "$study_root/$case_relpath" "$sparta_bin"
