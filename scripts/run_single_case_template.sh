#!/bin/bash
#SBATCH --nodes 1
#SBATCH --job-name condensation2d
#SBATCH --ntasks-per-node 32
#SBATCH --cpus-per-task 1
#SBATCH --time 7-00:00:00

set -euo pipefail

case_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
study_root="$(cd "$case_dir/../../.." && pwd)"

module purge
module load gcc
module load openmpi

if [[ -n "${SPARTA_BIN:-}" ]]; then
  sparta_bin="$SPARTA_BIN"
elif [[ -x "$case_dir/spa_mpi" ]]; then
  sparta_bin="$case_dir/spa_mpi"
elif [[ -x "$study_root/spa_mpi" ]]; then
  sparta_bin="$study_root/spa_mpi"
else
  echo "Could not find spa_mpi. Set SPARTA_BIN, place spa_mpi in the case directory, or place one shared spa_mpi at $study_root/spa_mpi." >&2
  exit 1
fi

cd "$case_dir"
srun "$sparta_bin" < in.condensation > log.txt
