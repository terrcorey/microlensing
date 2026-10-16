#!/bin/bash
# scratch/sensitivity.py: an event's planet sensitivity (Delta chi2) map, from its results/<short>/ run.
#   sbatch sensitivity.sbatch ../input/<short>.toml   (submit from slurm/)

#SBATCH --job-name=sensitivity
#SBATCH --partition=large-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=1G
#SBATCH --time=12:00:00
#SBATCH --output=output/slurm-%x-%j.out
#SBATCH --error=errors/slurm-%x-%j.err

set -euo pipefail
config=$(realpath "$1")
cd "$HOME/work/microlensing"

source .venv/bin/activate
export PYTHONWARNINGS="ignore::UserWarning:torch.cuda"
export PYTHONUNBUFFERED=1
export PYTHON_CPU_COUNT=$SLURM_CPUS_PER_TASK  # ProcessPoolExecutor's default worker count
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 scratch/sensitivity.py "$config"
