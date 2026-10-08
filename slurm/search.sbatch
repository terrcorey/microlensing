#!/bin/bash
# search.py: the config-driven FSPL vs 2L1S pipeline for one event.
#   sbatch search.sbatch ../input/<short>.toml   (submit from slurm/; path relative to the project root works too)
#   optional $2 = --stage (e.g. diagnose); override resources on the command line for light stages:
#   sbatch --cpus-per-task=1 --time=00:30:00 --dependency=afterok:<search job> search.sbatch ../input/<short>.toml diagnose
# Pool() sizing: see 2l1s_fit.sbatch's PYTHON_CPU_COUNT note.

#SBATCH --job-name=search
#SBATCH --partition=small-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=2G
#SBATCH --time=2-00:00:00
#SBATCH --output=output/slurm-%x-%j.out
#SBATCH --error=errors/slurm-%x-%j.err

set -euo pipefail
config=$(realpath "${1:-../input/O-03-BLG235.toml}")
cd "$HOME/work/microlensing"

source .venv/bin/activate
# lc_models.py probes torch.cuda on import; on these CPU nodes that emits a
# driver-too-old UserWarning per process (and per Pool worker) -- noise in .err.
export PYTHONWARNINGS="ignore::UserWarning:torch.cuda"
# stdout to a file is block-buffered -- without this, prints only reach .out in chunks.
export PYTHONUNBUFFERED=1
export PYTHON_CPU_COUNT=$SLURM_CPUS_PER_TASK
# one BLAS thread per worker: 24 workers x a BLAS pool each would oversubscribe the 24 cores
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 search.py --config "$config" ${2:+--stage=$2}
