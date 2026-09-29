#!/bin/bash
# scratch/cassan_caustic.py: both Cassan fits (from Bond's and session 11's
# Student-t (s, q)), then the finite-source MCMC around the better one
# (~0.8 s/chi2). Overwrites scratch/2l1s/cassan/'s chain with the 7-column one.
# Pool() sizing: see 2l1s_fit.sbatch's PYTHON_CPU_COUNT note.
# Submit from slurm/ (see 2l1s_mcmc.sbatch's header).

#SBATCH --job-name=cassan
#SBATCH --partition=small-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=2G
#SBATCH --time=12:00:00
#SBATCH --output=output/slurm-%x-%j.out
#SBATCH --error=errors/slurm-%x-%j.err

set -euo pipefail
cd "$HOME/work/microlensing"

source .venv/bin/activate
# lc_models.py probes torch.cuda on import; on these CPU nodes that emits a
# driver-too-old UserWarning per process (and per Pool worker) -- noise in .err.
export PYTHONWARNINGS="ignore::UserWarning:torch.cuda"
# stdout to a file is block-buffered -- without this, prints only reach .out in chunks.
export PYTHONUNBUFFERED=1
export PYTHON_CPU_COUNT=$SLURM_CPUS_PER_TASK
python3 scratch/cassan_caustic.py
