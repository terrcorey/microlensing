#!/bin/bash
# The four scratch/cross_check_mulensmodel*.py validations against MulensModel
# (needs `pip install MulensModel` in .venv). _fit reads pspl_235's
# fit_summary.dat, _fs reads cassan's chain -- submit_all.sh chains on both.
# Single-process, so 1 core. Submit from slurm/ (see 2l1s_mcmc.sbatch's header).

#SBATCH --job-name=cross-checks
#SBATCH --partition=small-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00
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
for f in scratch/cross_check_mulensmodel*.py; do python3 "$f"; done
