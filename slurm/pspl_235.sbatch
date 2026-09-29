#!/bin/bash
# O-03-BLG235's PSPL pipeline: preprocess_binary_data.py -> mcmc_fit_binary.py.
# Every 2L1S job depends on this one -- scratch/fit_2l1s.py reads its
# data/processed/O-03-BLG235_fit_summary.dat at import time. Serial emcee, 1 core.
# Submit from slurm/ (see 2l1s_mcmc.sbatch's header).

#SBATCH --job-name=pspl-235
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
python3 preprocess_binary_data.py
python3 mcmc_fit_binary.py --stage=raw
python3 mcmc_fit_binary.py
