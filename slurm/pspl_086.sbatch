#!/bin/bash
# O-05-BLG086's full PSPL pipeline (mcmc_fit.py): results/O-05-BLG086/pspl/,
# data/processed/O-05-BLG086_*. Serial emcee (no Pool), so 1 core.
# Submit from slurm/ (see 2l1s_mcmc.sbatch's header).

#SBATCH --job-name=pspl-086
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
python3 mcmc_fit.py --stage=raw
python3 mcmc_fit.py
