#!/bin/bash
# scratch/compare_pspl_2l1s.py: PSPL fit + Nelder-Mead polish of the Cassan chain's
# best sample (~0.8 s/chi2, serial), then BIC. Reads cassan.sbatch's chain.
# Submit from slurm/ (see 2l1s_mcmc.sbatch's header).

#SBATCH --job-name=compare
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
# see cassan.sbatch for these two
export PYTHONWARNINGS="ignore::UserWarning:torch.cuda"
export PYTHONUNBUFFERED=1
python3 scratch/compare_pspl_2l1s.py
