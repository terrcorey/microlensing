#!/bin/bash
# scratch/fit_2l1s.py's Nelder-Mead multi-start (joint then MOA-only, finite
# source), then scratch/compare_2l1s_fits.py's overlay plot.
# fit_2l1s.py's ProcessPoolExecutor() sizes itself from os.process_cpu_count(),
# which reads the whole node, not SLURM's allocation -- PYTHON_CPU_COUNT (Python
# 3.13+) caps it without touching the script.
# Submit from slurm/ (see 2l1s_mcmc.sbatch's header).

#SBATCH --job-name=2l1s-fit
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
python3 scratch/fit_2l1s.py
python3 scratch/compare_2l1s_fits.py
