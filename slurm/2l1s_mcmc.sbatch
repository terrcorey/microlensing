#!/bin/bash
# Submits one of scratch/mcmc_fit_2l1s.py's MCMCs, joint then MOA-only, to a CPU
# partition. Optional arg picks the likelihood: none = Student-t (run_mcmc, the
# script's own __main__), `chi2` = run_mcmc_chi2, `huber` = run_mcmc_huber.
#   sbatch 2l1s_mcmc.sbatch [chi2|huber]
#
# CPU, not GPU -- see CLAUDE.md / conversation: the Pool() parallelism is 48
# independent walker processes, each doing small (~hundreds-of-points) per-step
# work; a shared GPU across 48 processes without MPS would serialize rather
# than parallelize.
#
# small-short caps at 24 cores, below run_mcmc()'s nwalkers=48 -- no 1:1
# walker-to-core mapping possible, but that's fine: Pool(processes=24) (set via
# SLURM_CPUS_PER_TASK, see mcmc_fit_2l1s.py's run_mcmc()) just chunks the 48
# walkers into 2 batches of 24 per step instead of 1, still using every core
# with no oversubscription -- roughly half the per-step throughput of a
# hypothetical 48-core job, not broken.
#
# Lives in microlensing/slurm/; --output/--error resolve into slurm/output/ and slurm/errors/, so
# submit from slurm/ (submit_all.sh does this, and creates both dirs).

#SBATCH --job-name=2l1s-mcmc
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
fn="run_mcmc${1:+_$1}"
python3 -c "import sys; sys.path.insert(0, 'scratch'); import mcmc_fit_2l1s as m; m.$fn(); m.$fn(use_ogle=False)"
