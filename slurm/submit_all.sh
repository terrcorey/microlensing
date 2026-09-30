#!/bin/bash
# Regenerates every output: both PSPL pipelines at once, then every 2L1S job
# once pspl_235 succeeds (they read its fit_summary.dat).
#   bash slurm/submit_all.sh    (from anywhere)
# Skipped: scratch/fit_2l1s_moa19008.py (data/MOA-2019-BLG-008L.dat not downloaded).
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p output errors

sbatch pspl_086.sbatch
pspl_235=$(sbatch --parsable pspl_235.sbatch)
after="--dependency=afterok:$pspl_235"

sbatch "$after" 2l1s_fit.sbatch
cassan=$(sbatch --parsable "$after" cassan.sbatch)
sbatch "$after,afterok:$cassan" cross_checks.sbatch
sbatch "$after,afterok:$cassan" compare.sbatch
sbatch "$after" 2l1s_mcmc.sbatch
sbatch "$after" --job-name=2l1s-mcmc-chi2 2l1s_mcmc.sbatch chi2
sbatch "$after" --job-name=2l1s-mcmc-huber 2l1s_mcmc.sbatch huber
