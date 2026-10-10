"""One-off: are O-05-BLG169's fits consistent with Batista et al. 2015's Keck proper motion?
Light-curve side: t* = rho tE (the caustic crossing fixes it); Keck side: t* = theta* / mu_geo with
mu_geo = 7.0 +/- 0.2 mas/yr (their geocentric value, (N, E) = (4.79, 5.16)) and theta* = 0.44 +/- 0.04 uas
(Gould et al. 2006's colour-based source radius, quoted by Batista). Also piE vs the Keck-implied
piE = (pi_rel / theta_E) mu_geo/|mu_geo| = (0.11 / 0.838) (4.79, 5.16) / 7.04. Reads saved chains only
(login-node safe): `.venv/bin/python scratch/keck_mu_check.py`."""
from pathlib import Path

import numpy as np

MU, THETA_STAR = (7.0, 0.2), (0.44e-3, 0.04e-3)  # mas/yr, mas
TS = THETA_STAR[0] / MU[0] * 365.25
TS_ERR = TS * np.hypot(THETA_STAR[1] / THETA_STAR[0], MU[1] / MU[0])
PIE = 0.11 / 0.838 * np.array([4.79, 5.16]) / np.hypot(4.79, 5.16)  # (N, E)
RUNS = ["O-05-BLG169/pre24", "O-05-BLG169-fineq/session22", "O-05-BLG169-fineq/session23a",
        "O-05-BLG169-noAuckland/pre24", "O-05-BLG169-noFTN/pre24", "O-05-BLG169-noMDM/pre24"]

print(f"Keck: t* = {TS:.4f} +/- {TS_ERR:.4f} d, piE (N, E) = ({PIE[0]:.3f}, {PIE[1]:.3f})\n")
print(f"{'run / mode':40s} {'best chi2':>9s} {'t* 16/50/84 (d)':>24s} {'z(t*)':>6s} {'best t*':>7s} {'in 2sig':>7s} {'piE_N':>6s} {'piE_E':>6s} {'d(piE)':>6s}")
for run in RUNS:
    files = sorted(Path("results", run).glob("mcmc_mode*_chain.npz"), key=lambda f: int(f.stem.split("_")[1][4:]))
    for f in files:
        d = np.load(f)
        names, chain, lp = list(d["labels"]), d["chain"], d["log_prob"]
        s = chain[len(chain) // 4:].reshape(-1, chain.shape[2])
        col = lambda n: s[:, names.index(n)]
        ts = col("rho") * col("tE")
        lo, med, hi = np.percentile(ts, [16, 50, 84])
        z = (med - TS) / np.hypot(hi - med if med < TS else med - lo, TS_ERR)
        b = chain[np.unravel_index(np.argmax(lp), lp.shape)]
        tsb, inband = b[names.index("rho")] * b[names.index("tE")], np.mean(np.abs(ts - TS) < 2 * TS_ERR)
        pie = np.c_[col("piE_N"), col("piE_E")]
        dp = PIE - pie.mean(0)  # Mahalanobis distance of the Keck piE from the chain's piE cloud
        dist = np.sqrt(dp @ np.linalg.solve(np.cov(pie.T), dp))
        print(f"{run + ' ' + f.stem.split('_')[1]:40s} {-2 * lp.max():9.1f} {lo:7.4f} {med:7.4f} {hi:7.4f} "
              f"{z:+6.1f} {tsb:7.4f} {inband:7.2f} {np.median(pie[:, 0]):+6.2f} {np.median(pie[:, 1]):+6.2f} {dist:6.1f}")
    print()
