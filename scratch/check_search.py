"""Self-check for search.py's session-22 robustness fixes (~1 min; srun, not the login node):

    srun --partition=small-short --cpus-per-task=2 --mem=4G --time=00:10:00 .venv/bin/python scratch/check_search.py

1. chi2 guards: NaN / tiny-rho / out-of-box inputs return inf before reaching VBBL (which
   hangs on NaN and segfaulted in that regime).
2. pool(): a worker dying hard raises BrokenProcessPool instead of hanging.
3. run_mcmc(): checkpoints, and a second call resumes to nsteps instead of restarting.
"""
import os
import sys
import tempfile
from concurrent.futures.process import BrokenProcessPool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np

import search
from event import load_event

if __name__ == "__main__":  # spawned pool workers re-import this file
    event = load_event("input/O-03-BLG235.toml")
    good = search.Binary(2848.0, 0.1, 60.0, 1e-3, 0.0, 0.0, 1.12, 4e-3, 3.9)
    assert np.isfinite(search.chi2_binary(good, event))
    for bad in (good._replace(u0=np.nan), good._replace(rho=1e-9), good._replace(q=1e-7),
                good._replace(s=20.0), good._replace(piE_N=6.0), good._replace(tE=np.inf)):
        assert search.chi2_binary(bad, event) == np.inf, bad
    assert np.isfinite(search.chi2_binary(good._replace(rho=0.0), event))  # point source still allowed
    assert search.chi2_fspl((2848.0, 0.1, 60.0, 1e-15, 0, 0), event) == np.inf
    print("guards ok")

    with search.pool(event) as ex:
        try:
            list(ex.map(os._exit, [1]))
            raise AssertionError("dead worker did not raise")
        except BrokenProcessPool:
            print("dead worker raises ok")

    with tempfile.TemporaryDirectory() as d, search.pool(event) as ex:
        path = Path(d) / "mcmc_x_chain.npz"
        search.run_mcmc(ex, good, path, "x", nsteps=6, every=3)
        assert np.load(path)["chain"].shape == (6, 32, 9)
        search.run_mcmc(ex, good, path, "x", nsteps=10, every=3)
        c = np.load(path)
        assert c["chain"].shape == (10, 32, 9) and c["log_prob"].shape == (10, 32)
        search.save_mcmc(path, Path(d), "x")
        print("checkpoint + resume ok")
