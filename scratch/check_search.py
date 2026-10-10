"""Self-check for search.py's session-22 robustness fixes (~1 min; srun, not the login node):

    srun --partition=small-short --cpus-per-task=2 --mem=4G --time=00:10:00 .venv/bin/python scratch/check_search.py

1. chi2 guards: NaN / out-of-box inputs return inf before reaching VBBL (which hangs on NaN and
   segfaulted at tiny rho); in 2L1S, 0 <= rho < RHO_MIN is the exact point source (session 23).
2. pool(): a worker dying hard raises BrokenProcessPool instead of hanging.
3. run_mcmc(): checkpoints, and a second call resumes to nsteps instead of restarting.
4. caustic_crossing() (session 23): straight down the axis of a resonant caustic crosses it,
   0.5 thetaE off it is a near miss beyond CLOSE.
5. free LD (session 24): load_event's per-band consistency, with_ld, the prior term, LD bounds,
   the MCMC coordinate round trip and log_prob with LD appended.
6. session 25: tE above Event.tE_max is inf; run_grid() scores a cell past CELL_TIMEOUT inf and
   carries on (pebble kills the worker) instead of waiting on it forever.
"""
import os
import sys
import tempfile
from concurrent.futures.process import BrokenProcessPool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np

import search
from event import ld_penalty, load_event, with_ld

if __name__ == "__main__":  # spawned pool workers re-import this file
    event = load_event("input/O-03-BLG235.toml")
    good = search.Binary(2848.0, 0.1, 60.0, 1e-3, 0.0, 0.0, 1.12, 4e-3, 3.9)
    assert np.isfinite(search.chi2_binary(good, event))
    for bad in (good._replace(u0=np.nan), good._replace(rho=-1e-3), good._replace(rho=0.2), good._replace(q=1e-7),
                good._replace(s=20.0), good._replace(piE_N=6.0), good._replace(tE=np.inf)):
        assert search.chi2_binary(bad, event) == np.inf, bad
    assert search.chi2_binary(good, event._replace(tE_max=50.0)) == np.inf  # tE 60 > ceiling
    point = search.chi2_binary(good._replace(rho=0.0), event)
    assert np.isfinite(point) and search.chi2_binary(good._replace(rho=1e-9), event) == point  # below floor = point
    assert search.chi2_fspl((2848.0, 0.1, 60.0, 1e-15, 0, 0), event) == np.inf
    assert np.allclose(search.from_mcmc(search.to_mcmc(good)), good)
    search._init(event)  # log_prob reads the module's EVENT
    x = search.to_mcmc(good)
    assert np.isfinite(search.log_prob(x, good.alpha)) and search.log_prob(x, good.alpha - 4.0) == -np.inf  # alpha bound
    print("guards + MCMC coordinates ok")
    trk = search.track(event, 2848.0)
    through = search.Binary(2848.0, 0.0, 60.0, 1e-3, 0.0, 0.0, 1.0, 0.1, 0.0)
    assert search.caustic_crossing(trk, through)[0]
    miss = search.caustic_crossing(trk, through._replace(u0=0.5, q=1e-4))
    assert not miss[0] and miss[1] == np.inf, miss
    print("caustic_crossing ok")

    with tempfile.TemporaryDirectory() as d:  # ld_sigma on OGLE only: MDM (also I) disagrees -> raise
        toml = "\n".join(l for l in Path("input/O-05-BLG169.toml").read_text().splitlines()
                         if not l.startswith("ld_sigma"))  # the config sets it since session 24
        Path(d, "bad.toml").write_text(toml.replace("ld = 0.53", "ld = 0.53\nld_sigma = 0.1", 1))
        try:
            load_event(str(Path(d, "bad.toml")))
            raise AssertionError("inconsistent band prior did not raise")
        except ValueError:
            pass
        Path(d, "ok.toml").write_text(toml.replace("ld = 0.53", "ld = 0.53\nld_sigma = 0.1"))
        assert load_event(str(Path(d, "ok.toml"))).ld_prior == (("I", 0.53, 0.1),)
    ev = event._replace(ld_prior=(("I", 0.5, 0.1),))  # O-03-BLG235: OGLE and MOA are both I
    assert [i.ld for i in with_ld(ev, [0.3]).instruments] == [0.3, 0.3] and with_ld(ev, []) is ev
    assert [i.ld for i in with_ld(ev._replace(ld_prior=(("V", 0.5, 0.1),)), [0.3]).instruments] == [
        i.ld for i in event.instruments]  # a band no instrument has: nothing changes
    assert np.isclose(ld_penalty(ev, [0.7]), 4.0) and ld_penalty(ev, []) == 0
    assert search.chi2_binary((*good, 1.5), ev) == np.inf and search.chi2_binary((*good, -0.1), ev) == np.inf
    th = np.array([*good, 0.7])
    assert np.allclose(search.from_mcmc(search.to_mcmc(th)), th)
    search._init(ev)
    assert np.isclose(search.log_prob(search.to_mcmc(th), good.alpha),
                      -0.5 * (search.chi2_binary(th, ev) + 4.0) - 2 * np.log(good.tE))
    search._init(event)
    print("free LD ok")

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

    with tempfile.TemporaryDirectory() as d:  # 2 x 2 point-source cells take ~1 min each: all time out
        search.CELL_TIMEOUT = 2  # read by run_grid in this process; the workers never see it
        tiny = event._replace(grid={**event.grid, "log_s": [0.0, 0.05, 0.05], "log_q": [-3.0, -2.75, 0.25]})
        g = search.run_grid(tiny, search.FSPL(2848.0, 0.1, 60.0, 2e-4, 0.0, 0.0), Path(d) / "grid.npz", "k",
                            np.ones(len(event.instruments)))
        assert np.all(g["chi2"] == np.inf), g["chi2"]
        print("cell timeout ok")
