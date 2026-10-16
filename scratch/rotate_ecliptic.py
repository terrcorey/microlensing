"""One-off (session 31): search.py outputs from before the ecliptic switch, rotated equatorial -> ecliptic.

    .venv/bin/python scratch/rotate_ecliptic.py      (login-safe: np.load / np.savez / renames only)

Exact: (piE_N, piE_E) and the Earth offsets (dsN, dsE) rotate by the same angle, so piE . ds and piE x ds --
all the model sees -- are unchanged and the saved chi2 / log_prob stay valid. Per input/*.toml run with a
results/<short>/refined.npz: rotates refined*.npz and fspl.npz in place, moves the old top-level MCMC +
diagnose outputs (sampled flat in the piE vector) into results/<short>/pi_E/ and rotates those chains.
Archived subdirs (pre24, session22, ...) and scratch/injections stay equatorial. Rotated files get
frame="ecliptic" and are skipped on a rerun. Future runs compute in the ecliptic frame from the start.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from event import load_event
from lc_models import sky_basis
from search import Binary

PIE = slice(Binary._fields.index("piE_N"), Binary._fields.index("piE_E") + 1)  # same slots in FSPL
MOVED = ("summary.txt", "fit_lc.png", "fit_lc_noplanet.png", "fspl_parallax_map.png")

for toml in sorted(Path("input").glob("*.toml")):
    ev = load_event(str(toml))
    out = Path("results") / ev.short_name
    if not (out / "refined.npz").exists():
        continue
    (n, e), (n0, e0) = sky_basis(ev.coords), sky_basis(ev.coords, pole=(0, 0, 1))
    rot = np.array([[n @ n0, n @ e0], [e @ n0, e @ e0]])  # ecliptic (N, E) = rot @ equatorial (N, E)
    sub = out / "pi_E"
    sub.mkdir(exist_ok=True)
    for f in [*out.glob("mcmc_mode*"), *(out / m for m in MOVED if (out / m).exists())]:
        f.rename(sub / f.name)
    if (sub / "summary.txt").exists() and "session 31" not in (s := (sub / "summary.txt").read_text()):
        (sub / "summary.txt").write_text("# NOTE (session 31): written in equatorial N/E; refined.npz + chains since rotated "
                                         "to ecliptic -- rerun --param pi_E --stage=diagnose to refresh\n" + s)
    for f in [*out.glob("refined*.npz"), out / "fspl.npz", *sub.glob("mcmc_mode*_chain.npz")]:
        if not f.exists() or "frame" in (d := dict(np.load(f))):
            continue
        for k in ("theta", "chain", "start"):
            if k in d:
                d[k] = d[k].copy()
                d[k][..., PIE] = d[k][..., PIE] @ rot.T
        with open(f.with_suffix(".tmp"), "wb") as fh:
            np.savez(fh, **d, frame="ecliptic")
        f.with_suffix(".tmp").replace(f)
        print(f"rotated {f} by {np.degrees(np.arctan2(rot[0, 1], rot[0, 0])):+.3f} deg")
