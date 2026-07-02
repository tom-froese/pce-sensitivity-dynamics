#!/usr/bin/env python3
"""Within-trial neural-complexity intermediates for the decoupling Bayes factors.

Ports the two complexity measures from the master-loop's within-trial complexity
analysis (``atlas/within_trial_complexity.py``) that the neural<->phenomenal
decoupling test (``computeDecouplingBF.py``) consumes:

  * **multichannel Lempel-Ziv complexity slope** (``mc_slope``), per trial — the
    within-trial state-space *complexity-collapse rate*. Computed exactly as in
    the master ``trial_complexity()``: for each 4 s window (2 s step over the
    [2, 58] s segment), take each channel's analytic-signal (Hilbert) envelope,
    decimate by 5 (250 -> 50 Hz), binarize each channel at its within-window
    median, flatten the 64 x T binary matrix in C order, and take the
    normalized LZ76 complexity; ``mc_slope`` is the least-squares slope of that
    per-window series against window index.

  * **global spectral entropy** (``specent_global``), per window — computed as in
    the master ``window_measures()``: normalized Welch spectral entropy of the
    channel-mean signal (``nperseg = 2*sf`` samples).

Reads the cleaned continuous EEG (``data/preprocessed/EEG/pceXX/
pceXX_PY_task-raw.fif``, 250 Hz, produced by ``preprocessEEGForExponent.py`` or
downloaded from Zenodo) via the shared ``_eeg_io`` loaders. Writes two long-format
CSVs to ``results/``:

  * ``within_trial_complexity_pertrial.csv``  (one row per subject x trial;
    columns ``subj, dyad, participant, trial, mc_slope``)
  * ``within_trial_windows.csv``              (one row per subject x trial x
    window; columns ``subj, dyad, participant, trial, trial_time,
    specent_global``)

These are validated against the master's committed intermediates to be
bit-identical (mc_slope) / within 1e-5 (spectral entropy).

Usage
-----
    python code/analysis/computeWithinTrialComplexity.py            # all dyads
    python code/analysis/computeWithinTrialComplexity.py --dyads 1 2 3   # smoke
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import hilbert

_HERE = Path(__file__).resolve().parent
_CODE = _HERE.parent
for _p in (str(_CODE), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import _config as C  # noqa: E402
import _eeg_io as A  # noqa: E402

# Window grid — identical to the master (atlas/_common).
WIN_S, STEP_S = 4.0, 2.0
SEG = (2.0, 58.0)
DECIMATE = 5                    # envelope decimation for multichannel LZc (250 -> 50 Hz)


def _lz(bits) -> float:
    import antropy
    return float(antropy.lziv_complexity(bits, normalize=True))


def _mc_series(X: np.ndarray, sf: float) -> np.ndarray:
    """Per-window multichannel LZc series for one (64, n) trial segment."""
    L = int(round(WIN_S * sf))
    step = int(round(STEP_S * sf))
    env = np.abs(hilbert(X, axis=1))
    out = []
    for s in range(0, X.shape[1] - L + 1, step):
        w = env[:, s:s + L][:, ::DECIMATE]
        bmat = (w > np.median(w, axis=1, keepdims=True)).astype(np.int8)
        out.append(_lz(bmat.flatten(order="C")))
    return np.array(out)


def _specent_window(Xw: np.ndarray, sf: float) -> float:
    """Normalized Welch spectral entropy of the channel-mean signal in a window."""
    import antropy
    gm = Xw.mean(0)
    return float(antropy.spectral_entropy(
        gm, sf, method="welch", normalize=True,
        nperseg=min(len(gm), int(2 * sf))))


def compute(ds: A.Dataset, dyads=None):
    if dyads is None:
        dyads = sorted(int(p.name[3:5]) for p in C.PREPROC_EEG.glob("pce*")
                       if p.is_dir() and ds.raw_path(int(p.name[3:5]), 1).exists())
    pertrial, windows = [], []
    t0 = time.time()
    for d in dyads:
        try:
            prov = A.load_provenance(ds, d)
            raws = A.load_dyad_raws(ds, d)
        except Exception as exc:
            print(f"  [dyad {d:02d}] skipped: {type(exc).__name__}: {exc}", flush=True)
            continue
        loaded = {p: prov["trials_loaded"][f"P{p}"] for p in C.PARTICIPANTS}
        onsets = {p: A.segment_onsets(raws[p], loaded[p]) for p in C.PARTICIPANTS}
        sf = float(raws[C.PARTICIPANTS[0]].info["sfreq"])
        L = int(round(WIN_S * sf))
        step = int(round(STEP_S * sf))
        for p in C.PARTICIPANTS:
            subj = f"{d:02d}_{p}"
            for trial in loaded[p]:
                s0 = int(round((onsets[p][trial] + SEG[0]) * sf))
                s1 = int(round((onsets[p][trial] + SEG[1]) * sf))
                X = raws[p].get_data(start=s0, stop=s1)
                mc = _mc_series(X, sf)
                n_win = len(mc)
                pertrial.append(dict(
                    subj=subj, dyad=d, participant=p, trial=int(trial),
                    mc_slope=float(np.polyfit(np.arange(n_win), mc, 1)[0])))
                for i, st in enumerate(range(0, X.shape[1] - L + 1, step)):
                    ci = SEG[0] + i * STEP_S + WIN_S / 2
                    windows.append(dict(
                        subj=subj, dyad=d, participant=p, trial=int(trial),
                        trial_time=float(ci),
                        specent_global=_specent_window(X[:, st:st + L], sf)))
        print(f"  [complexity] dyad {d:02d} ({time.time()-t0:.0f}s)", flush=True)
    return pd.DataFrame(pertrial), pd.DataFrame(windows)


def run(dyads=None):
    ds = A.default_dataset()
    pt, w = compute(ds, dyads)
    C.RESULTS.mkdir(parents=True, exist_ok=True)
    pt_path = C.RESULTS / "within_trial_complexity_pertrial.csv"
    w_path = C.RESULTS / "within_trial_windows.csv"
    pt.to_csv(pt_path, index=False)
    w.to_csv(w_path, index=False)
    print(f"  {len(pt)} trials, {len(w)} windows, {pt.subj.nunique()} participants")
    print(f"  saved -> {pt_path.relative_to(C.REPO_ROOT)}, {w_path.relative_to(C.REPO_ROOT)}")
    return pt, w


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dyads", type=int, nargs="+", default=None,
                    help="restrict to these dyad numbers (default: all available)")
    args = ap.parse_args()
    print("\n=== within-trial complexity (mc_slope + spectral entropy) ===", flush=True)
    run(dyads=args.dyads)


if __name__ == "__main__":
    main()
