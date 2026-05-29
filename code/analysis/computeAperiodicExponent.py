#!/usr/bin/env python3
"""Per-(participant × within-trial bin) FOOOF aperiodic EXPONENT.

For each participant and each 4 s within-trial window (2 s step over [2, 58] s),
this computes the cohort-averaged Welch PSD across that participant's trials in
that window, then runs a fixed-mode FOOOF aperiodic fit on **2-40 Hz** (the
canonical exponent, ``E``) and on the **EMG-reduced 2-20 Hz** band (``E_lo``,
the high-frequency-clean confound check).

Outputs (in ``results/``):
  ``aperiodic_exponent_per_participant.csv``              — 2-40 Hz exponent
                                                            (subjects × bin centers)
  ``aperiodic_exponent_per_participant_band-2-20.csv``    — 2-20 Hz exponent
  ``aperiodic_exponent_per_participant.json``             — provenance sidecar

Reads the cleaned noICA .fif files produced by
:mod:`preprocessEEGForExponent`.

Usage
-----
    python code/analysis/computeAperiodicExponent.py             # all dyads
    python code/analysis/computeAperiodicExponent.py --dyads 1 2 3   # smoke
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
_CODE = _HERE.parent
for _p in (str(_CODE), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import _config as C       # noqa: E402
import _eeg_io as A       # noqa: E402
from _exponent_common import WIN_S, STEP_S, SEG, OUT  # noqa: E402


# FOOOF settings — pinned here so the JSON sidecar always agrees with the run.
# Match the master-loop's atlas/per_subj_exponent.FOOOF_SETTINGS exactly.
FOOOF_SETTINGS = {
    "wideband":    {"freq_range": [2, 40], "aperiodic_mode": "fixed", "max_n_peaks": 6},
    "emg_reduced": {"freq_range": [2, 20], "aperiodic_mode": "fixed", "max_n_peaks": 4},
}
WELCH_NPERSEG_S = 2.0   # Welch segment length in seconds; nperseg = int(2 * sf)


def compute(ds: A.Dataset, dyads=None):
    """Per-(participant × bin) FOOOF aperiodic exponent.

    Returns
    -------
    subjects : list[str]                — sorted ``"<dyad:02d>_<p>"`` IDs
    centers  : ndarray, shape (n_bin,)  — bin-centre times in seconds
    E        : ndarray, (n_subj, n_bin) — 2-40 Hz aperiodic exponent
    E_lo     : ndarray, (n_subj, n_bin) — 2-20 Hz exponent (EMG-reduced)
    """
    from scipy.signal import welch
    from fooof import FOOOF

    if dyads is None:
        dyads = sorted(int(p.name[3:5]) for p in C.PREPROC_EEG.glob("pce*")
                       if p.is_dir() and ds.raw_path(int(p.name[3:5]), 1).exists())

    L = int(round(WIN_S * C.SFREQ_TARGET))
    step = int(round(STEP_S * C.SFREQ_TARGET))
    n_centers = len(range(0, int(round((SEG[1] - SEG[0]) * C.SFREQ_TARGET)) - L + 1, step))
    centers = np.array([SEG[0] + i * STEP_S + WIN_S / 2 for i in range(n_centers)])

    acc = {}                              # subj -> [psum (n_centers, nf), cnt (n_centers,)]
    freqs = None
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
        for p in C.PARTICIPANTS:
            subj = f"{d:02d}_{p}"
            for trial in loaded[p]:
                s0 = int(round((onsets[p][trial] + SEG[0]) * sf))
                s1 = int(round((onsets[p][trial] + SEG[1]) * sf))
                X = raws[p].get_data(start=s0, stop=s1)
                for i, st in enumerate(range(0, X.shape[1] - L + 1, step)):
                    if i >= n_centers:
                        break
                    fr, P = welch(X[:, st:st + L], fs=sf,
                                  nperseg=int(WELCH_NPERSEG_S * sf), axis=1)
                    if freqs is None:
                        freqs = fr
                    if subj not in acc:
                        acc[subj] = [np.zeros((n_centers, len(fr))), np.zeros(n_centers)]
                    acc[subj][0][i] += P.mean(0)
                    acc[subj][1][i] += 1
        print(f"  [{ds.name}] dyad {d:02d} ({time.time()-t0:.0f}s)", flush=True)

    fmask = (freqs >= FOOOF_SETTINGS["wideband"]["freq_range"][0]) & \
            (freqs <= FOOOF_SETTINGS["wideband"]["freq_range"][1])
    fsel = freqs[fmask]
    fmask_lo = (freqs >= FOOOF_SETTINGS["emg_reduced"]["freq_range"][0]) & \
               (freqs <= FOOOF_SETTINGS["emg_reduced"]["freq_range"][1])
    fsel_lo = freqs[fmask_lo]

    subjects = sorted(acc)
    E = np.full((len(subjects), n_centers), np.nan)
    E_lo = np.full((len(subjects), n_centers), np.nan)
    for si, subj in enumerate(subjects):
        psum, cnt = acc[subj]
        for i in range(n_centers):
            if cnt[i] < 3:
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pm = psum[i] / cnt[i]
                fm = FOOOF(aperiodic_mode=FOOOF_SETTINGS["wideband"]["aperiodic_mode"],
                           max_n_peaks=FOOOF_SETTINGS["wideband"]["max_n_peaks"],
                           verbose=False)
                fm.fit(fsel, pm[fmask], FOOOF_SETTINGS["wideband"]["freq_range"])
                E[si, i] = float(fm.get_params("aperiodic_params", "exponent"))
                fl = FOOOF(aperiodic_mode=FOOOF_SETTINGS["emg_reduced"]["aperiodic_mode"],
                           max_n_peaks=FOOOF_SETTINGS["emg_reduced"]["max_n_peaks"],
                           verbose=False)
                fl.fit(fsel_lo, pm[fmask_lo], FOOOF_SETTINGS["emg_reduced"]["freq_range"])
                E_lo[si, i] = float(fl.get_params("aperiodic_params", "exponent"))
    return subjects, centers, E, E_lo


def save(subjects, centers, E, E_lo, out=OUT):
    """Write per-subj CSVs + provenance sidecar.

    File names omit a pipeline suffix because the spoke targets only the
    noICA pipeline.
    """
    out.mkdir(parents=True, exist_ok=True)
    csv_wb = out / "aperiodic_exponent_per_participant.csv"
    csv_lo = out / "aperiodic_exponent_per_participant_band-2-20.csv"
    json_sidecar = out / "aperiodic_exponent_per_participant.json"

    pd.DataFrame(E,    index=subjects, columns=centers).to_csv(csv_wb)
    pd.DataFrame(E_lo, index=subjects, columns=centers).to_csv(csv_lo)

    sidecar = {
        "producer": "code/analysis/computeAperiodicExponent.py",
        "produced_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "pipeline": "noICA (task-raw)",
        "input_eeg": "data/preprocessed/EEG/pce*/pce*_P{1,2}_task-raw.fif",
        "n_participants": len(subjects),
        "bin_centers_s": [float(c) for c in centers],
        "window_s": WIN_S,
        "step_s": STEP_S,
        "segment_s": list(SEG),
        "welch": {"nperseg_s": WELCH_NPERSEG_S, "method": "scipy.signal.welch"},
        "fooof": FOOOF_SETTINGS,
        "outputs": {
            "wideband_2-40Hz": csv_wb.name,
            "emg_reduced_2-20Hz": csv_lo.name,
        },
        "notes": (
            "Per-(participant × bin) FOOOF aperiodic exponent. PSDs are cohort-averaged "
            "across that participant's trials in that bin (channel-mean Welch). FOOOF "
            "fit in 'fixed' aperiodic mode on the indicated freq range."
        ),
    }
    json_sidecar.write_text(json.dumps(sidecar, indent=2))
    print(f"  saved -> {csv_wb.name}, {csv_lo.name}, {json_sidecar.name}", flush=True)


def run(dyads=None):
    ds = A.default_dataset()
    subjects, centers, E, E_lo = compute(ds, dyads)
    save(subjects, centers, E, E_lo)
    return subjects, centers, E, E_lo


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dyads", type=int, nargs="+", default=None,
                    help="restrict to these dyad numbers (default: all available)")
    args = ap.parse_args()
    print("\n=== aperiodic exponent (per participant × bin) ===", flush=True)
    run(dyads=args.dyads)


if __name__ == "__main__":
    main()
