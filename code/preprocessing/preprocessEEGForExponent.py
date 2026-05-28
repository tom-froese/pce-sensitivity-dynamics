#!/usr/bin/env python3
"""Raw PCE EEG (.mat) -> cleaned, bandpassed, dyad-aligned EEG dataset.

Spoke-side driver around :func:`hypyp_ext.preprocess_dyad`. For each dyad it
loads both participants' task trials, builds MNE Raws with the 64-ch actiCAP
montage, concatenates trials per participant, runs the moderate (no-ICA)
pipeline, and writes per-participant cleaned continuous data (``*-raw.fif``) and
cleaned epochs (``*-epo.fif``) plus a per-dyad provenance JSON.

The resulting ``data/preprocessed/EEG/pceXX/pceXX_PY_task-raw.fif`` files are
the input to :mod:`computeAperiodicExponent`. They are also archived on Zenodo
as part of ``preprocessed.zip`` for readers who do not want to download the
~15 GB raw EEG from OSF.

Usage
-----
    # pilot: one dyad
    python code/preprocessing/preprocessEEGForExponent.py --dyad 1

    # quick smoke test: few trials only
    python code/preprocessing/preprocessEEGForExponent.py --dyad 1 --limit-trials 3

    # full run, all dyads with EEG (~30 min)
    python code/preprocessing/preprocessEEGForExponent.py --all
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.io as sio
import mne

_HERE = Path(__file__).resolve().parent           # code/preprocessing
_CODE = _HERE.parent                              # code
for _p in (str(_CODE), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import _config as C  # noqa: E402
from hypyp_ext import preprocess_dyad  # noqa: E402

mne.set_log_level("ERROR")


# --------------------------------------------------------------------------- #
#  Raw .mat loading
# --------------------------------------------------------------------------- #
def load_trial_matrix(path: Path) -> np.ndarray:
    """Load a raw EEG .mat file as a ``(N_CHAN, n_times)`` float array."""
    md = sio.loadmat(str(path))
    arrays = [
        v for k, v in md.items()
        if not k.startswith("__") and isinstance(v, np.ndarray) and v.ndim == 2
    ]
    mat = None
    for v in arrays:                       # prefer the array with a channel-sized dim
        if C.N_CHAN in v.shape:
            mat = v
            break
    if mat is None:
        raise ValueError(f"No {C.N_CHAN}-channel matrix found in {path.name}")
    if mat.shape[0] != C.N_CHAN:
        mat = mat.T
    return np.asarray(mat, dtype=float)


def build_raw(mat: np.ndarray) -> mne.io.RawArray:
    """Wrap a ``(N_CHAN, n_times)`` matrix as an MNE Raw with montage + units."""
    info = mne.create_info(list(C.CH_NAMES), C.SFREQ_ORIG, ch_types="eeg")
    raw = mne.io.RawArray(mat * C.EEG_UNIT_SCALE, info, verbose=False)
    montage = mne.channels.make_standard_montage(C.MONTAGE)
    raw.set_montage(montage, on_missing="warn", verbose=False)
    return raw


def load_participant_task(folder: Path, dnum: int, p: int, limit_trials: int | None):
    """Concatenate a participant's task trials into one continuous Raw."""
    trials = C.TRIALS if limit_trials is None else C.TRIALS[:limit_trials]
    raws, loaded = [], []
    for t in trials:
        fp = folder / f"pce{dnum:02d}_P{p}_Trial{t}.mat"
        if not fp.exists():
            continue
        raws.append(build_raw(load_trial_matrix(fp)))
        loaded.append(t)
    if not raws:
        return None, []
    raw = mne.concatenate_raws(raws, verbose=False) if len(raws) > 1 else raws[0]
    return raw, loaded


def _jsonable(obj):
    """Recursively coerce numpy types so json.dump succeeds."""
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


# --------------------------------------------------------------------------- #
#  Per-dyad processing
# --------------------------------------------------------------------------- #
def process_dyad(dnum: int, limit_trials: int | None, run_ar: bool) -> bool:
    folder = next(iter(sorted(C.RAW_EEG.glob(f"pce{dnum:02d}*"))), None)
    if folder is None:
        print(f"  [dyad {dnum:02d}] no EEG folder found — skipping")
        return False

    print(f"[dyad {dnum:02d}] {folder.name}: loading task trials ...")
    t0 = time.time()
    raws, loaded = [], []
    for p in C.PARTICIPANTS:
        raw, ld = load_participant_task(folder, dnum, p, limit_trials)
        if raw is None:
            print(f"  P{p}: no trial files — skipping dyad")
            return False
        rms_uv = float(raw.get_data().std() * 1e6)
        print(f"  P{p}: {len(ld)} trials, {raw.n_times} samples, RMS ~{rms_uv:.1f} uV")
        raws.append(raw)
        loaded.append(ld)

    print("  running preprocess_dyad (notch->bandpass->bad-ch->ref->resample->epochs->AR) ...")
    result = preprocess_dyad(
        raws,
        line_freq=C.LINE_FREQ,
        bandpass=C.BANDPASS,
        sfreq_target=C.SFREQ_TARGET,
        reref=C.REREF,
        bad_ch_z=C.BAD_CHANNEL_ZSCORE,
        bad_ch_lof_threshold=C.BAD_CHANNEL_LOF_THRESHOLD,
        bad_ch_max_frac=C.BAD_CHANNEL_MAX_FRAC,
        epoch_len_s=C.EPOCH_LEN_S,
        epoch_overlap_s=C.EPOCH_OVERLAP_S,
        run_autoreject=run_ar,
        ar_strategy="union",
    )

    outdir = C.PREPROC_EEG / f"pce{dnum:02d}"
    outdir.mkdir(parents=True, exist_ok=True)
    for i, p in enumerate(C.PARTICIPANTS):
        result.raws[i].save(outdir / f"pce{dnum:02d}_P{p}_task-raw.fif", overwrite=True)
        result.epochs[i].save(outdir / f"pce{dnum:02d}_P{p}_task-epo.fif", overwrite=True)

    prov = {
        "dyad": dnum,
        "folder": folder.name,
        "participants": list(C.PARTICIPANTS),
        "trials_loaded": {f"P{C.PARTICIPANTS[i]}": loaded[i] for i in range(len(loaded))},
        "limit_trials": limit_trials,
        "rms_uv_raw": {
            f"P{C.PARTICIPANTS[i]}": float(raws[i].get_data().std() * 1e6)
            for i in range(len(raws))
        },
        "elapsed_s": round(time.time() - t0, 1),
        **result.provenance,
    }
    with open(outdir / f"pce{dnum:02d}_provenance.json", "w") as fh:
        json.dump(_jsonable(prov), fh, indent=2)

    ar = result.provenance.get("autoreject", {})
    print(
        f"  done in {prov['elapsed_s']}s | epochs {result.provenance.get('n_epochs_pre_ar')}"
        f"->{result.provenance.get('n_epochs_post_ar')} | bad ch {result.provenance.get('bad_channels')}"
        f" | AR {ar}"
    )
    print(f"  saved to {outdir}")
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dyad", type=int, action="append", help="dyad number (repeatable)")
    ap.add_argument("--all", action="store_true", help="process every dyad with EEG")
    ap.add_argument("--limit-trials", type=int, default=None, help="use only the first N trials (smoke test)")
    ap.add_argument("--no-autoreject", action="store_true", help="skip the autoreject step")
    args = ap.parse_args()

    C.ensure_output_dirs()

    if args.all:
        dyads = sorted(
            int(p.name[3:5]) for p in C.RAW_EEG.glob("pce*") if p.is_dir()
        )
        dyads = [d for d in dyads if d not in C.EXCLUDE_DYADS]
    elif args.dyad:
        dyads = args.dyad
    else:
        dyads = [1]
        print("No --dyad/--all given; defaulting to pilot dyad 1.")

    print(f"Dyads to process: {dyads}")
    ok = 0
    for d in dyads:
        try:
            ok += int(process_dyad(d, args.limit_trials, run_ar=not args.no_autoreject))
        except Exception as exc:
            print(f"  [dyad {d:02d}] ERROR: {type(exc).__name__}: {exc}")
    print(f"\nProcessed {ok}/{len(dyads)} dyads successfully.")


if __name__ == "__main__":
    main()
