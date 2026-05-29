#!/usr/bin/env python3
"""Raw PCE EEG (.mat) → minimally cleaned, 250 Hz, dyad-aligned *-raw.fif.

The recording-level chain (Lerique et al. 2024 ECSU-PCE Dataset tech report,
p.7) already provides: Brain Vision Brainamp DC, 0.016-1000 Hz analog cutoff,
1000 Hz sample rate, FCz common reference, 64-ch actiCAP, **60 Hz notch
applied at the BrainVision Recorder**. So a software notch is redundant. Per
``docs/2026-05-29-spoke-preprocessing-recon.md`` in the master loop, the
minimum we need on top of that, for the per-(participant × bin) FOOOF
aperiodic-exponent analysis (Figure 2D), is:

    bandpass 1-40 Hz   → Nyquist alias safety (analog cutoff is at fs, so
                         500-1000 Hz folds into 0-500 Hz at acquisition) +
                         slow-drift control
    resample 250 Hz    → FOOOF-compatible rate (Welch nperseg = 2 × sf)
    bad-channel LOF + interpolate
    average reference  → the single most important step
                         (reproduces byte-identical to the master)

That's the whole chain. Layer-2 byte-identical to the master-loop's audit
reference (recon's ``minimal+reref+lof`` variant). No need for HyPyP,
autoreject, notch, or fixed-length epochs.

Outputs ``data/preprocessed/EEG/pceXX/pceXX_PY_task-raw.fif`` and a per-dyad
provenance JSON, also archived on Zenodo as part of ``preprocessed.zip`` for
readers who do not want to download the ~15 GB raw EEG from OSF.

Usage
-----
    # pilot: one dyad
    python code/preprocessing/preprocessEEGForExponent.py --dyad 1

    # quick smoke test: few trials only
    python code/preprocessing/preprocessEEGForExponent.py --dyad 1 --limit-trials 3

    # full run, all dyads with EEG
    python code/preprocessing/preprocessEEGForExponent.py --all
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import List, Sequence

import numpy as np
import scipy.io as sio
import mne

_HERE = Path(__file__).resolve().parent           # code/preprocessing
_CODE = _HERE.parent                              # code
for _p in (str(_CODE), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import _config as C  # noqa: E402

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


# --------------------------------------------------------------------------- #
#  Bad-channel detection + dyad-consistent interpolation
#  (small MNE wrappers — kept inline so the spoke has no preprocessing
#  subpackage; cap at BAD_CHANNEL_MAX_FRAC of channels to avoid the LOF
#  detector over-flagging on no-ICA broadband data).
# --------------------------------------------------------------------------- #
def detect_bad_channels(
    raw: mne.io.BaseRaw,
    z_thresh: float = 4.0,
    lof_threshold: float = 2.5,
    max_bad_frac: float = 0.15,
    verbose: bool = False,
) -> List[str]:
    """Flag outlier EEG channels (LOF; statistical fallback). Capped at
    ``max_bad_frac`` of the EEG channels."""
    picks = mne.pick_types(raw.info, eeg=True, exclude=[])
    eeg_names = [raw.ch_names[i] for i in picks]
    max_bad_n = max(1, int(round(max_bad_frac * len(picks))))

    def _cap(bads, severity):
        bads = list(bads)
        if len(bads) > max_bad_n:
            bads = sorted(bads, key=lambda c: severity.get(c, 0.0), reverse=True)[:max_bad_n]
        return bads

    try:
        res = mne.preprocessing.find_bad_channels_lof(
            raw, picks="eeg", threshold=lof_threshold,
            return_scores=True, verbose=verbose,
        )
        bads, scores = res if isinstance(res, tuple) else (res, None)
        severity = (
            {n: float(s) for n, s in zip(eeg_names, np.asarray(scores).ravel())}
            if scores is not None else {}
        )
        return _cap(bads, severity)
    except Exception:
        pass  # fall through to the statistical detector

    data = raw.get_data(picks=picks)  # (n_ch, n_times)

    def _robust_z(x: np.ndarray) -> np.ndarray:
        med = np.median(x)
        mad = np.median(np.abs(x - med)) or 1e-30
        return 0.6745 * (x - med) / mad

    z_var = _robust_z(np.log(np.var(data, axis=1) + 1e-30))
    corr = np.corrcoef(data)
    np.fill_diagonal(corr, np.nan)
    z_corr = _robust_z(np.nanmean(np.abs(corr), axis=1))
    bad_mask = (np.abs(z_var) > z_thresh) | (z_corr < -z_thresh)
    severity = {eeg_names[i]: max(abs(z_var[i]), abs(z_corr[i])) for i in range(len(eeg_names))}
    return _cap([eeg_names[i] for i in np.where(bad_mask)[0]], severity)


def interpolate_bad_channels_dyad(
    raw_S: List[mne.io.BaseRaw],
    bads_S: Sequence[Sequence[str]],
) -> List[mne.io.BaseRaw]:
    """Interpolate each participant's bad channels (spherical splines), keeping
    the channel set intact so the dyad stays aligned."""
    out: List[mne.io.BaseRaw] = []
    for raw, bads in zip(raw_S, bads_S):
        r = raw.copy()
        r.info["bads"] = list(bads)
        if bads:
            r.interpolate_bads(reset_bads=True)
        out.append(r)
    return out


# --------------------------------------------------------------------------- #
#  Minimal preprocessing — the chain
# --------------------------------------------------------------------------- #
def minimal_preprocess(raws: List[mne.io.BaseRaw]) -> tuple[List[mne.io.BaseRaw], dict]:
    """Apply the minimal preprocessing chain to one dyad and return
    (raws_out, provenance dict)."""
    prov = {
        "chain": [
            f"bandpass {C.BANDPASS[0]:.1f}-{C.BANDPASS[1]:.1f} Hz (FIR)",
            f"resample {C.SFREQ_TARGET:.0f} Hz",
            f"bad-channel LOF + interpolate (cap {C.BAD_CHANNEL_MAX_FRAC:.0%})",
            f"average reference",
        ],
        "params": {
            "bandpass_hz": list(C.BANDPASS),
            "sfreq_target_hz": C.SFREQ_TARGET,
            "bad_ch_lof_threshold": C.BAD_CHANNEL_LOF_THRESHOLD,
            "bad_ch_max_frac": C.BAD_CHANNEL_MAX_FRAC,
            "reref": C.REREF,
        },
        "recording_chain_baked_in": (
            "Brain Vision Brainamp DC, 0.016-1000 Hz analog cutoff, 1000 Hz sample "
            "rate, FCz common reference, 64-ch actiCAP, 60 Hz notch applied at the "
            "BrainVision Recorder (Lerique et al. 2024, ECSU-PCE Dataset tech report, p.7). "
            "No software notch needed."
        ),
        "sfreq_in_hz": float(raws[0].info["sfreq"]),
    }

    work: List[mne.io.BaseRaw] = []
    for r in raws:
        r2 = r.copy().filter(C.BANDPASS[0], C.BANDPASS[1],
                             fir_design="firwin", verbose=False)
        r2.resample(C.SFREQ_TARGET, verbose=False)
        work.append(r2)

    bads_S = [
        detect_bad_channels(
            r,
            z_thresh=C.BAD_CHANNEL_ZSCORE,
            lof_threshold=C.BAD_CHANNEL_LOF_THRESHOLD,
            max_bad_frac=C.BAD_CHANNEL_MAX_FRAC,
            verbose=False,
        )
        for r in work
    ]
    prov["bad_channels"] = {f"P{i+1}": list(b) for i, b in enumerate(bads_S)}
    work = interpolate_bad_channels_dyad(work, bads_S)

    for r in work:
        r.set_eeg_reference(C.REREF, projection=False, verbose=False)

    prov["sfreq_out_hz"] = float(work[0].info["sfreq"])
    prov["n_channels_out"] = len(work[0].ch_names)
    return work, prov


def _jsonable(obj):
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
def process_dyad(dnum: int, limit_trials: int | None) -> bool:
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

    print("  running minimal preprocessing (bandpass → resample → LOF → reref) ...")
    raws_out, chain_prov = minimal_preprocess(raws)

    outdir = C.PREPROC_EEG / f"pce{dnum:02d}"
    outdir.mkdir(parents=True, exist_ok=True)
    for i, p in enumerate(C.PARTICIPANTS):
        raws_out[i].save(outdir / f"pce{dnum:02d}_P{p}_task-raw.fif", overwrite=True)

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
        **chain_prov,
    }
    with open(outdir / f"pce{dnum:02d}_provenance.json", "w") as fh:
        json.dump(_jsonable(prov), fh, indent=2)

    print(
        f"  done in {prov['elapsed_s']}s | "
        f"sfreq {prov['sfreq_in_hz']:.0f}→{prov['sfreq_out_hz']:.0f} Hz | "
        f"bad ch {prov['bad_channels']}"
    )
    print(f"  saved to {outdir}")
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dyad", type=int, action="append", help="dyad number (repeatable)")
    ap.add_argument("--all", action="store_true", help="process every dyad with EEG")
    ap.add_argument("--limit-trials", type=int, default=None, help="use only the first N trials (smoke test)")
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
            ok += int(process_dyad(d, args.limit_trials))
        except Exception as exc:
            print(f"  [dyad {d:02d}] ERROR: {type(exc).__name__}: {exc}")
    print(f"\nProcessed {ok}/{len(dyads)} dyads successfully.")


if __name__ == "__main__":
    main()
