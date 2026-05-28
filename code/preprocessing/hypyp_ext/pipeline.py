"""Continuous-Raw, dyad-aware preprocessing orchestrator with provenance.

HyPyP's worked examples start from data that has already been loaded, filtered,
epoched, and (often) ICA-cleaned by hand. This module provides the missing
front-end: a single call that takes a dyad's two *continuous* recordings and
returns analysis-ready, channel-aligned cleaned Raw **and** epochs, plus a
provenance dict recording every parameter and outcome. It reuses HyPyP for the
steps HyPyP already does well (:func:`hypyp.prep.filt`,
:func:`hypyp.prep.AR_local`) and ``hypyp_ext`` for the gaps.

This is the "moderate, no-ICA" pipeline:
    notch (conditional) -> band-pass -> bad-channel interpolation ->
    average reference -> resample -> fixed-length epochs -> autoreject.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import mne
from hypyp import prep

from .line_noise import notch_filter
from .channels import (
    detect_bad_channels,
    interpolate_bad_channels_dyad,
    harmonize_channels,
)


@dataclass
class DyadPreprocResult:
    """Output of :func:`preprocess_dyad`."""

    raws: List[mne.io.BaseRaw]          # cleaned continuous data, one per participant
    epochs: List[mne.Epochs]            # cleaned, channel-aligned epochs (AR'd if requested)
    provenance: dict = field(default_factory=dict)


def preprocess_dyad(
    raw_S: List[mne.io.BaseRaw],
    *,
    line_freq: Optional[float] = 60.0,
    bandpass: Tuple[Optional[float], Optional[float]] = (1.0, 40.0),
    sfreq_target: Optional[float] = 250.0,
    reref: str = "average",
    bad_ch_z: float = 4.0,
    bad_ch_method: str = "lof",
    bad_ch_lof_threshold: float = 2.5,
    bad_ch_max_frac: float = 0.15,
    epoch_len_s: float = 2.0,
    epoch_overlap_s: float = 0.0,
    run_autoreject: bool = True,
    ar_strategy: str = "union",
    ar_threshold: float = 50.0,
    verbose: bool = False,
) -> DyadPreprocResult:
    """Run the moderate (no-ICA) preprocessing pipeline on one dyad.

    Parameters
    ----------
    raw_S : list of mne.io.BaseRaw
        Exactly two continuous recordings (``[P1, P2]``) with a montage set.
        Inputs are never mutated.
    line_freq : float or None
        Mains frequency for the notch. Applied **only** when it is not already
        removed by the low-pass (i.e. when ``bandpass[1]`` is None or
        ``>= line_freq``). Pass None to disable.
    bandpass : (float or None, float or None)
        ``(l_freq, h_freq)`` for the FIR band-pass (via ``hypyp.prep.filt``).
    sfreq_target : float or None
        Resample target in Hz; None keeps the original rate.
    reref : str
        Reference scheme passed to ``set_eeg_reference`` (default ``"average"``).
    bad_ch_z, bad_ch_method : float, str
        Bad-channel detector settings (see :func:`hypyp_ext.channels.detect_bad_channels`).
    epoch_len_s, epoch_overlap_s : float
        Fixed-length epoching parameters.
    run_autoreject : bool
        If True, run ``hypyp.prep.AR_local`` (dyad-consistent autoreject).
    ar_strategy, ar_threshold : str, float
        Forwarded to ``AR_local``.
    verbose : bool
        Verbosity for MNE/HyPyP calls.

    Returns
    -------
    DyadPreprocResult
        ``raws`` (cleaned continuous), ``epochs`` (cleaned, aligned), and
        ``provenance`` (parameters + per-step outcomes).
    """
    if len(raw_S) != 2:
        raise ValueError("preprocess_dyad expects exactly two recordings (one dyad).")

    prov: dict = {
        "params": {
            "line_freq": line_freq,
            "bandpass": list(bandpass),
            "sfreq_in": float(raw_S[0].info["sfreq"]),
            "sfreq_target": sfreq_target,
            "reref": reref,
            "bad_ch_method": bad_ch_method,
            "bad_ch_z": bad_ch_z,
            "bad_ch_lof_threshold": bad_ch_lof_threshold,
            "bad_ch_max_frac": bad_ch_max_frac,
            "epoch_len_s": epoch_len_s,
            "epoch_overlap_s": epoch_overlap_s,
            "run_autoreject": run_autoreject,
            "ar_strategy": ar_strategy,
            "ar_threshold": ar_threshold,
        },
        "n_channels": len(raw_S[0].ch_names),
    }

    # Never mutate caller data.
    work = [r.copy() for r in raw_S]

    # 1) Line-noise notch — only if the low-pass doesn't already remove it.
    hp, lp = bandpass
    notch_applied = line_freq is not None and (lp is None or lp >= line_freq)
    if notch_applied:
        work = notch_filter(
            work, line_freq, include_harmonics_below=lp, verbose=verbose
        )
    prov["notch_applied"] = bool(notch_applied)
    if line_freq is not None and not notch_applied:
        prov["notch_note"] = (
            f"skipped: {line_freq} Hz is above the {lp} Hz low-pass "
            f"(already attenuated)"
        )

    # 2) Band-pass (HyPyP).
    work = prep.filt(work, freqs=(hp, lp))

    # 3) Bad-channel detection + interpolation (dyad-consistent channel set).
    bads_S = [
        detect_bad_channels(
            r,
            z_thresh=bad_ch_z,
            method=bad_ch_method,
            lof_threshold=bad_ch_lof_threshold,
            max_bad_frac=bad_ch_max_frac,
            verbose=verbose,
        )
        for r in work
    ]
    prov["bad_channels"] = {f"P{i+1}": list(b) for i, b in enumerate(bads_S)}
    work = interpolate_bad_channels_dyad(work, bads_S, reset_bads=True)

    # 4) Reference.
    for r in work:
        r.set_eeg_reference(reref, projection=False, verbose=verbose)

    # 5) Resample.
    if sfreq_target is not None and sfreq_target != work[0].info["sfreq"]:
        for r in work:
            r.resample(sfreq_target, verbose=verbose)
    prov["sfreq_out"] = float(work[0].info["sfreq"])

    # 6) Guarantee identical channel set/order across the dyad.
    work = harmonize_channels(work)
    prov["n_channels_out"] = len(work[0].ch_names)

    # 7) Fixed-length epochs (boundary-crossing epochs dropped by annotation).
    epochs_S = [
        mne.make_fixed_length_epochs(
            r,
            duration=epoch_len_s,
            overlap=epoch_overlap_s,
            preload=True,
            reject_by_annotation=True,
            verbose=verbose,
        )
        for r in work
    ]
    # Equalize epoch counts so AR_local's union indices are valid for both.
    n_eq = min(len(e) for e in epochs_S)
    epochs_S = [e[:n_eq] for e in epochs_S]
    prov["n_epochs_pre_ar"] = int(n_eq)

    # 8) Autoreject (HyPyP, dyad-consistent).
    if run_autoreject and n_eq > 0:
        try:
            epochs_S, dic_ar = prep.AR_local(
                epochs_S, strategy=ar_strategy, threshold=ar_threshold, verbose=False
            )
            prov["autoreject"] = {k: v for k, v in dic_ar.items()}
            prov["n_epochs_post_ar"] = int(len(epochs_S[0]))
        except Exception as exc:  # threshold exceeded or AR failure: keep pre-AR epochs
            prov["autoreject"] = {"error": f"{type(exc).__name__}: {exc}"}
            prov["n_epochs_post_ar"] = int(n_eq)
    else:
        prov["n_epochs_post_ar"] = int(n_eq)

    return DyadPreprocResult(raws=work, epochs=epochs_S, provenance=prov)
