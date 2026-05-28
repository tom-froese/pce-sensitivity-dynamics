"""Dyad-consistent bad-channel handling.

Inter-brain connectivity (HyPyP ``analyses``) requires the two members of a dyad
to carry an **identical channel set in identical order**. HyPyP assumes this is
already true; in practice each recording has its own bad channels. These helpers
detect bad channels automatically, interpolate them (which preserves the channel
set), and guarantee the dyad ends channel-aligned.
"""
from __future__ import annotations

from typing import List, Sequence

import numpy as np
import mne


def detect_bad_channels(
    raw: mne.io.BaseRaw,
    z_thresh: float = 4.0,
    method: str = "lof",
    lof_threshold: float = 2.5,
    max_bad_frac: float = 0.15,
    verbose: bool = False,
) -> List[str]:
    """Automatically flag bad EEG channels in a single recording.

    The result is **capped** at ``max_bad_frac`` of the EEG channels: if more
    are flagged, only the most extreme are kept. A high raw count usually means
    the detector is reacting to a global artifact (e.g. ocular activity in
    no-ICA pipelines) rather than truly broken sensors — interpolating a whole
    region would do more harm than good, so that case is reined in and left to
    autoreject instead.

    Parameters
    ----------
    raw : mne.io.BaseRaw
        Continuous recording (montage already set, for the LOF method).
    z_thresh : float, optional
        Robust z-score cutoff for the fallback statistical detector (default 4.0).
    method : {"lof", "zscore"}, optional
        ``"lof"`` uses MNE's Local Outlier Factor detector when available and
        falls back to ``"zscore"`` otherwise.
    lof_threshold : float, optional
        LOF score cutoff (default 2.5). MNE's own default (1.5) over-flags
        broadband data; 2.5 keeps only clearly extreme channels.
    max_bad_frac : float, optional
        Maximum fraction of EEG channels that may be returned (default 0.15).
    verbose : bool, optional
        Passed through to MNE where applicable.

    Returns
    -------
    list of str
        Names of channels judged bad (possibly empty), worst-first when capped.
    """
    picks = mne.pick_types(raw.info, eeg=True, exclude=[])
    eeg_names = [raw.ch_names[i] for i in picks]
    max_bad_n = max(1, int(round(max_bad_frac * len(picks))))

    def _cap(bads, severity):
        bads = list(bads)
        if len(bads) > max_bad_n:
            bads = sorted(bads, key=lambda c: severity.get(c, 0.0), reverse=True)[:max_bad_n]
        return bads

    if method == "lof":
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
    reset_bads: bool = True,
) -> List[mne.io.BaseRaw]:
    """Interpolate each participant's bad channels (spherical splines).

    Interpolation rebuilds flagged channels from their neighbours and keeps the
    full channel set, so the dyad stays channel-aligned regardless of which
    channels were bad in whom.

    Parameters
    ----------
    raw_S : list of mne.io.BaseRaw
        One dyad. A montage with sensor positions must be set.
    bads_S : sequence of sequence of str
        Per-participant bad-channel name lists (same length/order as ``raw_S``).
    reset_bads : bool, optional
        Clear ``info['bads']`` after interpolation (default True).

    Returns
    -------
    list of mne.io.BaseRaw
        Interpolated copies.
    """
    out: List[mne.io.BaseRaw] = []
    for raw, bads in zip(raw_S, bads_S):
        r = raw.copy()
        r.info["bads"] = list(bads)
        if bads:
            r.interpolate_bads(reset_bads=reset_bads)
        out.append(r)
    return out


def harmonize_channels(raw_S: List[mne.io.BaseRaw]) -> List[mne.io.BaseRaw]:
    """Force every Raw in the dyad onto the same channel set and order.

    Drops any channel not present in *all* members, then reorders to a common
    layout. After interpolation this is normally a no-op, but it is a cheap
    guarantee that downstream inter-brain estimators receive aligned data.

    Parameters
    ----------
    raw_S : list of mne.io.BaseRaw
        One dyad.

    Returns
    -------
    list of mne.io.BaseRaw
        Channel-harmonized copies.
    """
    common = [ch for ch in raw_S[0].ch_names if all(ch in r.ch_names for r in raw_S)]
    return [r.copy().reorder_channels(common) for r in raw_S]
