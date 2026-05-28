"""Power-line noise removal for a dyad of Raw recordings.

HyPyP's :func:`hypyp.prep.filt` applies a band-pass but offers no dedicated
line-noise step. When the analysis band extends above the mains frequency
(e.g. broadband or gamma work), an explicit notch is needed. This helper mirrors
``prep.filt``'s "operate on a list of two Raws" convention and leaves the inputs
untouched (returns copies).
"""
from __future__ import annotations

from typing import List, Optional, Union

import numpy as np
import mne


def notch_filter(
    raw_S: List[mne.io.BaseRaw],
    line_freq: float,
    include_harmonics_below: Optional[float] = None,
    picks: Union[str, list, None] = "eeg",
    verbose: bool = False,
) -> List[mne.io.BaseRaw]:
    """Notch out mains noise (and optional harmonics) from each Raw in a dyad.

    Parameters
    ----------
    raw_S : list of mne.io.BaseRaw
        One dyad: ``[participant_1, participant_2]`` (any length works).
    line_freq : float
        Mains frequency in Hz (e.g. 50 in most of the world, 60 in the
        Americas / Okinawa).
    include_harmonics_below : float or None, optional
        If given, also notch harmonics ``line_freq, 2*line_freq, ...`` strictly
        below this frequency (capped at the Nyquist). If None, only the
        fundamental is notched.
    picks : str | list | None, optional
        Channels to filter (default ``"eeg"``).
    verbose : bool, optional
        Passed through to MNE.

    Returns
    -------
    list of mne.io.BaseRaw
        Notch-filtered copies, in the same order as ``raw_S``.

    Notes
    -----
    With a low-pass cutoff *below* ``line_freq`` the mains component is already
    attenuated, making this step redundant; :func:`hypyp_ext.pipeline.preprocess_dyad`
    therefore applies it conditionally.
    """
    out: List[mne.io.BaseRaw] = []
    for raw in raw_S:
        r = raw.copy()
        nyq = r.info["sfreq"] / 2.0
        if include_harmonics_below is None:
            freqs = [line_freq]
        else:
            top = min(include_harmonics_below, nyq)
            freqs = list(np.arange(line_freq, top, line_freq))
        freqs = [f for f in freqs if 0.0 < f < nyq]
        if freqs:
            r.notch_filter(freqs=freqs, picks=picks, verbose=verbose)
        out.append(r)
    return out
