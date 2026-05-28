"""hypyp_ext — extensions to HyPyP for dyad-aware continuous-EEG preprocessing.

HyPyP's :mod:`hypyp.prep` is epochs-centric: it assumes already-loaded data and
provides bandpass filtering, ICA, and dyad-consistent autoreject (``AR_local``).
A few steps needed for a rigorous, reproducible *continuous-data* pipeline are
missing. This package fills those gaps with functions written to HyPyP's API
conventions (operate on a 2-element list ``raw_S`` / ``epochs_S`` = one dyad,
numpydoc docstrings, no side effects on inputs) so they can be proposed upstream.

Public API
----------
- :func:`line_noise.notch_filter`            — power-line notch (gap in ``prep``)
- :func:`channels.detect_bad_channels`       — automatic bad-channel detection
- :func:`channels.interpolate_bad_channels_dyad`
- :func:`channels.harmonize_channels`        — dyad-consistent channel set
- :func:`pipeline.preprocess_dyad`           — continuous-Raw orchestrator + provenance
"""
from . import line_noise, channels, pipeline
from .line_noise import notch_filter
from .channels import (
    detect_bad_channels,
    interpolate_bad_channels_dyad,
    harmonize_channels,
)
from .pipeline import preprocess_dyad, DyadPreprocResult

__all__ = [
    "line_noise",
    "channels",
    "pipeline",
    "notch_filter",
    "detect_bad_channels",
    "interpolate_bad_channels_dyad",
    "harmonize_channels",
    "preprocess_dyad",
    "DyadPreprocResult",
]
