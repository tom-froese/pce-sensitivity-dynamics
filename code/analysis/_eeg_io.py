"""Dataset + loaders for the spoke's Python EEG analyses.

Mirrors the loader API used by the master-loop's ``awareness.py``:

    Dataset(eeg_root=..., raw_suffix="task-raw")
      .dyad_dir(dyad) -> Path                       # preprocessed/EEG/pceXX
      .raw_path(dyad, p) -> Path                    # ..._task-raw.fif
      .provenance_path(dyad) -> Path                # ..._provenance.json

    load_dyad_raws(ds, dyad) -> {p: mne.Raw}        # channel-ordered to CH_NAMES
    load_provenance(ds, dyad) -> dict
    segment_onsets(raw, trials_loaded) -> {trial: onset_s}

The spoke targets only the noICA pipeline (task-raw), so unlike the master we
do not parameterise across pipeline variants — the Dataset defaults are the
only thing this module needs to do.
"""
from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import mne

_HERE = Path(__file__).resolve().parent
_CODE = _HERE.parent
for _p in (str(_CODE), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import _config as C  # noqa: E402


@dataclass
class Dataset:
    """Paths for the noICA EEG dataset; defaults to this repo."""

    eeg_root: Path = C.PREPROC_EEG
    name: str = "spoke"
    raw_suffix: str = "task-raw"

    def dyad_dir(self, dyad: int) -> Path:
        return self.eeg_root / f"pce{dyad:02d}"

    def raw_path(self, dyad: int, p: int) -> Path:
        return self.dyad_dir(dyad) / f"pce{dyad:02d}_P{p}_{self.raw_suffix}.fif"

    def provenance_path(self, dyad: int) -> Path:
        return self.dyad_dir(dyad) / f"pce{dyad:02d}_provenance.json"


def default_dataset() -> Dataset:
    return Dataset()


def _reorder_to_canonical(raw: mne.io.BaseRaw) -> mne.io.BaseRaw:
    """Force channel order to ``C.CH_NAMES`` so brains are cross-aligned."""
    missing = [ch for ch in C.CH_NAMES if ch not in raw.ch_names]
    if missing:
        raise ValueError(f"raw missing channels {missing}")
    return raw.copy().reorder_channels(list(C.CH_NAMES))


def load_dyad_raws(ds: Dataset, dyad: int) -> Dict[int, mne.io.BaseRaw]:
    """Load both participants' cleaned continuous Raws for one dyad."""
    raws = {}
    for p in C.PARTICIPANTS:
        raw = mne.io.read_raw_fif(ds.raw_path(dyad, p), preload=True, verbose=False)
        raws[p] = _reorder_to_canonical(raw)
    return raws


def load_provenance(ds: Dataset, dyad: int) -> dict:
    return json.loads(ds.provenance_path(dyad).read_text())


def segment_onsets(raw: mne.io.BaseRaw, trials_loaded: List[int]) -> Dict[int, float]:
    """Map trial number -> onset (s) on the concatenated timeline.

    Uses position in ``trials_loaded`` * TRIAL_DUR_S, validated against the .fif
    boundary annotations (exact 60 s segments).
    """
    onsets = {t: i * C.TRIAL_DUR_S for i, t in enumerate(trials_loaded)}
    # validate against boundary annotations
    bnd = sorted({round(o, 3) for o, d in
                  zip(raw.annotations.onset, raw.annotations.description)
                  if "boundary" in str(d)})
    expected = [i * C.TRIAL_DUR_S for i in range(1, len(trials_loaded))]
    for e in expected:
        if not any(abs(e - b) < 0.5 for b in bnd):
            raise AssertionError(f"missing expected boundary near {e}s (got {bnd})")
    dur = raw.times[-1]
    if abs(dur + 1.0 / raw.info["sfreq"] - len(trials_loaded) * C.TRIAL_DUR_S) > 1.0:
        raise AssertionError(
            f"recording duration {dur:.1f}s != {len(trials_loaded)} x {C.TRIAL_DUR_S}s"
        )
    return onsets
