"""Shared Python configuration for the pce-sensitivity-dynamics spoke.

Paths are resolved relative to the repo root so scripts work regardless of the
working directory. Acquisition + preprocessing constants live here so the EEG
preprocessing, FOOOF analysis, and figure code stay in sync.

The aperiodic-exponent pipeline (Figure 2 Panel D) uses every constant in this
module. The existing MATLAB preprocessing scripts (computeGSPStats etc.) carry
their own constants — keep this module focused on the Python pipeline.
"""
from __future__ import annotations

from pathlib import Path

# --- Paths --------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = REPO_ROOT / "code"

DATA_RAW = REPO_ROOT / "data" / "raw"
DATA_PREPROC = REPO_ROOT / "data" / "preprocessed"
RESULTS = REPO_ROOT / "results"

RAW_EEG = DATA_RAW / "EEG"
PREPROC_EEG = DATA_PREPROC / "EEG"
PREPROC_CLICKPAS = DATA_PREPROC / "ClickTimes"

# --- Acquisition --------------------------------------------------------
SFREQ_ORIG = 1000.0
N_CHAN = 64
TRIAL_DUR_S = 60
PARTICIPANTS = (1, 2)
TRIALS = tuple(range(1, 19))
EXCLUDE_DYADS = (31,)

# 64-channel actiCAP label order = row order of the raw .mat matrices.
# Same list as the existing MATLAB scripts (preprocessGSP.m, extractAllChannels.m).
CH_NAMES = [
    "Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8", "FC5",
    "FC1", "FC2", "FC6", "T7", "C3", "Cz", "C4", "T8",
    "TP9", "CP5", "CP1", "CP2", "CP6", "TP10", "P7", "P3",
    "Pz", "P4", "P8", "PO9", "O1", "Oz", "O2", "PO10",
    "AF7", "AF3", "AF4", "AF8", "F5", "F1", "F2", "F6",
    "FT9", "FT7", "FC3", "FC4", "FT8", "FT10", "C5", "C1",
    "C2", "C6", "TP7", "CP3", "CPz", "CP4", "TP8", "P5",
    "P1", "P2", "P6", "PO7", "PO3", "POz", "PO4", "PO8",
]
assert len(CH_NAMES) == N_CHAN

# --- Preprocessing (moderate cleaning, no ICA) -------------------------
MONTAGE = "standard_1005"
EEG_UNIT_SCALE = 1e-6        # raw .mat values are microvolts -> Volts
LINE_FREQ = 60.0             # Okinawa 60 Hz grid
BANDPASS = (1.0, 40.0)       # FIR via MNE Raw.filter (Nyquist alias safety + drift)
SFREQ_TARGET = 250.0
REREF = "average"

BAD_CHANNEL_ZSCORE = 4.0
BAD_CHANNEL_LOF_THRESHOLD = 2.5
BAD_CHANNEL_MAX_FRAC = 0.15

EPOCH_LEN_S = 2.0
EPOCH_OVERLAP_S = 0.0


def ensure_output_dirs() -> None:
    """Create the preprocessed/results output directories if missing."""
    for d in (PREPROC_EEG, RESULTS):
        d.mkdir(parents=True, exist_ok=True)
