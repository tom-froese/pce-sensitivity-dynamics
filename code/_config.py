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
# This is the dataset's .ced ANATOMICAL montage (AFz/Iz) — NOT the BrainProducts
# manufacturer-default template (PO9/PO10). The default was wrong in ORDER (61/64 rows)
# and MEMBERSHIP (PO9/PO10 vs the cap's real AFz/Iz), which scrambled every per-channel /
# topographic label built from these names (e.g. the Zenodo *_task-raw.fif). Verified end
# to end against the raw BrainVision montage CMA-64_REF_HS2.bvef + blink physiology (AFz).
# Kept in lockstep with the MATLAB scripts (preprocessGSP.m, extractAllChannels.m), which
# were corrected 2026-07-13 (commit 5e18466). See
# pce-master-loop/docs/2026-07-13-montage-cap-provenance.md.
CH_NAMES = [
    "Fp1", "Fp2", "AF7", "AF3", "AFz", "AF4", "AF8", "F7",
    "F5", "F3", "F1", "Fz", "F2", "F4", "F6", "F8",
    "FT9", "FT7", "FC5", "FC3", "FC1", "FC2", "FC4", "FC6",
    "FT8", "FT10", "T7", "C5", "C3", "C1", "Cz", "C2",
    "C4", "C6", "T8", "TP9", "TP7", "CP5", "CP3", "CP1",
    "CPz", "CP2", "CP4", "CP6", "TP8", "TP10", "P7", "P5",
    "P3", "P1", "Pz", "P2", "P4", "P6", "P8", "PO7",
    "PO3", "POz", "PO4", "PO8", "O1", "Oz", "O2", "Iz",
]
assert len(CH_NAMES) == N_CHAN
# GUARDRAIL: this cap recorded AFz/Iz, never PO9/PO10 — fail loud if the order ever
# regresses to a manufacturer-default template (mirrors the extractAllChannels.m assert).
assert (
    "AFz" in CH_NAMES and "Iz" in CH_NAMES
    and "PO9" not in CH_NAMES and "PO10" not in CH_NAMES
), (
    "CH_NAMES must be the .ced anatomical order with AFz/Iz (not PO9/PO10). "
    "See pce-master-loop/docs/2026-07-13-montage-cap-provenance.md."
)

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
