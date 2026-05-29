"""Shared constants for the aperiodic-exponent analysis (Figure 2 Panel D).

Centralises the within-trial bin grid + 1/e-framework reference constants used
by :mod:`computeAperiodicExponent` and :mod:`fitExponentSensitivity`.

Mirrors the master-loop's ``atlas/_common.py`` so the spoke pipeline produces
byte-identical numbers when given byte-identical input data.
"""
from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_CODE = _HERE.parent
for _p in (str(_CODE), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import _config as C  # noqa: E402

# Within-trial bin grid: 4 s windows, 2 s step, over the [2, 58] s segment.
WIN_S = 4.0
STEP_S = 2.0
SEG = (2.0, 58.0)

# Population boot-up onset (from the EDA / sensitivity-dynamics framework)
# + nominal trial duration. The canonical 1/e marker is TAU + (T-TAU)/e.
TAU = 3.9
T_TRIAL = 60.0

# Spoke output directory for the exponent CSVs + sidecars.
OUT = C.RESULTS
