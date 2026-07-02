#!/usr/bin/env python3
"""EDA exponential-versus-linear decisive test (free-lambda robustness).

The main behavioral figure (``plotFigure1_Behavioral.m``, Panel C) fits the
within-trial electrodermal-activity (EDA) trajectory to the reliability decay
``R(x) = A0*exp(-e*x) + B`` with the rate **fixed** at lambda = e and reports
R^2 ~ 0.97. This script adds the manuscript's *decisive* robustness test on the
same trajectory:

  * free the rate: fit ``A*exp(-lambda*x) + B`` with lambda a free parameter,
  * compare the fixed-``lambda=e`` exponential against a straight line.

The trajectory is built exactly as in the master-loop's within-trial analysis
(``atlas/within_trial_complexity.py::eda_analysis``): the task EDA is windowed
into 4 s windows at a 2 s step over the [2, 58] s segment (bin centers 4..56 s),
each trial's windowed trace is z-scored across windows, and the traces are
cohort-averaged. The fixed-``lambda=e`` and free-``lambda`` exponentials are fit
in the reduced coordinate ``x = (t - tau)/(T - tau)`` with tau = 3.9 s, T = 60 s.

Reported (matches the manuscript, SI Fig. 1C robustness):
  * free-``lambda`` ~ 2.62 (~ e), i.e. freeing the rate recovers e,
  * ``dR2 = R2[exp@e] - R2[linear] = +0.07`` (the exponential decisively beats a
    straight line).

Reads the bundled ``data/preprocessed/EDA/EDA_Task_Preprocessed.csv`` (25 Hz,
1500 samples per 60 s trial); writes ``results/eda_free_lambda.csv``.

Usage
-----
    python code/analysis/computeEDAFreeLambda.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

_HERE = Path(__file__).resolve().parent
_CODE = _HERE.parent
for _p in (str(_CODE), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import _config as C  # noqa: E402

# Trajectory grid (identical to the within-trial complexity/EDA analysis).
WIN_S, STEP_S = 4.0, 2.0
SEG = (2.0, 58.0)
TAU, T_TRIAL = 3.9, 60.0
FS_EDA = 25.0                        # preprocessed task EDA sample rate (Hz)


def _r2(y, pred):
    return 1.0 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)


def eda_trajectory():
    """Cohort-averaged, per-trial-z-scored within-trial EDA trajectory.

    Returns (centers_s, traj_z, n_trials).
    """
    csv = C.DATA_PREPROC / "EDA" / "EDA_Task_Preprocessed.csv"
    if not csv.is_file():
        raise SystemExit(f"missing {csv} — run preprocessEDA.m first (see README)")
    eda = pd.read_csv(csv)
    L = int(round(WIN_S * FS_EDA))
    step = int(round(STEP_S * FS_EDA))
    seg_lo, seg_hi = int(SEG[0] * FS_EDA), int(SEG[1] * FS_EDA)
    acc = []
    for _, g in eda.groupby(["DyadID", "ParticipantID", "TrialNum"], sort=False):
        v = g.sort_values("Time_s").EDA_uS.values
        if len(v) < seg_hi:
            continue
        seg = v[seg_lo:seg_hi]
        starts = range(0, len(seg) - L + 1, step)
        acc.append(np.array([seg[s:s + L].mean() for s in starts]))
    arr = np.asarray(acc)
    n_win = arr.shape[1]
    centers = np.array([SEG[0] + i * STEP_S + WIN_S / 2 for i in range(n_win)])
    # per-trial z-score across windows, then cohort-average (as in the master).
    traj = ((arr - arr.mean(1, keepdims=True)) / (arr.std(1, keepdims=True) + 1e-9)).mean(0)
    return centers, traj, arr.shape[0]


def analyze():
    centers, traj, n_trials = eda_trajectory()
    x = (centers - TAU) / (T_TRIAL - TAU)

    # fixed lambda = e (linear least squares on the exp basis)
    c = np.where(x > 0, np.exp(-np.e * x), 1.0)
    Ae = np.column_stack([c, np.ones_like(c)])
    coef_e, *_ = np.linalg.lstsq(Ae, traj, rcond=None)
    r2_e = _r2(traj, Ae @ coef_e)

    # straight line (in trial time)
    Al = np.column_stack([centers, np.ones_like(centers)])
    coef_l, *_ = np.linalg.lstsq(Al, traj, rcond=None)
    r2_lin = _r2(traj, Al @ coef_l)

    # free lambda: A*exp(-lambda*x) + B on the z-scored trajectory
    tz = (traj - traj.mean()) / traj.std()

    def model(xx, A, lam, B):
        return A * np.exp(-lam * xx) + B

    popt, pcov = curve_fit(model, x, tz, p0=[2.0, np.e, -1.0], maxfev=10000)
    lam_free = float(popt[1])
    lam_se = float(np.sqrt(np.diag(pcov))[1])
    r2_free = _r2(tz, model(x, *popt))

    d_r2 = r2_e - r2_lin
    print("=== EDA exponential-versus-linear decisive test ===")
    print(f"  n trials              = {n_trials}")
    print(f"  R^2  exp @ lambda=e   = {r2_e:.3f}   (Panel C reliability decay)")
    print(f"  R^2  linear           = {r2_lin:.3f}")
    print(f"  dR^2 (exp@e - linear) = {d_r2:+.3f}   (manuscript: +0.07)")
    print(f"  free lambda           = {lam_free:.2f} +/- {lam_se:.2f}   "
          f"(manuscript: 2.62 ~ e = {np.e:.3f})")

    out = C.RESULTS / "eda_free_lambda.csv"
    C.RESULTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([dict(
        n_trials=n_trials,
        tau_s=TAU, T_s=T_TRIAL, win_s=WIN_S, step_s=STEP_S,
        seg_lo_s=SEG[0], seg_hi_s=SEG[1],
        R2_exp_at_e=r2_e, R2_linear=r2_lin, delta_R2_exp_minus_linear=d_r2,
        free_lambda=lam_free, free_lambda_se=lam_se, R2_free=r2_free,
    )]).to_csv(out, index=False)
    print(f"  saved -> {out.relative_to(C.REPO_ROOT)}")
    return d_r2, lam_free


if __name__ == "__main__":
    analyze()
