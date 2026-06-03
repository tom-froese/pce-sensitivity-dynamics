#!/usr/bin/env python3
"""Reproduction-validation gate for Froese (2026) PNAS Brief Report.

Unlike a smoke test, this *gates*: it asserts that the values in the
committed result files match the published numbers within tolerance, and
exits non-zero (failing CI / an interactive run) if any check fails.

It reads only existing outputs and never regenerates or modifies anything:
  - results/Figure2_perchannel_fits.csv          (per-channel S(x) fits)
  - results/Figure3_crossover_stats.csv          (PAS 4/3 logistic crossover)
  - results/aperiodic_exponent_within_trial.json (1/f exponent S(x) fit)
  - data/preprocessed/EEG/globalScalpPotential_stats.mat  (GSP grand fit)

Usage:
    python3 validate_reproduction.py            # validate this checkout
    python3 validate_reproduction.py --root DIR # validate a Zenodo unzip elsewhere

Exit status: 0 if every check passes, 1 otherwise.

Design notes
------------
The expected values below are the *published* numbers. Reproduction is
seeded and deterministic, so tolerances are tight; they exist only to absorb
cross-platform floating-point and BLAS differences, not genuine drift. Two
checks are deliberately discriminating rather than cosmetic:
  * the rest-block control R^2 must be LOW (the S(x) fit is task-specific, not
    an artifact that fits any time series), and
  * the bootstrap lambda CI must bracket e (the paper's central claim).
A reproduction that silently broke one of these would fail here.
"""

import argparse
import json
import math
import os
import sys

import h5py
import numpy as np
import pandas as pd

E = math.e  # 2.718281828...


# --------------------------------------------------------------------------
#  Tiny check harness
# --------------------------------------------------------------------------
class Checks:
    def __init__(self):
        self.rows = []  # (name, ok, detail)

    def near(self, name, actual, expected, atol):
        ok = actual is not None and abs(float(actual) - expected) <= atol
        self.rows.append((name, ok,
                          f"{actual} vs {expected} (±{atol})"))

    def cond(self, name, ok, detail):
        self.rows.append((name, bool(ok), detail))

    def report(self):
        width = max(len(n) for n, _, _ in self.rows)
        print("\n  Reproduction validation — Froese (2026) sensitivity dynamics")
        print("  " + "=" * (width + 30))
        n_fail = 0
        for name, ok, detail in self.rows:
            tag = "PASS" if ok else "FAIL"
            if not ok:
                n_fail += 1
            print(f"  [{tag}] {name.ljust(width)}  {detail}")
        print("  " + "=" * (width + 30))
        total = len(self.rows)
        if n_fail:
            print(f"  RESULT: {n_fail}/{total} checks FAILED — reproduction does NOT match.\n")
        else:
            print(f"  RESULT: all {total} checks passed — reproduction matches the paper.\n")
        return n_fail == 0


# --------------------------------------------------------------------------
#  Loaders
# --------------------------------------------------------------------------
def _scalar(f, key):
    """Read a top-level scalar dataset from a MATLAB v7.3 (HDF5) file."""
    return float(np.array(f[key]).flatten()[0])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", default=os.path.dirname(os.path.abspath(__file__)),
                    help="repo/bundle root (default: this script's directory)")
    args = ap.parse_args()
    root = args.root

    results = os.path.join(root, "results")
    mat_path = os.path.join(root, "data", "preprocessed", "EEG",
                            "globalScalpPotential_stats.mat")

    required = [
        os.path.join(results, "Figure2_perchannel_fits.csv"),
        os.path.join(results, "Figure3_crossover_stats.csv"),
        os.path.join(results, "aperiodic_exponent_within_trial.json"),
        mat_path,
    ]
    missing = [p for p in required if not os.path.isfile(p)]
    if missing:
        print("  MISSING required result files — cannot validate:")
        for p in missing:
            print(f"    {os.path.relpath(p, root)}")
        print("\n  Download from https://doi.org/10.5281/zenodo.19425014 "
              "or run reproduce.m first.")
        return 1

    c = Checks()

    # --- Global scalp potential: the grand S(x) fit + permutation test ------
    with h5py.File(mat_path, "r") as f:
        grand_r2 = _scalar(f, "grandR2")
        p_perm = _scalar(f, "pPerm")
        lam = _scalar(f, "lambda")
        n_part = _scalar(f, "nPart")
        trough = _scalar(f, "grandTrough")
        rest_r2 = _scalar(f, "grandRestR2")

    c.near("GSP grand-average R^2", grand_r2, 0.8864, 0.005)
    c.cond("GSP permutation p significant", p_perm <= 0.005, f"pPerm={p_perm}")
    c.near("GSP fixed lambda == e", lam, E, 1e-4)
    c.cond("GSP n participants == 62", int(round(n_part)) == 62, f"nPart={n_part:g}")
    c.near("GSP trough time (s)", trough, 26.3, 1.5)
    # Control: the same S(x) form must NOT fit the rest blocks well.
    c.cond("CONTROL rest-block R^2 is low",
           rest_r2 < 0.40 and rest_r2 < grand_r2 - 0.3,
           f"rest R^2={rest_r2:.3f} << task R^2={grand_r2:.3f}")

    # --- Per-channel fits ---------------------------------------------------
    fits = pd.read_csv(os.path.join(results, "Figure2_perchannel_fits.csv"))
    c.cond("Per-channel fit count == 64", len(fits) == 64, f"{len(fits)} channels")
    c.near("Per-channel median R^2_free", fits["R2_free"].median(), 0.670, 0.02)

    # --- PAS 4/3 crossover --------------------------------------------------
    xo = pd.read_csv(os.path.join(results, "Figure3_crossover_stats.csv")).iloc[0]
    c.near("PAS crossover median (s)", xo["within_part_median_s"], 29.0, 3.0)
    c.cond("PAS early-vs-late significant", xo["early_late_p"] < 0.05,
           f"p={xo['early_late_p']:.4g}")
    c.cond("PAS effect direction (d_z > 0)", xo["early_late_dz"] > 0.15,
           f"d_z={xo['early_late_dz']:.3f}")
    c.cond("PAS within-trial slope significant", xo["trial_beta1_p"] < 1e-3,
           f"beta1 p={xo['trial_beta1_p']:.2g}")

    # --- Aperiodic 1/f exponent S(x) fit ------------------------------------
    with open(os.path.join(results, "aperiodic_exponent_within_trial.json")) as fh:
        ap_j = json.load(fh)
    free = ap_j["free_lambda_fit"]
    boot = ap_j["bootstrap_lambda_ci"]
    curv = ap_j["per_participant"]["curvature_t_test"]
    c.near("Aperiodic free-fit R^2", free["R2"], 0.876, 0.02)
    c.cond("Aperiodic peak in bootstrap CI", 18.0 <= free["peak_s"] <= 38.0,
           f"peak={free['peak_s']:.1f}s")
    c.cond("Aperiodic lambda CI brackets e", boot.get("brackets_e") is True,
           f"CI={boot['ci95']}")
    c.cond("Aperiodic population inverted-U", curv["t"] < -2.0 and curv["p"] < 0.01,
           f"t={curv['t']:.2f}, p={curv['p']:.2g}")

    ok = c.report()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
