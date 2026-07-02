#!/usr/bin/env python3
"""Build the per-click ClickPAS table (click time + PAS rating per click).

Merges the bundled click-response times
(``data/preprocessed/ClickTimes/ClickResponseTimes.csv``) with the per-trial
Perceptual Awareness Scale (PAS) ratings from the raw behavioral
questionnaires (``data/raw/Behavior/pceXX*/questionnaires/
pair_XX_PY_PAS_confidence_absence.csv``, the ``click_presence`` question), and
writes ``data/preprocessed/ClickPAS/ClickPAS.csv``.

This is the shared input for the two PAS-dependent statistics that the paper
reports beyond the crossover: the within-participant click-time/PAS Spearman
correlation (``computePASSpearman.py``) and the neural<->phenomenal decoupling
Bayes factors (``computeDecouplingBF.py``). It mirrors the ``merged`` table that
the MATLAB ``computePASCrossover.m`` builds internally, but persisted as a small
derived CSV so the Python analyses need only the raw *behavior* (not the raw
EEG) to regenerate.

Output columns: ``DyadID, ParticipantID, TrialNum, ClickTime_s, PAS`` (one row
per clicked trial with a PAS rating; dyad 31 excluded, matching the rest of the
pipeline). The PAS files use 0-indexed trial ids; click data is 1-indexed — the
merge reconciles this (``trial_id + 1``).

Usage
-----
    python code/analysis/computeClickPAS.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
_CODE = _HERE.parent
for _p in (str(_CODE), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import _config as C  # noqa: E402

_DYAD_RE = re.compile(r"^pce(\d{2})\d{6}$")


def build() -> pd.DataFrame:
    ct_path = C.PREPROC_CLICKPAS / "ClickResponseTimes.csv"
    if not ct_path.is_file():
        raise SystemExit(f"missing {ct_path} — run preprocessClicks.m first (see README)")
    ct = pd.read_csv(ct_path)
    # DyadID here is the raw folder id (e.g. 1230807); the dyad number is its
    # leading pair of digits (1..32).
    ct["DyadNum"] = (ct.DyadID // 1_000_000).astype(int)
    ct = ct[~ct.DyadNum.isin(C.EXCLUDE_DYADS)].copy()

    raw_beh = C.DATA_RAW / "Behavior"
    if not raw_beh.is_dir():
        raise SystemExit(
            f"missing {raw_beh} — the raw behavioral questionnaires carry the PAS "
            "ratings; download data/raw/Behavior/ from OSF (see README)."
        )

    pas_rows = []
    for d in sorted(raw_beh.glob("pce*")):
        m = _DYAD_RE.match(d.name)
        if not m:
            continue
        dN = int(m.group(1))
        if dN in C.EXCLUDE_DYADS:
            continue
        for p in C.PARTICIPANTS:
            fp = d / "questionnaires" / f"pair_{dN:02d}_P{p}_PAS_confidence_absence.csv"
            if not fp.is_file():
                continue
            tbl = pd.read_csv(fp)
            sub = tbl[tbl.question_id == "click_presence"]
            for _, r in sub.iterrows():
                pas_rows.append((dN, p, int(r.trial_id) + 1, int(r.answer)))
    pas = pd.DataFrame(pas_rows, columns=["DyadNum", "ParticipantID", "TrialNum", "PAS"])

    ct["Key"] = ct.DyadNum * 1000 + ct.ParticipantID * 100 + ct.TrialNum
    pas["Key"] = pas.DyadNum * 1000 + pas.ParticipantID * 100 + pas.TrialNum
    merged = ct.merge(pas[["Key", "PAS"]], on="Key", how="inner")
    merged = merged[(merged.Clicked == 1) & merged.ClickTime_s.notna()].copy()

    out = (merged[["DyadNum", "ParticipantID", "TrialNum", "ClickTime_s", "PAS"]]
           .rename(columns={"DyadNum": "DyadID"})
           .sort_values(["DyadID", "ParticipantID", "TrialNum"])
           .reset_index(drop=True))
    return out


def main():
    df = build()
    out_dir = C.DATA_PREPROC / "ClickPAS"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "ClickPAS.csv"
    df.to_csv(out, index=False)
    n_subj = (df.DyadID * 10 + df.ParticipantID).nunique()
    print(f"  built {len(df)} click-PAS rows across {n_subj} participants "
          f"-> {out.relative_to(C.REPO_ROOT)}")


if __name__ == "__main__":
    main()
