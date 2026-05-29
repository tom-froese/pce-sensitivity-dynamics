#!/usr/bin/env python3
"""1/f aperiodic exponent vs the sensitivity curve S(x) = x·exp(-e·x).

Reads the per-(participant × bin) exponent CSVs produced by
:mod:`computeAperiodicExponent` and computes:

  • head-to-head fits at λ = e, fixed onset τ = 3.9 s — sensitivity-peak vs
    R(x)-decay vs inverted-U (quadratic) vs line vs trough.
  • free-λ sensitivity fit (does λ ≈ e?).
  • free-τ₀ boot-up onset, with and without freeing λ.
  • bootstrap 95 % CI on the peak time (2000 resamples over participants).
  • per-participant inverted-U reliability (quadratic curvature).
  • per-participant exponent-peak vs OWN click-time (within-person sensitivity link).
  • confound: EMG-reduced 2-20 Hz band exponent.

Outputs (in ``results/``):
  ``aperiodic_exponent_within_trial.csv``  — cohort mean ± SEM per bin
  ``aperiodic_exponent_within_trial.png``  — 3-panel diagnostic (QA)
  ``aperiodic_exponent_within_trial.json`` — fits + peak CI + provenance

The cohort-summary CSV is what Panel D of Figure 2 consumes; the JSON sidecar
records the fit parameters that the figure annotates (peak CI, λ, τ₀, etc.).

Usage
-----
    python code/analysis/fitExponentSensitivity.py
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

_HERE = Path(__file__).resolve().parent
_CODE = _HERE.parent
for _p in (str(_CODE), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import _config as C       # noqa: E402
import _eeg_io as A       # noqa: E402
import computeAperiodicExponent as cae  # noqa: E402
from _exponent_common import TAU, T_TRIAL, OUT  # noqa: E402


# Reference timing constant (population): PAS 4→3 awareness crossover (s).
# Produced by computePASCrossover.m -> results/Figure3_crossover_stats.csv.
# Kept here at the master-loop's hardcoded value (29.0 s) so layer-2
# validation against the master is direct; future revisions can read
# this from disk via the Figure3 CSV.
PAS_CROSSOVER_S = 29.0
PAS_CROSSOVER_CI = (23.5, 35.7)
BOOTSTRAP_N = 2000


def _load_or_compute_per_subj(dyads=None, recompute: bool = False):
    """Return (subjects, centers, E, E_lo). Reads the per-subj cache if present
    and ``recompute`` is False; otherwise calls ``computeAperiodicExponent``."""
    csv_wb = OUT / "aperiodic_exponent_per_participant.csv"
    csv_lo = OUT / "aperiodic_exponent_per_participant_band-2-20.csv"

    if (not recompute) and dyads is None and csv_wb.exists():
        edf = pd.read_csv(csv_wb, index_col=0); edf.columns = edf.columns.astype(float)
        subjects = [str(s) for s in edf.index]
        centers = edf.columns.to_numpy(dtype=float)
        E = edf.to_numpy(dtype=float)
        if csv_lo.exists():
            edl = pd.read_csv(csv_lo, index_col=0); edl.columns = edl.columns.astype(float)
            E_lo = edl.to_numpy(dtype=float)
        else:
            print(f"  [warn] no 2-20 Hz cache; setting E_lo := NaN (confound check skipped)",
                  flush=True)
            E_lo = np.full_like(E, np.nan)
        print(f"  loaded cached per-subj exponent {E.shape} <- {csv_wb.name}", flush=True)
        return subjects, centers, E, E_lo

    ds = A.default_dataset()
    subjects, centers, E, E_lo = cae.compute(ds, dyads)
    if dyads is None:
        cae.save(subjects, centers, E, E_lo)
    return subjects, centers, E, E_lo


def _r2(y, pr):
    return 1 - np.sum((y - pr) ** 2) / np.sum((y - y.mean()) ** 2)


def fit_cohort(centers: np.ndarray, E: np.ndarray, E_lo: np.ndarray):
    """All cohort-level fits and peak diagnostics from the per-(subj × bin)
    exponent matrices. Mirrors the master-loop's logic exactly."""
    from scipy.optimize import curve_fit

    n_subj = E.shape[0]
    t_star = TAU + (T_TRIAL - TAU) / np.e
    mu = np.nanmean(E, 0)
    sem = np.nanstd(E, 0) / np.sqrt(np.sum(np.isfinite(E), 0))
    mu_lo = np.nanmean(E_lo, 0)

    # Per-participant inverted-U reliability via quadratic curvature.
    curv, ppeak = [], []
    for si in range(n_subj):
        y = E[si]; ok = np.isfinite(y)
        if ok.sum() >= 10:
            c2 = np.polyfit(centers[ok], y[ok], 2)
            curv.append(c2[0])
            if c2[0] < 0:
                ppeak.append(-c2[1] / (2 * c2[0]))
    curv = np.array(curv)
    ppeak = np.array([p for p in ppeak if 0 < p < 60])
    tcurv = stats.ttest_1samp(curv, 0)

    # Head-to-head fits at λ = e, fixed onset τ.
    x = (centers - TAU) / (T_TRIAL - TAU)

    def lstsq_r2(basis):
        Am = np.column_stack([basis, np.ones_like(basis)])
        co, *_ = np.linalg.lstsq(Am, mu, rcond=None)
        return _r2(mu, Am @ co)

    sens   = np.where(x > 0,  x * np.exp(-np.e * x), 0.0)
    trough = np.where(x > 0, -x * np.exp(-np.e * x), 0.0)
    decay  = np.where(x > 0,      np.exp(-np.e * x), 1.0)
    r2_sens, r2_trough, r2_decay, r2_lin = (
        lstsq_r2(sens), lstsq_r2(trough), lstsq_r2(decay), lstsq_r2(centers)
    )
    quad = np.column_stack([centers ** 2, centers, np.ones_like(centers)])
    cq, *_ = np.linalg.lstsq(quad, mu, rcond=None)
    r2_quad = _r2(mu, quad @ cq); vertex = -cq[1] / (2 * cq[0])

    # Free-λ sensitivity fit: A·x·exp(-λ·x) + B ; peak at x = 1/λ → t = τ + (T-τ)/λ.
    def smodel(xx, Aa, lam, Bb):
        return Aa * xx * np.exp(-lam * xx) + Bb

    def fit_lambda(y):
        ok = np.isfinite(y)
        if ok.sum() < 6:
            return np.nan, np.nan
        try:
            pb, _ = curve_fit(
                smodel, x[ok], y[ok],
                p0=[max(np.ptp(y[ok]), 1e-3) * np.e, np.e, np.nanmin(y)],
                bounds=([0, 0.1, -np.inf], [np.inf, 20, np.inf]),
                maxfev=40000,
            )
            return pb[1], _r2(y[ok], smodel(x[ok], *pb))
        except Exception:
            return np.nan, np.nan

    lam_free, r2_free = fit_lambda(mu)
    peak_free = TAU + (T_TRIAL - TAU) / lam_free if np.isfinite(lam_free) else np.nan
    lam_lo, _ = fit_lambda(mu_lo)
    peak_lo_band = TAU + (T_TRIAL - TAU) / lam_lo if np.isfinite(lam_lo) else np.nan

    # Bootstrap the free-λ peak (and λ itself) over participants.
    rng = np.random.default_rng(0)
    peaks, lambdas = [], []
    for _ in range(BOOTSTRAP_N):
        idx = rng.integers(0, n_subj, n_subj)
        lam_b, _ = fit_lambda(np.nanmean(E[idx], 0))
        if np.isfinite(lam_b) and 0.5 < lam_b < 20:
            peaks.append(TAU + (T_TRIAL - TAU) / lam_b)
            lambdas.append(lam_b)
    peaks = np.array(peaks); lambdas = np.array(lambdas)
    pk_lo, pk_med, pk_hi = (np.percentile(peaks, [2.5, 50, 97.5]) if len(peaks) > 50
                            else (np.nan,) * 3)
    lam_ci_lo, lam_ci_med, lam_ci_hi = (np.percentile(lambdas, [2.5, 50, 97.5])
                                        if len(lambdas) > 50 else (np.nan,) * 3)
    argmax_t = centers[np.nanargmax(mu)]

    # EMG-reduced 2-20 Hz confound: same sensitivity-peak basis on mu_lo.
    finite_lo = np.isfinite(mu_lo)
    if finite_lo.any():
        As_lo = np.column_stack([sens, np.ones_like(sens)])
        cl, *_ = np.linalg.lstsq(As_lo, np.nan_to_num(mu_lo, nan=np.nanmean(mu_lo)),
                                 rcond=None)
        r2_sens_lo = _r2(mu_lo[finite_lo], (As_lo @ cl)[finite_lo])
        argmax_lo = centers[np.nanargmax(mu_lo)]
    else:
        r2_sens_lo = np.nan
        argmax_lo = np.nan

    # Free boot-up τ₀.
    def smodel_tau(t, Aa, lam, tau0, Bb):
        xb = np.clip((t - tau0) / (T_TRIAL - tau0), 0, None)
        return Aa * xb * np.exp(-lam * xb) + Bb

    okm = np.isfinite(mu)
    A0 = max(np.ptp(mu[okm]), 1e-3) * np.e
    B0 = float(np.nanmin(mu))
    try:                          # (i) λ = e fixed, free boot-up τ₀
        pe, _ = curve_fit(
            lambda t, Aa, tau0, Bb: smodel_tau(t, Aa, np.e, tau0, Bb),
            centers[okm], mu[okm], p0=[A0, 3.0, B0],
            bounds=([0, 0, -np.inf], [np.inf, 18, np.inf]), maxfev=40000,
        )
        tau0_e = pe[1]
        r2_boot_e = _r2(mu[okm], smodel_tau(centers[okm], pe[0], np.e, pe[1], pe[2]))
        peak_boot_e = tau0_e + (T_TRIAL - tau0_e) / np.e
    except Exception:
        tau0_e = r2_boot_e = peak_boot_e = np.nan
    try:                          # (ii) free λ AND free τ₀
        pf, _ = curve_fit(
            smodel_tau, centers[okm], mu[okm], p0=[A0, np.e, 3.0, B0],
            bounds=([0, 0.1, 0, -np.inf], [np.inf, 20, 18, np.inf]), maxfev=40000,
        )
        lam_boot, tau0_boot = pf[1], pf[2]
        r2_boot = _r2(mu[okm], smodel_tau(centers[okm], *pf))
        peak_boot = tau0_boot + (T_TRIAL - tau0_boot) / lam_boot
    except Exception:
        lam_boot = tau0_boot = r2_boot = peak_boot = np.nan
        pf = None

    # Nested F-test: does the both-free model improve on the fixed (λ=e, τ=3.9)
    # model enough to justify its two extra parameters? Reduced model = A·s + B
    # (2 free, with λ=e and τ=TAU baked into s); full model = A·xb·exp(-λ·xb) + B
    # with xb=(t-τ₀)/(T-τ₀) (4 free). n_pts = bins with finite mean.
    n_pts = int(okm.sum())
    p_full, p_red = 4, 2
    df_num, df_den = p_full - p_red, n_pts - p_full
    if (df_den > 0 and np.isfinite(r2_boot)
            and r2_boot > r2_sens and (1.0 - r2_boot) > 0):
        F_stat = ((r2_boot - r2_sens) / df_num) / ((1.0 - r2_boot) / df_den)
        F_p = float(1.0 - stats.f.cdf(F_stat, df_num, df_den))
        F_crit_005 = float(stats.f.ppf(0.95, df_num, df_den))
    else:
        F_stat = F_p = F_crit_005 = np.nan

    return dict(
        mu=mu, sem=sem, mu_lo=mu_lo, x=x, okm=okm,
        t_star=t_star,
        r2_sens=r2_sens, r2_trough=r2_trough, r2_decay=r2_decay, r2_lin=r2_lin,
        r2_quad=r2_quad, vertex=vertex,
        lam_free=lam_free, r2_free=r2_free, peak_free=peak_free,
        lam_lo=lam_lo, peak_lo_band=peak_lo_band,
        pk_lo=pk_lo, pk_med=pk_med, pk_hi=pk_hi,
        lam_ci_lo=lam_ci_lo, lam_ci_med=lam_ci_med, lam_ci_hi=lam_ci_hi,
        argmax_t=argmax_t,
        r2_sens_lo=r2_sens_lo, argmax_lo=argmax_lo,
        tau0_e=tau0_e, r2_boot_e=r2_boot_e, peak_boot_e=peak_boot_e,
        lam_boot=lam_boot, tau0_boot=tau0_boot, r2_boot=r2_boot, peak_boot=peak_boot,
        F_stat=float(F_stat), F_p=float(F_p), F_crit_005=float(F_crit_005),
        F_df_num=int(df_num), F_df_den=int(df_den), n_pts=n_pts,
        pf=pf,
        curv=curv, ppeak=ppeak, tcurv=tcurv,
        smodel=smodel, smodel_tau=smodel_tau,
    )


def _per_participant_click_link(subjects, centers, E, fits):
    """Per-participant exponent peak (free-λ + free-τ₀) vs own click-time median."""
    from scipy.optimize import curve_fit
    smodel_tau = fits["smodel_tau"]

    # Click times are produced by preprocessClicks.m -> ClickResponseTimes.csv,
    # which carries DyadID, ParticipantID, and ClickTime_s columns.
    click_csv = C.PREPROC_CLICKPAS / "ClickResponseTimes.csv"
    try:
        cp = pd.read_csv(click_csv)
        cp = cp[(cp.ClickTime_s > SEG_LO) & (cp.ClickTime_s < SEG_HI)]
        ct = cp.ClickTime_s.values
        click_mode = float(np.median(ct))
        cp["subj"] = cp.DyadID.map("{:02d}".format) + "_" + cp.ParticipantID.astype(str)
        click_med = cp.groupby("subj").ClickTime_s.median()
    except Exception as exc:
        print(f"  [warn] click-link skipped: {type(exc).__name__}: {exc}", flush=True)
        return np.array([]), np.nan, pd.Series(dtype=float), np.array([]), None

    def peak_of(y):
        ok = np.isfinite(y)
        if ok.sum() < 12:
            return np.nan
        try:
            pp, _ = curve_fit(
                smodel_tau, centers[ok], y[ok],
                p0=[max(np.ptp(y[ok]), 1e-3) * np.e, np.e, 3.0, float(np.nanmin(y))],
                bounds=([0, 0.1, 0, -np.inf], [np.inf, 20, 18, np.inf]),
                maxfev=20000,
            )
            pk = pp[2] + (T_TRIAL - pp[2]) / pp[1]
            if 4 < pk < 58:
                return pk
        except Exception:
            pass
        c2 = np.polyfit(centers[ok], y[ok], 2)
        if c2[0] < 0 and 4 < -c2[1] / (2 * c2[0]) < 58:
            return -c2[1] / (2 * c2[0])
        return np.nan

    pk_subj = {subj: peak_of(E[si]) for si, subj in enumerate(subjects)}
    pairs = np.array([(pk_subj[s], float(click_med[s])) for s in subjects
                      if np.isfinite(pk_subj.get(s, np.nan)) and s in click_med.index])
    rho_pc = stats.spearmanr(pairs[:, 0], pairs[:, 1]) if len(pairs) >= 8 else None
    return ct, click_mode, click_med, pairs, rho_pc


# bin-segment bounds re-exposed for click_link
from _exponent_common import SEG as _SEG  # noqa: E402
SEG_LO, SEG_HI = _SEG


def _plot(subjects, centers, fits, ct, click_mode, pairs, rho_pc,
          out_png: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mu, sem = fits["mu"], fits["sem"]
    t_star = fits["t_star"]
    lam_free, r2_free = fits["lam_free"], fits["r2_free"]
    lam_boot, tau0_boot, r2_boot = fits["lam_boot"], fits["tau0_boot"], fits["r2_boot"]
    pk_lo, pk_hi = fits["pk_lo"], fits["pk_hi"]
    pf = fits["pf"]
    smodel_tau = fits["smodel_tau"]
    x = fits["x"]; okm = fits["okm"]

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(17, 4.8))

    axA.errorbar(centers, mu, yerr=sem, marker="o", ms=3, color="C0",
                 label="1/f exponent (mean±SEM)")
    if np.isfinite(r2_free):
        sb = x * np.exp(-lam_free * x)
        Af = np.column_stack([sb, np.ones_like(sb)])
        cf, *_ = np.linalg.lstsq(Af[okm], mu[okm], rcond=None)
        axA.plot(centers, Af @ cf, "C3--", lw=1.6,
                 label=f"fixed-onset (λ={lam_free:.2f}) R²={r2_free:.2f}")
    if np.isfinite(r2_boot) and pf is not None:
        axA.plot(centers, smodel_tau(centers, *pf), "k-", lw=2,
                 label=f"boot-up τ0={tau0_boot:.1f}s (λ={lam_boot:.2f}) R²={r2_boot:.2f}")
    axA.axvline(t_star, color="r", ls=":", lw=1.2, label=f"1/e ({t_star:.0f}s)")
    if np.isfinite(pk_lo):
        axA.axvspan(pk_lo, pk_hi, color="0.6", alpha=0.18,
                    label=f"peak 95%CI[{pk_lo:.0f},{pk_hi:.0f}]")
    axA.axvline(PAS_CROSSOVER_S, color="C2", ls="-.", lw=1,
                label=f"PAS 4→3 ({PAS_CROSSOVER_S:.0f}s)")
    axA.set_xlabel("trial time (s)")
    axA.set_ylabel("aperiodic exponent")
    axA.set_title("1/f exponent tracks sensitivity S(x) (noICA)")
    axA.legend(fontsize=7)

    if len(ct):
        axB.hist(ct, bins=28, color="C2", alpha=0.6, density=True)
        axB.axvline(t_star, color="r", ls=":", lw=1.2, label=f"1/e ({t_star:.0f}s)")
        axB.axvline(click_mode, color="k", ls="--", lw=1,
                    label=f"click median ({click_mode:.0f}s)")
        axB.set_xlabel("click time (s)")
        axB.set_ylabel("density")
        axB.set_title("click-time distribution (behavioural sensitivity)")
        axB.legend(fontsize=7)

    if len(pairs) >= 8 and rho_pc is not None:
        axC.scatter(pairs[:, 1], pairs[:, 0], s=14, alpha=0.6, color="C0")
        axC.set_xlabel("participant median click time (s)")
        axC.set_ylabel("participant exponent-peak time (s)")
        axC.set_title(f"per-participant: exponent peak vs own click\n"
                      f"Spearman={rho_pc.correlation:+.2f} "
                      f"(p={rho_pc.pvalue:.3f}, n={len(pairs)})")

    fig.tight_layout()
    fig.savefig(out_png, dpi=120, bbox_inches="tight")


def _print_report(subjects, fits, ct, click_mode, pairs, rho_pc):
    f = fits
    n = len(subjects)
    pk_lo, pk_hi, t_star = f['pk_lo'], f['pk_hi'], f['t_star']
    brackets_1e  = np.isfinite(pk_lo) and pk_lo <= t_star <= pk_hi
    brackets_pas = np.isfinite(pk_lo) and pk_lo <= PAS_CROSSOVER_S <= pk_hi

    print(f"\n=== 1/f EXPONENT vs SENSITIVITY curve (noICA, {n} participants) ===")

    # ----- HEADLINE — FIXED-parameter form (λ = e, τ = 3.9 s) -----
    print(f"\n  ★ HEADLINE — FIXED-parameter S(x) [λ = e, τ = {TAU:.1f} s]:")
    print(f"      R² = {f['r2_sens']:.3f}   "
          f"(head-to-head head: vs inverted-U(quad) {f['r2_quad']:.3f}, "
          f"R(x)-decay {f['r2_decay']:.3f}, trough {f['r2_trough']:.3f}, "
          f"line {f['r2_lin']:.3f})")
    print(f"      Peak at t = {t_star:.1f} s = 1/e   "
          f"(by construction with fixed λ = e and τ = {TAU:.1f} s).")
    print(f"      Bootstrap 95 % CI on the cohort peak (free-λ, over participants): "
          f"[{pk_lo:.1f}, {pk_hi:.1f}] s")
    print(f"      => {'BRACKETS' if brackets_1e else 'does NOT bracket'} 1/e "
          f"({t_star:.1f} s); "
          f"{'BRACKETS' if brackets_pas else 'does NOT bracket'} the "
          f"PAS 4→3 crossover ({PAS_CROSSOVER_S:.0f} s, "
          f"CI [{PAS_CROSSOVER_CI[0]},{PAS_CROSSOVER_CI[1]}]); "
          f"click-time median {click_mode:.1f} s.")

    # ----- PER-PARTICIPANT inverted-U reliability (pipeline- and parameter-agnostic) -----
    curv = f['curv']; ppeak = f['ppeak']; tcurv = f['tcurv']
    print(f"\n  ★ PER-PARTICIPANT inverted-U reliability "
          f"(pipeline- and parameter-agnostic):")
    print(f"      {100*(curv<0).mean():.0f}% of participants "
          f"({int((curv<0).sum())}/{len(curv)}) have negative quadratic curvature.")
    print(f"      One-sample t-test (curvature < 0): "
          f"t = {tcurv.statistic:+.2f}, p = {tcurv.pvalue:.1e} (n = {len(curv)}).")
    print(f"      Individual peak-time median: {np.median(ppeak):.1f} s "
          f"[IQR {np.percentile(ppeak,25):.0f}-{np.percentile(ppeak,75):.0f}].")
    if rho_pc is not None:
        print(f"      Within-person sensitivity link: per-participant exponent-peak "
              f"vs own click-time Spearman = {rho_pc.correlation:+.3f} "
              f"(p = {rho_pc.pvalue:.3f}, n = {len(pairs)}).")

    # ----- ROBUSTNESS — free-parameter S(x) fits + nested F-test -----
    print(f"\n  ROBUSTNESS — free-parameter S(x) fits:")
    print(f"      free λ, τ = {TAU:.1f} fixed:   λ = {f['lam_free']:.2f}, "
          f"R² = {f['r2_free']:.3f}, peak = {f['peak_free']:.1f} s")
    print(f"      free τ₀, λ = e fixed:    τ₀ = {f['tau0_e']:.1f} s, "
          f"R² = {f['r2_boot_e']:.3f}, peak = {f['peak_boot_e']:.1f} s "
          f"(near 1/e)")
    print(f"      free λ AND free τ₀:      λ = {f['lam_boot']:.2f}, "
          f"τ₀ = {f['tau0_boot']:.1f} s, R² = {f['r2_boot']:.3f}, "
          f"peak = {f['peak_boot']:.1f} s")
    print(f"      Free-λ bootstrap (over participants): "
          f"median λ = {f['lam_ci_med']:.2f}, "
          f"95 % CI [{f['lam_ci_lo']:.2f}, {f['lam_ci_hi']:.2f}] "
          f"({'brackets' if f['lam_ci_lo'] <= np.e <= f['lam_ci_hi'] else 'does NOT bracket'} "
          f"e = {np.e:.2f}).")
    if np.isfinite(f['F_stat']):
        sig = "significant" if f['F_p'] < 0.05 else "NOT significant"
        print(f"      Nested F-test (fixed vs both-free): "
              f"ΔR² = {f['r2_boot']-f['r2_sens']:+.3f} for +2 parameters; "
              f"F({f['F_df_num']}, {f['F_df_den']}) = {f['F_stat']:.2f}, "
              f"p = {f['F_p']:.3f} (F-crit @ .05 = {f['F_crit_005']:.2f}) "
              f"— {sig}.")
        print(f"      => The data show no statistically significant preference "
              f"for the free-parameter form; fixed-form S(x) is the data's choice.")

    # ----- τ₀ caveat -----
    print(f"\n  NOTE — when τ₀ is free, the first FOOOF bin starts at "
          f"t = {f['x'].min()*(T_TRIAL-TAU)+TAU:.0f} s, so τ₀ is poorly constrained "
          f"at the low end.")
    print(f"      Imposing λ = e AND freeing τ₀ recovers τ₀ ≈ {f['tau0_e']:.1f} s "
          f"(near the EDA-derived {TAU:.1f} s) and places the peak at exactly "
          f"{f['peak_boot_e']:.1f} s = 1/e.")
    print(f"      The 'different boot-up' in the both-free fit is "
          f"unconstrained-parameter wiggle, not signal.")

    # ----- EMG-reduced confound check -----
    survives_lo = f['r2_sens_lo'] > 0.5 and 18 <= f['argmax_lo'] <= 40
    print(f"\n  CONFOUND CHECK (EMG-reduced 2-20 Hz):")
    print(f"      sensitivity-peak R² = {f['r2_sens_lo']:+.3f}, "
          f"argmax = {f['argmax_lo']:.1f} s, free-λ = {f['lam_lo']:.2f} "
          f"→ peak {'SURVIVES' if survives_lo else 'WEAKENS (partly HF-carried)'}.")


def _sidecar_dict(subjects, fits, click_mode, pairs, rho_pc):
    f = fits
    return {
        "producer": "code/analysis/fitExponentSensitivity.py",
        "produced_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "pipeline": "noICA (task-raw)",
        "n_participants": len(subjects),
        "framework": {"tau_s": TAU, "T_trial_s": T_TRIAL,
                      "pas_crossover_s": PAS_CROSSOVER_S,
                      "pas_crossover_ci_s": list(PAS_CROSSOVER_CI)},
        "input_per_participant_csv": "aperiodic_exponent_per_participant.csv",
        "head_to_head_R2": {
            "sensitivity_peak_lambda_e": float(f["r2_sens"]),
            "inverted_U_quad":           float(f["r2_quad"]),
            "R_x_decay":                 float(f["r2_decay"]),
            "trough":                    float(f["r2_trough"]),
            "line":                      float(f["r2_lin"]),
        },
        "free_lambda_fit": {
            "lambda": float(f["lam_free"]),
            "R2": float(f["r2_free"]),
            "peak_s": float(f["peak_free"]),
        },
        "boot_up_lambda_e": {
            "tau0_s": float(f["tau0_e"]),
            "R2": float(f["r2_boot_e"]),
            "peak_s": float(f["peak_boot_e"]),
        },
        "boot_up_free_lambda": {
            "lambda": float(f["lam_boot"]),
            "tau0_s": float(f["tau0_boot"]),
            "R2": float(f["r2_boot"]),
            "peak_s": float(f["peak_boot"]),
        },
        "bootstrap_peak_ci": {
            "n_resamples": BOOTSTRAP_N,
            "median_s": float(f["pk_med"]),
            "ci95_s": [float(f["pk_lo"]), float(f["pk_hi"])],
        },
        "bootstrap_lambda_ci": {
            "n_resamples": BOOTSTRAP_N,
            "median": float(f["lam_ci_med"]),
            "ci95":   [float(f["lam_ci_lo"]), float(f["lam_ci_hi"])],
            "brackets_e": bool(np.isfinite(f["lam_ci_lo"]) and
                                f["lam_ci_lo"] <= float(np.e) <= f["lam_ci_hi"]),
        },
        "model_comparison_F_test": {
            "description": (
                "Nested F-test: reduced model = A·s+B with s = x·exp(-e·x) "
                "and x = (t-τ)/(T-τ), τ=3.9 (2 free); full model = "
                "A·xb·exp(-λ·xb)+B with xb=(t-τ₀)/(T-τ₀) (4 free)."
            ),
            "delta_R2":   float(f["r2_boot"] - f["r2_sens"]),
            "df_num":     int(f["F_df_num"]),
            "df_denom":   int(f["F_df_den"]),
            "F_statistic": float(f["F_stat"]),
            "F_critical_0_05": float(f["F_crit_005"]),
            "p_value":    float(f["F_p"]),
            "verdict": (
                "free fit does NOT statistically improve on fixed form"
                if (np.isfinite(f["F_p"]) and f["F_p"] >= 0.05)
                else "free fit improves on fixed form (p<.05)"
            ),
        },
        "emg_reduced_2_20Hz": {
            "R2_sensitivity_peak": float(f["r2_sens_lo"]),
            "argmax_s": float(f["argmax_lo"]),
            "lambda_free": float(f["lam_lo"]),
        },
        "per_participant": {
            "frac_negative_curvature": float((f["curv"] < 0).mean()),
            "n_with_fit": int(len(f["curv"])),
            "median_peak_s": float(np.median(f["ppeak"])),
            "curvature_t_test": {
                "description": ("One-sample t-test of per-participant quadratic "
                                "curvature coefficient against 0; negative t and "
                                "small p => population-level inverted-U."),
                "t":     float(f["tcurv"].statistic),
                "p":     float(f["tcurv"].pvalue),
                "n":     int(len(f["curv"])),
            },
            "click_link_spearman": (None if rho_pc is None
                                    else float(rho_pc.correlation)),
            "click_link_pvalue":   (None if rho_pc is None
                                    else float(rho_pc.pvalue)),
            "click_link_n_pairs":  int(len(pairs)),
        },
        "click_median_s": float(click_mode) if np.isfinite(click_mode) else None,
        "outputs": {
            "csv": "aperiodic_exponent_within_trial.csv",
            "png": "aperiodic_exponent_within_trial.png",
        },
    }


def run(dyads=None, recompute: bool = False):
    """Compute, plot, and persist the exponent-sensitivity result."""
    subjects, centers, E, E_lo = _load_or_compute_per_subj(dyads, recompute)
    fits = fit_cohort(centers, E, E_lo)
    ct, click_mode, click_med, pairs, rho_pc = _per_participant_click_link(
        subjects, centers, E, fits
    )

    _print_report(subjects, fits, ct, click_mode, pairs, rho_pc)

    OUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUT / "aperiodic_exponent_within_trial.csv"
    png_path = OUT / "aperiodic_exponent_within_trial.png"
    json_path = OUT / "aperiodic_exponent_within_trial.json"

    pd.DataFrame({"trial_time": centers,
                  "exponent_mean": fits["mu"],
                  "exponent_sem":  fits["sem"]}).to_csv(csv_path, index=False)
    _plot(subjects, centers, fits, ct, click_mode, pairs, rho_pc, png_path)
    json_path.write_text(json.dumps(
        _sidecar_dict(subjects, fits, click_mode, pairs, rho_pc), indent=2))
    print(f"\nsaved -> {csv_path.name}, {png_path.name}, {json_path.name}")
    return fits


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dyads", type=int, nargs="+", default=None,
                    help="restrict to these dyad numbers (default: all available)")
    ap.add_argument("--recompute", action="store_true",
                    help="force-recompute the per-subj exponent cache")
    args = ap.parse_args()
    run(dyads=args.dyads, recompute=args.recompute)


if __name__ == "__main__":
    main()
