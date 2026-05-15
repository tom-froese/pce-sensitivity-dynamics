"""
analyseHapticEEGCoupling.py
=========================================================================
Systematic per-participant test of whether the haptic residual
(Figure 1 Panel B) is coupled to the bandpassed parietal-hemisphere
signal (Figure 2 Panel C), and whether that coupling is regime-specific
(pre- vs. post-sensitivity-trough).

Layers of analysis:
  (1) REPLICATION
      For each participant matched across haptic + EEG:
        - compute smoothed haptic mean time course (across their trials)
        - subtract the GRAND exp-rise fit to form the haptic residual
        - bandpass-filter their own L and R parietal time courses (0.05-0.20 Hz)
        - compute per-participant zero-lag Pearson r(haptic_resid, R-hem)
          and r(haptic_resid, L-hem) in three windows:
            * full stable  (10-50 s)
            * pre-trough   (t95 - t_trough)
            * post-trough  (t_trough - 50 s)
        - compute per-participant xcorr peak lag in each window.
        - compute per-participant L-R PLV in each window.

  (2) SURROGATE CONTROLS
      Shuffle surrogate: for each participant i, pair haptic_i with
      R-hem_j (j != i) and recompute r.  Build a null distribution of
      "mismatched" r values and compare against the matched distribution.
      This test is specific to within-subject coupling and controls for
      filtering/spectral artefacts (both legs of the null share spectra
      with the real data).

  (3) REGIME TRANSITION
      Compare per-participant r, lag, and PLV between pre- and post-
      trough windows with paired tests (Wilcoxon signed-rank).

Outputs:
  results/Figure2_HapticEEGCoupling_summary.png
  results/Figure2_HapticEEGCoupling_perParticipant.csv

AUTHOR: Embodied Cognitive Science Unit, OIST
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.signal import firwin, filtfilt, hilbert, correlate
from scipy.stats import pearsonr, wilcoxon, ttest_1samp

# ---------------------------------------------------------------------------
# Paths + constants
# ---------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
EEG_NPZ = HERE / '_scratch' / 'parietal_hemisphere_data.npz'
HAPTIC_CSV = ROOT / 'data' / 'preprocessed' / 'Haptics' / 'HapticFeedback.csv'
OUT_FIG = ROOT / 'results' / 'Figure2_HapticEEGCoupling_summary.png'
OUT_CSV = ROOT / 'results' / 'Figure2_HapticEEGCoupling_perParticipant.csv'

FS_EEG       = 10            # EEG sampled at 10 Hz after block-avg
FS_HAPTIC    = 100           # 100 Hz
SMOOTH_EEG_S = 5.0           # boxcar smoothing (already used in Figure 2)
SMOOTH_HAP_S = 5.0           # same kernel, keep signals comparable
BP_LO, BP_HI = 0.05, 0.20    # bandpass for EEG hemispheric signal
FIR_N        = 101
TAU_LOCKED   = 3.90
T_TRIAL      = 60.0
E            = np.e
T_TROUGH     = TAU_LOCKED + (1.0 / E) * (T_TRIAL - TAU_LOCKED)   # 24.54 s
POST_END     = 50.0

RNG = np.random.default_rng(42)

# ---------------------------------------------------------------------------
# Load EEG hemisphere data (already per-participant grand-median per trial)
# ---------------------------------------------------------------------------
d = np.load(EEG_NPZ)
tTask   = d['tTask']          # (600,)  time axis at 10 Hz
left    = d['left_TC']        # (nPart, 600)  per-participant mean across trials
right   = d['right_TC']
dyad_id = d['dyad_id']
pid     = d['pid']
n_part  = left.shape[0]
print(f'EEG participants: {n_part}')

# Smooth each participant's trace (same kernel as grand-median in Figure 2)
def smooth_ts(x, fs, win_s):
    w = int(round(win_s * fs))
    kernel = np.ones(w) / w
    out = np.array([np.convolve(x[i], kernel, mode='valid') for i in range(x.shape[0])])
    hK = w // 2
    t_out = tTask[hK: hK + out.shape[1]]  # only valid for EEG here
    return t_out, out

t_eeg_sm, left_sm  = smooth_ts(left,  FS_EEG, SMOOTH_EEG_S)
_,        right_sm = smooth_ts(right, FS_EEG, SMOOTH_EEG_S)

# Zero-phase FIR bandpass per participant
taps = firwin(FIR_N, [BP_LO, BP_HI], pass_zero=False, fs=FS_EEG, window='hamming')
bp_L = filtfilt(taps, 1.0, left_sm,  axis=1, padlen=3*FIR_N)
bp_R = filtfilt(taps, 1.0, right_sm, axis=1, padlen=3*FIR_N)
print(f'Bandpassed per-participant L and R: shape = {bp_L.shape}')

# ---------------------------------------------------------------------------
# Load haptic data and build per-participant mean trace
# ---------------------------------------------------------------------------
print('Loading haptic CSV (~7M rows)...')
haptic_df = pd.read_csv(HAPTIC_CSV)
# All trials share the same Time_s grid (5999 samples @ 100 Hz, 0.01-59.98 s)
# Average across trials WITHIN a participant -> per-participant haptic trace
print('Building per-participant haptic traces...')
pp_mean = haptic_df.groupby(
    ['DyadID', 'ParticipantID', 'Time_s'])['HapticFeedback'].mean().reset_index()
# Pivot to (DyadID, ParticipantID) x Time_s
pivot = pp_mean.pivot_table(
    index=['DyadID', 'ParticipantID'], columns='Time_s',
    values='HapticFeedback', aggfunc='mean')
t_haptic = pivot.columns.to_numpy()
haptic_pp_full = pivot.to_numpy()                # (64, 5999)
hkey_all = np.array(list(pivot.index))           # (64, 2)

# Smooth with the same 5-s boxcar
w = int(round(SMOOTH_HAP_S * FS_HAPTIC))
kernel = np.ones(w) / w
haptic_sm = np.array([np.convolve(haptic_pp_full[i], kernel, mode='valid')
                      for i in range(haptic_pp_full.shape[0])])
hK = w // 2
t_haptic_sm = t_haptic[hK: hK + haptic_sm.shape[1]]

# Match haptic participants to EEG participants by (dyad, pid).
# Haptic DyadID packs dyad_number * 1_000_000 + date_code (e.g. 1230807 = dyad 1).
# EEG dyad_id is just the dyad number 1-32.
def haptic_dyad_to_num(dyad_id_long):
    return int(dyad_id_long) // 1_000_000

eeg_keys = {(int(dyad_id[i]), int(pid[i])): i for i in range(n_part)}
matched_eeg_idx = []
matched_hap_idx = []
for j, (D, P) in enumerate(hkey_all):
    key = (haptic_dyad_to_num(D), int(P))
    if key in eeg_keys:
        matched_hap_idx.append(j)
        matched_eeg_idx.append(eeg_keys[key])
matched_hap_idx = np.asarray(matched_hap_idx, dtype=int)
matched_eeg_idx = np.asarray(matched_eeg_idx, dtype=int)
print(f'Matched haptic-EEG participants: {len(matched_hap_idx)}')

haptic_sm_m = haptic_sm[matched_hap_idx]         # (nMatch, nHapticSm)
bp_L_m      = bp_L[matched_eeg_idx]
bp_R_m      = bp_R[matched_eeg_idx]
n_match     = len(matched_hap_idx)

# ---------------------------------------------------------------------------
# Grand exp-rise fit on haptic grand mean (this is the Figure 1B model)
# ---------------------------------------------------------------------------
grand_haptic = haptic_sm_m.mean(axis=0)

def exp_rise(t, A, tau_h, onset):
    x = np.maximum((t - onset) / tau_h, 0.0)
    return A * (1 - np.exp(-x))

popt, _ = curve_fit(exp_rise, t_haptic_sm, grand_haptic, p0=[0.9, 2.0, 1.0], maxfev=20000)
A_h, tau_h, onset_h = popt
t95 = onset_h - tau_h * np.log(1 - 0.95)
print(f'Grand haptic fit: A={A_h:.3f}, tau={tau_h:.2f}, onset={onset_h:.2f}, t95={t95:.2f} s')

# Per-participant haptic residual (their smoothed - grand fit)
grand_fit = exp_rise(t_haptic_sm, *popt)
haptic_resid_pp = haptic_sm_m - grand_fit[None, :]   # (nMatch, nHapticSm)

# Interpolate each residual onto the EEG time grid
def interp_to(t_src, y, t_dst):
    f = interp1d(t_src, y, axis=-1, bounds_error=False, fill_value=np.nan)
    return f(t_dst)

hap_on_eeg = interp_to(t_haptic_sm, haptic_resid_pp, t_eeg_sm)   # (nMatch, nEegSm)

# ---------------------------------------------------------------------------
# Per-participant windowed statistics
# ---------------------------------------------------------------------------
def window_mask(t, lo, hi):
    return (t >= lo) & (t <= hi)

m_pre  = window_mask(t_eeg_sm, max(t95, 4.0), T_TROUGH)
m_post = window_mask(t_eeg_sm, T_TROUGH, POST_END)
m_full = window_mask(t_eeg_sm, 10.0,       POST_END)

def corr_xcorr(a, b, fs=FS_EEG, max_lag_s=10.0):
    """Zero-lag Pearson r, plus peak r / lag within +-max_lag_s."""
    a = a - a.mean(); b = b - b.mean()
    r0 = np.corrcoef(a, b)[0, 1]
    xc = correlate(a, b, mode='full')
    xc /= (np.std(a) * np.std(b) * len(a))
    lags = np.arange(-len(a) + 1, len(a)) / fs
    sel = np.abs(lags) <= max_lag_s
    xc_s, lags_s = xc[sel], lags[sel]
    pk = np.argmax(np.abs(xc_s))
    return r0, xc_s[pk], lags_s[pk]

def plv(bL, bR):
    phi_L = np.angle(hilbert(bL))
    phi_R = np.angle(hilbert(bR))
    dphi = np.angle(np.exp(1j * (phi_R - phi_L)))
    return np.abs(np.mean(np.exp(1j * dphi))), np.angle(np.mean(np.exp(1j * dphi)))

rows = []
for i in range(n_match):
    hL = hap_on_eeg[i]
    bL = bp_L_m[i]
    bR = bp_R_m[i]
    rec = {'dyad': haptic_dyad_to_num(hkey_all[matched_hap_idx[i]][0]),
           'pid' : int(hkey_all[matched_hap_idx[i]][1])}
    for wname, mask in [('full', m_full), ('pre', m_pre), ('post', m_post)]:
        # Drop any NaN samples (can happen at haptic interp edges)
        valid = mask & ~np.isnan(hL)
        if valid.sum() < 20:
            for col in ['r_hR', 'r_hL', 'xc_peak', 'xc_lag', 'PLV', 'dphi_mu']:
                rec[f'{col}_{wname}'] = np.nan
            continue
        a, bL_w, bR_w = hL[valid], bL[valid], bR[valid]
        r_hR, xc_pk, xc_lag = corr_xcorr(a, bR_w)
        r_hL, _, _ = corr_xcorr(a, bL_w)
        P, dphi_mu = plv(bL_w, bR_w)
        rec[f'r_hR_{wname}']    = r_hR
        rec[f'r_hL_{wname}']    = r_hL
        rec[f'xc_peak_{wname}'] = xc_pk
        rec[f'xc_lag_{wname}']  = xc_lag
        rec[f'PLV_{wname}']     = P
        rec[f'dphi_mu_{wname}'] = dphi_mu
    rows.append(rec)

res = pd.DataFrame(rows)
res.to_csv(OUT_CSV, index=False)
print(f'\nPer-participant results saved: {OUT_CSV}')

# ---------------------------------------------------------------------------
# Summaries + statistical tests
# ---------------------------------------------------------------------------
def summarise(col):
    x = res[col].dropna().to_numpy()
    if len(x) < 2:
        return dict(n=len(x))
    med = np.median(x); q25, q75 = np.percentile(x, [25, 75])
    t_stat, t_p = ttest_1samp(x, 0.0)
    try:
        w_stat, w_p = wilcoxon(x)
    except ValueError:
        w_stat, w_p = np.nan, np.nan
    frac_pos = (x > 0).mean()
    return dict(n=len(x), mean=x.mean(), std=x.std(), median=med,
                q25=q25, q75=q75, frac_pos=frac_pos,
                t=t_stat, t_p=t_p, w=w_stat, w_p=w_p)

print('\n' + '='*78)
print('LAYER 1 — per-participant correlations (haptic residual vs. parietal bandpass)')
print('='*78)
for w in ['pre', 'post', 'full']:
    for side in ['R', 'L']:
        col = f'r_h{side}_{w}'
        s = summarise(col)
        print(f"  {col:<14s} n={s['n']:2d}  median={s['median']:+.3f}  "
              f"mean={s['mean']:+.3f} +/- {s['std']:.3f}   "
              f"frac>0={s['frac_pos']:.2f}   t-test p={s['t_p']:.2e}   "
              f"signed-rank p={s['w_p']:.2e}")

print('\n' + '='*78)
print('LAYER 1 — xcorr peak lags (negative = R-hem lags haptic residual)')
print('='*78)
for w in ['pre', 'post', 'full']:
    col = f'xc_lag_{w}'
    s = summarise(col)
    print(f"  {col:<14s} median={s['median']:+.2f} s   "
          f"IQR=[{s['q25']:+.2f}, {s['q75']:+.2f}]   "
          f"signed-rank vs 0 p={s['w_p']:.2e}")

print('\n' + '='*78)
print('LAYER 1 — L-R phase locking (in-band)')
print('='*78)
for w in ['pre', 'post', 'full']:
    col = f'PLV_{w}'
    s = summarise(col)
    print(f"  {col:<14s} median={s['median']:.3f}  "
          f"IQR=[{s['q25']:.3f}, {s['q75']:.3f}]")

# LAYER 2 — SURROGATES ----------------------------------------------------
# Pair haptic_i with R-hem_j (j != i) and recompute r.  Build null
# distribution from all n*(n-1) mismatches per window.
print('\n' + '='*78)
print('LAYER 2 — mismatched-participant surrogate null '
      '(haptic_i vs R_j, j != i)')
print('='*78)
for w, mask in [('pre', m_pre), ('post', m_post), ('full', m_full)]:
    matched_rs, null_rs = [], []
    for i in range(n_match):
        hLi = hap_on_eeg[i]
        valid = mask & ~np.isnan(hLi)
        if valid.sum() < 20:
            continue
        a = hLi[valid] - hLi[valid].mean()
        # matched
        bR_i = bp_R_m[i][valid]
        bR_i = bR_i - bR_i.mean()
        if np.std(a) > 0 and np.std(bR_i) > 0:
            matched_rs.append(np.corrcoef(a, bR_i)[0, 1])
        # mismatched
        for j in range(n_match):
            if j == i:
                continue
            bR_j = bp_R_m[j][valid] - bp_R_m[j][valid].mean()
            if np.std(bR_j) > 0:
                null_rs.append(np.corrcoef(a, bR_j)[0, 1])
    matched_rs = np.array(matched_rs)
    null_rs = np.array(null_rs)
    m_mean, m_med = matched_rs.mean(), np.median(matched_rs)
    n_mean, n_med = null_rs.mean(), np.median(null_rs)
    # One-sided empirical p-value: probability that a null r >= matched median
    p_emp = (null_rs >= m_med).mean()
    # Alternative: does matched distribution shift vs null distribution?
    from scipy.stats import mannwhitneyu
    u, p_u = mannwhitneyu(matched_rs, null_rs, alternative='greater')
    print(f'  [{w:4s}]  matched  n={len(matched_rs)} mean={m_mean:+.3f} median={m_med:+.3f}')
    print(f'           null     n={len(null_rs)} mean={n_mean:+.3f} median={n_med:+.3f}')
    print(f'           Mann-Whitney U (matched > null) p = {p_u:.2e}   '
          f'empirical p(null>=matched median) = {p_emp:.3f}')
    # stash for plotting
    res.attrs.setdefault('surrogate', {})[w] = dict(matched=matched_rs, null=null_rs)

# LAYER 3 — regime-transition paired tests --------------------------------
print('\n' + '='*78)
print('LAYER 3 — paired pre vs post trough')
print('='*78)
for stem in ['r_hR', 'r_hL', 'xc_lag', 'PLV']:
    pre  = res[f'{stem}_pre'].to_numpy()
    post = res[f'{stem}_post'].to_numpy()
    valid = ~(np.isnan(pre) | np.isnan(post))
    try:
        w_stat, w_p = wilcoxon(pre[valid], post[valid])
    except ValueError:
        w_stat, w_p = np.nan, np.nan
    print(f'  {stem:<8s} pre median={np.median(pre[valid]):+.3f}   '
          f'post median={np.median(post[valid]):+.3f}   '
          f'signed-rank p={w_p:.2e}')

# ---------------------------------------------------------------------------
# Summary figure
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(13, 7.0))

# (a) per-subject r, R-hem, across windows
ax = axes[0, 0]
data = [res['r_hR_pre'].dropna(), res['r_hR_post'].dropna(), res['r_hR_full'].dropna()]
labels = ['pre', 'post', 'full']
ax.boxplot(data, tick_labels=labels, showmeans=True, showfliers=True)
ax.axhline(0, color='0.4', lw=0.8, ls='--')
ax.set_ylabel('Per-subject r (haptic residual, R-hem bp)')
ax.set_title('(a) Haptic → R-hem coupling by window')

# (b) same but L-hem
ax = axes[0, 1]
data = [res['r_hL_pre'].dropna(), res['r_hL_post'].dropna(), res['r_hL_full'].dropna()]
ax.boxplot(data, tick_labels=labels, showmeans=True, showfliers=True)
ax.axhline(0, color='0.4', lw=0.8, ls='--')
ax.set_ylabel('Per-subject r (haptic residual, L-hem bp)')
ax.set_title('(b) Haptic → L-hem coupling by window')

# (c) matched vs null surrogate (full window)
ax = axes[0, 2]
sd = res.attrs['surrogate']['pre']
ax.hist(sd['null'], bins=40, density=True, alpha=0.5, label='null (mismatched)',
        color='0.6')
ax.hist(sd['matched'], bins=20, density=True, alpha=0.7, label='matched',
        color='#d62728')
ax.axvline(0, color='0.4', lw=0.8, ls='--')
ax.set_xlabel('r(haptic resid, R-hem bp)')
ax.set_ylabel('density')
ax.set_title('(c) Pre-trough: matched vs. mismatched surrogates')
ax.legend(frameon=False, fontsize=9)

# (d) xcorr peak lags
ax = axes[1, 0]
data = [res['xc_lag_pre'].dropna(), res['xc_lag_post'].dropna(),
        res['xc_lag_full'].dropna()]
ax.boxplot(data, tick_labels=labels, showmeans=True, showfliers=True)
ax.axhline(0, color='0.4', lw=0.8, ls='--')
ax.set_ylabel('xcorr peak lag (s)\n(neg = R-hem lags haptic)')
ax.set_title('(d) Lead–lag by window')

# (e) L-R PLV per window with per-subject pairing
ax = axes[1, 1]
pre, post = res['PLV_pre'].to_numpy(), res['PLV_post'].to_numpy()
for i in range(n_match):
    ax.plot([0, 1], [pre[i], post[i]], '-', color='0.7', lw=0.6, alpha=0.6)
ax.boxplot([pre[~np.isnan(pre)], post[~np.isnan(post)]],
           positions=[0, 1], tick_labels=['pre', 'post'], widths=0.35,
           showmeans=True)
ax.set_ylabel('L-R PLV (in-band)')
ax.set_title('(e) Hemispheric phase locking pre vs post trough')

# (f) surrogate in post-trough
ax = axes[1, 2]
sd = res.attrs['surrogate']['post']
ax.hist(sd['null'], bins=40, density=True, alpha=0.5, label='null (mismatched)',
        color='0.6')
ax.hist(sd['matched'], bins=20, density=True, alpha=0.7, label='matched',
        color='#1f77b4')
ax.axvline(0, color='0.4', lw=0.8, ls='--')
ax.set_xlabel('r(haptic resid, R-hem bp)')
ax.set_ylabel('density')
ax.set_title('(f) Post-trough: matched vs. mismatched surrogates')
ax.legend(frameon=False, fontsize=9)

fig.suptitle(f'Per-participant haptic residual ↔ parietal bandpass coupling  '
             f'(N = {n_match})', y=1.00, fontweight='bold')
fig.tight_layout()
fig.savefig(OUT_FIG, dpi=300, bbox_inches='tight')
print(f'\nSaved figure: {OUT_FIG}')
