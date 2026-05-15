"""
plotFigure2_Neural.py
=========================================================================
PNAS Figure 2 — Neural evidence, rebuilt from scratch (April 2026).

Two-row figure:

  TOP ROW (A) — Scalp-map triptych on the strict-filter set
             (tau* in [2, 6] s AND R^2_locked >= 0.70).
             Panels:  (i)  tau*  (s)
                      (ii) R^2 at locked tau = 3.90 s
                      (iii) Signed trough magnitude at tau-locked fit:
                            A * s(1/e) = A * (1/e) * e^{-1}.
                            Plotted with a negative-valued colormap to
                            honour the fact that the minimum of S(x) at
                            x = 1/e is a trough (A is negative on the
                            canonical posterior-lateral belt).

  BOTTOM ROW (B) — ROI-mean scalp potential of the left and right
             parietal hemispheres (L: P1/P3/P5/P7; R: P2/P4/P6/P8),
             overlaid with locked-tau S(x) = A*x*exp(-e*x) + B fits
             (tau = 3.90 s, fixed from the full-parietal free-tau
             optimum).

NOTE: an earlier bottom row reporting a 0.05-0.2 Hz hemispheric rhythm
and its Hilbert relative phase was removed after per-participant
reproducibility tests (split-halves across participants, split-halves
across trials within participant, phase-randomised surrogates, and
inter-participant PLV) all failed to distinguish the apparent
"~0.1 Hz rhythm" from phase-scrambled per-participant noise.  See
code/analysis/_scratch/haptic_rhythm_test.py for the diagnostic.

INPUTS:
  code/analysis/_scratch/parietal_hemisphere_data.npz
      keys: tTask, left_TC, mid_TC, right_TC
  results/FigureS2_GSP_TopoMap_FreeTau_perchannel.csv
      per-channel free-tau and locked-tau S(x) fits
  data/preprocessed/EEG/globalScalpPotential_data.mat  (reference / 9-ch parietal)

OUTPUTS:
  results/Figure2_Neural.png  (300 dpi)
  results/Figure2_Neural.pdf

AUTHOR: Embodied Cognitive Science Unit, OIST
"""
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

# ---------------------------------------------------------------------------
# Paths — resolve relative to this file so the script is location-independent
# ---------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent                       # .../pce-sensitivity-dynamics/
HEM  = HERE / '_scratch' / 'parietal_hemisphere_data.npz'
GSP  = ROOT / 'data' / 'preprocessed' / 'EEG' / 'globalScalpPotential_data.mat'
CSV  = ROOT / 'results' / 'FigureS2_GSP_TopoMap_FreeTau_perchannel.csv'
OUT_DIR = ROOT / 'results'
OUT_PNG = OUT_DIR / 'Figure2_Neural.png'
OUT_PDF = OUT_DIR / 'Figure2_Neural.pdf'

# ---------------------------------------------------------------------------
# Constants (matching the joint analysis)
# ---------------------------------------------------------------------------
E            = np.e
T_TRIAL      = 60.0
TAU_LOCKED   = 3.90       # joint-tau lock (seconds)
SMOOTH_WIN_S = 5.0        # grand-median boxcar smoothing window (seconds)
FS           = 10         # sampling rate of the task time course (Hz)

# Strict-filter criteria for the top-row topomaps
TAU_LO, TAU_HI = 2.0, 6.0
R2_MIN         = 0.70

# ---------------------------------------------------------------------------
# Load hemisphere data
# ---------------------------------------------------------------------------
d = np.load(HEM)
tTask = d['tTask']
left  = d['left_TC']          # (nPart, nSamples)
right = d['right_TC']
mid   = d['mid_TC']
n_part = left.shape[0]
print(f'Loaded hemisphere data: n_part = {n_part}, samples = {left.shape[1]}')


def smooth_median(mat):
    """Boxcar-smoothed grand median of (nPart, nSamples)."""
    med = np.median(mat, axis=0)
    win = int(round(SMOOTH_WIN_S * FS))
    halfK = win // 2
    sm = np.convolve(med, np.ones(win) / win, mode='valid')
    t_sm = tTask[halfK: halfK + sm.size]
    return t_sm, sm


def fit_locked(t_sm, sm, tau=TAU_LOCKED):
    """Least-squares fit of S(x) = A*x*exp(-e*x) + B past `tau`."""
    Teff = T_TRIAL - tau
    m = t_sm >= tau
    t_fit = t_sm[m]
    y_fit = sm[m]
    x = (t_fit - tau) / Teff
    s = x * np.exp(-E * x)
    X = np.column_stack([s, np.ones_like(s)])
    coef, *_ = np.linalg.lstsq(X, y_fit, rcond=None)
    A, B = float(coef[0]), float(coef[1])
    yhat = A * s + B
    ss_res = float(np.sum((y_fit - yhat) ** 2))
    ss_tot = float(np.sum((y_fit - y_fit.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot
    trough_x = 1.0 / E
    t_trough = tau + trough_x * Teff
    return dict(A=A, B=B, R2=r2, t_trough=t_trough, yhat=yhat, t_fit=t_fit)


# Load full-parietal reference (for the trough landmark)
def _cellstr(f, ref):
    arr = np.array(f[ref]).squeeze()
    return ''.join(chr(int(c)) for c in arr.ravel())


with h5py.File(GSP, 'r') as f:
    name_refs = f['roi/names']
    rpc_refs  = f['roiParticipantTC']
    roi_names = [_cellstr(f, name_refs[i, 0]) for i in range(name_refs.shape[0])]
    par_idx = roi_names.index('Parietal')
    par_full = np.asarray(f[rpc_refs[0, par_idx]]).T  # (nPart, 600)

# Grand-median, smoothed time courses + locked-tau fits per hemisphere
t_sm, sm_left  = smooth_median(left)
_,    sm_right = smooth_median(right)
_,    sm_full  = smooth_median(par_full)
fit_full  = fit_locked(t_sm, sm_full,  TAU_LOCKED)
fit_left  = fit_locked(t_sm, sm_left,  TAU_LOCKED)
fit_right = fit_locked(t_sm, sm_right, TAU_LOCKED)
trough_full = fit_full['t_trough']
for name, fit in [('Full parietal', fit_full),
                  ('Left  parietal', fit_left),
                  ('Right parietal', fit_right)]:
    print(f'{name} (locked tau = {TAU_LOCKED:.2f} s):  '
          f'A = {fit["A"]:+.2f}, B = {fit["B"]:+.2f}, '
          f'R^2 = {fit["R2"]:.3f},  trough at t = {fit["t_trough"]:.2f} s')

# ---------------------------------------------------------------------------
# Scalp-map data (strict filter)
# ---------------------------------------------------------------------------
df = pd.read_csv(CSV)
mask_tau_ok = (df['tau_star'] >= TAU_LO) & (df['tau_star'] <= TAU_HI)
mask_r2_ok  = df['R2_lock'] >= R2_MIN
keep = mask_tau_ok & mask_r2_ok
kept = df[keep].reset_index(drop=True)

N_TOTAL      = len(df)
N_DROP_TAU   = int((~mask_tau_ok).sum())
N_AFTER_TAU  = N_TOTAL - N_DROP_TAU
N_DROP_R2    = int((mask_tau_ok & ~mask_r2_ok).sum())
N_KEPT       = int(keep.sum())
print(
    f'Strict filter:  start = {N_TOTAL},  '
    f'dropped by tau*-range  = {N_DROP_TAU}  -> {N_AFTER_TAU} remain,  '
    f'dropped by R^2 cutoff = {N_DROP_R2}  -> {N_KEPT} retained.'
)

ch_names = kept['channel'].tolist()
montage = mne.channels.make_standard_montage('standard_1020')
info = mne.create_info(ch_names=ch_names, sfreq=1000., ch_types='eeg')
info.set_montage(montage, match_case=False, on_missing='raise')

tau_vals = kept['tau_star'].to_numpy()
R2_lock  = kept['R2_lock'].to_numpy()
A_lock   = kept['A_lock'].to_numpy()

# Signed trough magnitude at locked tau: A(ch) * x_trough * exp(-e*x_trough)
# with x_trough = 1/e  ->  shape factor = (1/e) * e^{-1} = e^{-2}
trough_shape = (1.0 / E) * np.exp(-1.0)
trough_signed = A_lock * trough_shape          # uV, signed (negative on posterior belt)

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
plt.rcParams.update({
    'font.size': 10,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

fig = plt.figure(figsize=(9.0, 7.2))
gs = GridSpec(
    2, 3,
    height_ratios=[0.95, 1.05],
    hspace=0.55, wspace=0.32,
    left=0.08, right=0.96, top=0.92, bottom=0.09,
)

# ---- TOP ROW — three scalp maps ----
ax_topo1 = fig.add_subplot(gs[0, 0])
ax_topo2 = fig.add_subplot(gs[0, 1])
ax_topo3 = fig.add_subplot(gs[0, 2])

# Panel (i): tau* — first filter stage (restrict to the canonical lag range)
im0, _ = mne.viz.plot_topomap(
    tau_vals, info, axes=ax_topo1, show=False, cmap='coolwarm',
    vlim=(TAU_LO, TAU_HI), contours=6, sensors=True, names=ch_names,
)
ax_topo1.set_title(
    rf'$\tau^{{*}}$  (s)' + '\n'
    rf'filter: $\tau^{{*}} \in [{TAU_LO:.0f},\,{TAU_HI:.0f}]$ s  '
    rf'(drops {N_DROP_TAU}/{N_TOTAL} ch.)',
    pad=6, fontsize=9,
)
cb0 = fig.colorbar(im0, ax=ax_topo1, shrink=0.75, pad=0.05)
cb0.set_label(r'$\tau^{*}$ (s)')

# Panel (ii): R^2 at locked tau — second filter stage (quality cutoff at 3.9 s)
im1, _ = mne.viz.plot_topomap(
    R2_lock, info, axes=ax_topo2, show=False, cmap='viridis',
    vlim=(R2_MIN, float(np.nanmax(R2_lock))), contours=6, sensors=True,
    names=ch_names,
)
ax_topo2.set_title(
    rf'R$^2$ at locked $\tau$ = {TAU_LOCKED:.2f} s' + '\n'
    rf'filter: R$^2 \geq {R2_MIN:.2f}$  '
    rf'(drops {N_DROP_R2}/{N_AFTER_TAU} remaining)',
    pad=6, fontsize=9,
)
cb1 = fig.colorbar(im1, ax=ax_topo2, shrink=0.75, pad=0.05)
cb1.set_label(r'R$^2$')

# Panel (iii): signed trough magnitude — final set of retained channels
tm_max = float(np.nanmax(np.abs(trough_signed)))
im2, _ = mne.viz.plot_topomap(
    trough_signed, info, axes=ax_topo3, show=False, cmap='magma',
    vlim=(-tm_max, 0.0), contours=6, sensors=True, names=ch_names,
)
ax_topo3.set_title(
    r'Signed trough magnitude  $A \cdot s(1/e)$' + '\n'
    rf'({N_KEPT}/{N_TOTAL} channels retained; $\mu$V)',
    pad=6, fontsize=9,
)
cb2 = fig.colorbar(im2, ax=ax_topo3, shrink=0.75, pad=0.05, extend='min')
cb2.set_label(r'$\mu$V (trough, $<$0)')

# ---- BOTTOM ROW — parietal L/R ROI-mean scalp potential + locked-tau S(x) fits ----
ax_mid = fig.add_subplot(gs[1, :])
col_L, col_R = '#1f77b4', '#d62728'

# Raw grand-median time courses (de-emphasised)
ax_mid.plot(t_sm, sm_left,  color=col_L, lw=1.0, alpha=0.35,
            label='Left parietal  (P1/P3/P5/P7)')
ax_mid.plot(t_sm, sm_right, color=col_R, lw=1.0, alpha=0.35,
            label='Right parietal (P2/P4/P6/P8)')

# Locked-tau S(x) fits on top (only defined for t >= TAU_LOCKED)
ax_mid.plot(fit_left['t_fit'],  fit_left['yhat'],  color=col_L, lw=2.2,
            label=rf'L fit: $A={fit_left["A"]:+.1f}\ \mu$V,  $R^2={fit_left["R2"]:.2f}$')
ax_mid.plot(fit_right['t_fit'], fit_right['yhat'], color=col_R, lw=2.2,
            label=rf'R fit: $A={fit_right["A"]:+.1f}\ \mu$V,  $R^2={fit_right["R2"]:.2f}$')

# Per-hemisphere trough markers (minimum of each locked-tau fit).
# NOTE: under a locked tau the minimum of S(x) sits at x = 1/e irrespective
# of A and B, so t_trough is identical for L and R by construction.  We show
# the markers as the natural visual anchors of each fit but do NOT annotate a
# shared "t_trough" callout — that would misleadingly suggest an empirically
# discovered temporal alignment rather than a model constraint.
for name, fit, col in [('L', fit_left, col_L), ('R', fit_right, col_R)]:
    y_trough = fit['A'] * (1.0 / E) * np.exp(-1.0) + fit['B']
    ax_mid.plot(fit['t_trough'], y_trough, 'o',
                mfc=col, mec='white', mew=1.2, ms=7, zorder=5)

# Locked-tau landmark (same convention as Figure 1)
ax_mid.axvline(TAU_LOCKED, color='0.5', ls=':', lw=1.0)
ax_mid.axhline(0, color='k', lw=0.4, alpha=0.3)

ax_mid.set_xlim(t_sm.min(), t_sm.max())
ax_mid.set_xlabel('Time (s)')
ax_mid.set_ylabel(r'ROI-mean scalp potential ($\mu$V)')
ax_mid.set_title(rf'Parietal hemispheres — grand median with '
                 rf'$S(x) = A\,x\,e^{{-e\,x}} + B$ fits  '
                 rf'($\tau = {TAU_LOCKED:.2f}$ s, fixed from full-parietal '
                 rf'free-$\tau$ optimum)', pad=6, fontsize=10)
ax_mid.legend(loc='upper center', frameon=False, fontsize=8.5, ncol=2,
              columnspacing=1.4, handlelength=1.8)

# Suptitle
fig.suptitle(
    rf'Figure 2 — Neural evidence  (N = {n_part} participants)',
    y=0.98, fontsize=12, fontweight='bold',
)

# Panel letters
for ax, letter in [(ax_topo1, 'A'), (ax_mid, 'B')]:
    ax.text(-0.08, 1.08, letter, transform=ax.transAxes,
            fontsize=16, fontweight='bold', va='bottom', ha='left')

fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight')
fig.savefig(OUT_PDF, bbox_inches='tight')
print(f'Saved: {OUT_PNG}')
print(f'Saved: {OUT_PDF}')
