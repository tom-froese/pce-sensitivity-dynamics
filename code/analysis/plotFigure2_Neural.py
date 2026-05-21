"""
plotFigure2_Neural.py
=========================================================================
PNAS Figure 2 — Neural Evidence

Two-row, three-column figure:

  TOP LEFT (A) — Sensitivity schematic: R(x) and dR/dlambda ∝ S(x),
             with rejection/selection phases and trough at x = 1/e.

  TOP RIGHT (B) — Global scalp potential (64-channel mean) grand median,
             overlaid with free-tau S(x) = A*x*exp(-e*x) + B fit.

  BOTTOM ROW (C) — Scalp-map triptych (all 64 channels, uncorrected):
             (i)   R^2_free — where does the sensitivity model fit?
             (ii)  tau* — best-fit boot-up lag per channel
             (iii) Signed trough value — effect size and polarity

INPUTS:
  data/preprocessed/EEG/allchannel_data.mat  (scipy format)
  results/FigureS2_GSP_TopoMap_FreeTau_perchannel.csv

OUTPUTS:
  results/Figure2_Neural.pdf  (vector)
  results/Figure2_Neural.png  (300 dpi raster for inspection)

AUTHOR: Embodied Cognitive Science Unit, OIST
DATE:   May 2026
=========================================================================
"""
from pathlib import Path

import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import scipy.io as sio

# -----------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------
HERE    = Path(__file__).resolve().parent
ROOT    = HERE.parent.parent
CSV     = ROOT / 'results' / 'FigureS2_GSP_TopoMap_FreeTau_perchannel.csv'
OUT_PNG = ROOT / 'results' / 'Figure2_Neural.png'
OUT_PDF = ROOT / 'results' / 'Figure2_Neural.pdf'

# -----------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------
E            = np.e
T_TRIAL      = 60.0
TAU_LOCKED   = 3.90
SMOOTH_WIN_S = 5.0
FS           = 10


# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------

def smooth_median(mat, tTask):
    med = np.median(mat, axis=0)
    win = int(round(SMOOTH_WIN_S * FS))
    halfK = win // 2
    sm = np.convolve(med, np.ones(win) / win, mode='valid')
    t_sm = tTask[halfK: halfK + sm.size]
    return t_sm, sm


def fit_locked(t_sm, sm, tau=TAU_LOCKED):
    Teff = T_TRIAL - tau
    m = t_sm >= tau
    t_fit, y_fit = t_sm[m], sm[m]
    x = (t_fit - tau) / Teff
    s = x * np.exp(-E * x)
    X = np.column_stack([s, np.ones_like(s)])
    coef, *_ = np.linalg.lstsq(X, y_fit, rcond=None)
    A, B = float(coef[0]), float(coef[1])
    yhat = A * s + B
    ss_res = float(np.sum((y_fit - yhat) ** 2))
    ss_tot = float(np.sum((y_fit - y_fit.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot
    t_trough = tau + (1.0 / E) * Teff
    return dict(A=A, B=B, R2=r2, t_trough=t_trough, yhat=yhat, t_fit=t_fit)


def make_info(ch_names):
    montage = mne.channels.make_standard_montage('standard_1020')
    info = mne.create_info(ch_names=ch_names, sfreq=1000., ch_types='eeg')
    info.set_montage(montage, match_case=False, on_missing='warn')
    return info


# -----------------------------------------------------------------------
# 1. Load all-channel data → global scalp potential
# -----------------------------------------------------------------------
ALLCH = ROOT / 'data' / 'preprocessed' / 'EEG' / 'allchannel_data.mat'
d_all  = sio.loadmat(str(ALLCH))
tTask  = d_all['tTask'].squeeze()
allTC  = d_all['allTC']                       # (n_part, 64, 600)
gsp    = allTC.mean(axis=1)                   # (n_part, 600)
n_part = gsp.shape[0]
print(f'Global scalp potential: n_part = {n_part}, samples = {gsp.shape[1]}')

t_sm, sm_gsp = smooth_median(gsp, tTask)

# Free-tau grid search for global SP
best_r2, best_tau = -np.inf, None
for tau_try in np.arange(0.0, 20.01, 0.10):
    Teff = T_TRIAL - tau_try
    if Teff <= 0:
        continue
    m = t_sm >= tau_try
    if m.sum() < 5:
        continue
    t_f, y_f = t_sm[m], sm_gsp[m]
    x = (t_f - tau_try) / Teff
    s = x * np.exp(-E * x)
    X = np.column_stack([s, np.ones_like(s)])
    coef, *_ = np.linalg.lstsq(X, y_f, rcond=None)
    yhat = float(coef[0]) * s + float(coef[1])
    ss_res = np.sum((y_f - yhat) ** 2)
    ss_tot = np.sum((y_f - y_f.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot
    if r2 > best_r2:
        best_r2, best_tau = r2, tau_try

TAU_GSP = best_tau
fit_gsp = fit_locked(t_sm, sm_gsp, tau=TAU_GSP)
print(f'Global SP (free tau = {TAU_GSP:.2f} s):  '
      f'A = {fit_gsp["A"]:+.2f}, B = {fit_gsp["B"]:+.2f}, '
      f'R² = {fit_gsp["R2"]:.3f},  trough at t = {fit_gsp["t_trough"]:.2f} s')

# -----------------------------------------------------------------------
# 2. Scalp-map data — all 64 channels, no filtering
# -----------------------------------------------------------------------
df = pd.read_csv(str(CSV))
N_TOTAL = len(df)

ch_all      = df['channel'].tolist()
R2_free_all = df['R2_free'].to_numpy()
tau_all     = df['tau_star'].to_numpy()
A_lock_all  = df['A_lock'].to_numpy()
R2_lock_all = df['R2_lock'].to_numpy()
info_all    = make_info(ch_all)

trough_shape  = (1.0 / E) * np.exp(-1.0)
trough_raw    = A_lock_all * trough_shape + df['B_lock'].to_numpy()

print(f'\nScalp-map data: all {N_TOTAL} channels')
print(f'  tau* range:   [{tau_all.min():.2f}, {tau_all.max():.2f}] s')
print(f'  R2_free range: [{R2_free_all.min():.3f}, {R2_free_all.max():.3f}]')
print(f'  R2_lock range: [{R2_lock_all.min():.3f}, {R2_lock_all.max():.3f}]')
print(f'  trough range:  [{trough_raw.min():.3f}, {trough_raw.max():.3f}] µV')

# -----------------------------------------------------------------------
# 3. Figure
# -----------------------------------------------------------------------
plt.rcParams.update({
    'font.size': 10,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

fig = plt.figure(figsize=(13.0, 12.0))
gs_outer = GridSpec(
    2, 1,
    height_ratios=[1.0, 1.3],
    hspace=0.175,
    left=0.06, right=0.96, top=0.90, bottom=0.05,
)
gs_top = GridSpecFromSubplotSpec(
    1, 2, subplot_spec=gs_outer[0], wspace=0.45,
    width_ratios=[1, 1],
)
gs_bot = GridSpecFromSubplotSpec(
    1, 3, subplot_spec=gs_outer[1], wspace=0.35,
)

# ---- TOP LEFT (A): Sensitivity schematic ----

ax_a = fig.add_subplot(gs_top[0, 0])
col_decay = '#808080'
col_sens  = '#d93320'
col_early = '#2673bf'
col_late  = '#d98c19'

x_dim = np.linspace(0, 1.0, 300)
Rx     = np.exp(-E * x_dim)
dRdlam = -x_dim * np.exp(-E * x_dim)
x_trough = 1.0 / E
y_trough_schematic = dRdlam[np.searchsorted(x_dim, x_trough)]

ax_a.fill_betweenx([-1, 0], 0, x_trough,
                    color=col_early, alpha=0.06)
ax_a.fill_betweenx([-1, 0], x_trough, 1.0,
                    color=col_late, alpha=0.06)

ax_a.plot(x_dim, Rx, '-', color=col_decay, lw=2.5, label=r'$R(x)$: reliability')
ax_a_r = ax_a.twinx()
ax_a_r.plot(x_dim, dRdlam, '-', color=col_sens, lw=2.5,
            label=r'$\partial R/\partial\lambda \propto S(x)$: rate sensitivity')
ax_a_r.plot(x_trough, y_trough_schematic, '^', ms=9,
            mfc=col_sens, mec='k', mew=0.5, zorder=5)
ax_a_r.annotate(
    rf'$x = 1/e$',
    xy=(x_trough, y_trough_schematic),
    xytext=(12, -8), textcoords='offset points',
    fontsize=8, fontweight='bold', color=col_sens,
)

ax_a.set_xlim(0, 1.0)
ax_a.set_ylim(-0.05, 1.05)
ax_a_r.set_ylim(-1, 0)
ax_a.set_xlabel(r'Dimensionless time  $x = k\,t\,/\,\lambda$', fontsize=9)
ax_a.set_ylabel(r'$R(x) = e^{-\lambda x}$', fontsize=9, color=col_decay)
ax_a_r.set_ylabel(r'$\partial R / \partial\lambda$', fontsize=9, color=col_sens)
ax_a.tick_params(axis='y', colors=col_decay)
ax_a_r.tick_params(axis='y', colors=col_sens)
ax_a_r.spines['right'].set_visible(True)
ax_a_r.spines['right'].set_color(col_sens)
ax_a.spines['left'].set_color(col_decay)

ax_a.text(0.08, 0.08, 'Rejection\nphase', transform=ax_a.transAxes,
          fontsize=8, color=col_early, fontstyle='italic')
ax_a.text(0.62, 0.08, 'Selection\nphase', transform=ax_a.transAxes,
          fontsize=8, color=col_late, fontstyle='italic')

lines_a = [l for l in ax_a.get_lines() + ax_a_r.get_lines()
           if not l.get_label().startswith('_')]
labels_a = [l.get_label() for l in lines_a]
ax_a.legend(lines_a, labels_a, loc='center right', frameon=True,
            fontsize=7.5, fancybox=False, edgecolor='0.7')

ax_a.set_title('Sensitivity of reliability to rate perturbations',
               pad=6, fontsize=10)

# ---- TOP RIGHT (B): Global scalp potential + S(x) fit ----

ax_b = fig.add_subplot(gs_top[0, 1])
col_gsp = '#c44e52'

ax_b.plot(t_sm, sm_gsp, color=col_gsp, lw=1.0, alpha=0.35,
          label='Global scalp potential (64-ch mean)')
ax_b.plot(fit_gsp['t_fit'], fit_gsp['yhat'], color=col_gsp, lw=2.5,
          label=rf'$\tau^*$={TAU_GSP:.1f} s,  '
                rf'A={fit_gsp["A"]:+.1f} $\mu$V,  '
                rf'$R^2$={fit_gsp["R2"]:.2f}')

y_trough_raw = fit_gsp['A'] * (1.0 / E) * np.exp(-1.0) + fit_gsp['B']
ax_b.plot(fit_gsp['t_trough'], y_trough_raw, 'o',
          mfc=col_gsp, mec='white', mew=1.2, ms=8, zorder=5)
ax_b.annotate(
    rf't = {fit_gsp["t_trough"]:.1f} s',
    xy=(fit_gsp['t_trough'], y_trough_raw),
    xytext=(0, 10), textcoords='offset points',
    ha='center', va='bottom', fontsize=9, fontweight='bold',
    color=col_gsp,
)

ax_b.axvline(TAU_GSP, color='0.5', ls=':', lw=1.0)
ax_b.axhline(0, color='k', lw=0.4, alpha=0.3)

ax_b.set_xlim(t_sm.min(), t_sm.max())
ax_b.set_xlabel('Time (s)')
ax_b.set_ylabel(r'Scalp potential ($\mu$V)')
ax_b.set_title(
    rf'Global scalp potential — $S(x)$ fit  '
    rf'(free $\tau^*$ = {TAU_GSP:.1f} s)',
    pad=6, fontsize=10,
)
ax_b.legend(loc='upper right', frameon=False, fontsize=7.5)

# ---- BOTTOM ROW (C): three scalp maps ----

ax_topo1 = fig.add_subplot(gs_bot[0, 0])
ax_topo2 = fig.add_subplot(gs_bot[0, 1])
ax_topo3 = fig.add_subplot(gs_bot[0, 2])

im1, _ = mne.viz.plot_topomap(
    R2_free_all, info_all, axes=ax_topo1, show=False, cmap='viridis',
    vlim=(0, 1), contours=6, sensors=True, extrapolate='local',
)
ax_topo1.set_title(r'$R^2_{\mathrm{free}}$  (best-fit quality)',
                   pad=6, fontsize=9)
cb1 = fig.colorbar(im1, ax=ax_topo1, shrink=0.75, pad=0.05)
cb1.set_label(r'$R^2$')

im2, _ = mne.viz.plot_topomap(
    tau_all, info_all, axes=ax_topo2, show=False, cmap='coolwarm',
    contours=6, sensors=True, extrapolate='local',
)
ax_topo2.set_title(r'$\tau^{*}$  (best-fit boot-up lag)',
                   pad=6, fontsize=9)
cb2 = fig.colorbar(im2, ax=ax_topo2, shrink=0.75, pad=0.05)
cb2.set_label(r'$\tau^{*}$ (s)')

tm_max = float(np.nanmax(np.abs(trough_raw)))
im3, _ = mne.viz.plot_topomap(
    trough_raw, info_all, axes=ax_topo3, show=False, cmap='RdBu',
    vlim=(-tm_max, tm_max), contours=6, sensors=True, extrapolate='local',
)
ax_topo3.set_title(r'Signed trough value  $A \cdot s(1/e) + B$',
                   pad=6, fontsize=9)
cb3 = fig.colorbar(im3, ax=ax_topo3, shrink=0.75, pad=0.05)
cb3.set_label(r'$\mu$V')

# ---- Panel letters and suptitle ----

fig.suptitle(
    rf'Figure 2 — Neural evidence  (N = {n_part} participants)',
    y=0.96, fontsize=12, fontweight='bold',
)

for ax, letter in [(ax_a, 'A'), (ax_b, 'B'), (ax_topo1, 'C')]:
    ax.text(-0.08, 1.08, letter, transform=ax.transAxes,
            fontsize=16, fontweight='bold', va='bottom', ha='left')

# ---- Save ----
fig.savefig(str(OUT_PNG), dpi=300, bbox_inches='tight', facecolor='w')
fig.savefig(str(OUT_PDF), bbox_inches='tight', facecolor='w')
print(f'Saved: {OUT_PNG}')
print(f'Saved: {OUT_PDF}')
plt.close(fig)
