"""
exploreHemisphereFreeTau.py
=========================================================================
Diagnostic: fit S(x) = A*x*exp(-e*x) + B with tau free (per hemisphere)
and compare against the locked-tau (3.90 s) version used in Figure 2 row B.

The minimum of S lives at x* = 1/e regardless of A, B, so
  t_trough = tau + (1/e) * (T - tau),
which means the "shared 24.5 s trough" under the locked-tau model is
purely a model artefact.  When tau is allowed to vary, each hemisphere's
trough moves with its own optimal tau.

OUTPUTS:
  results/Figure2_mid_freeTauVariant.png
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar

# ---------------------------------------------------------------------------
# Constants / paths (matching plotFigure2_Neural.py)
# ---------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
HEM  = HERE / '_scratch' / 'parietal_hemisphere_data.npz'
OUT  = ROOT / 'results' / 'Figure2_mid_freeTauVariant.png'

E            = np.e
T_TRIAL      = 60.0
TAU_LOCKED   = 3.90
SMOOTH_WIN_S = 5.0
FS           = 10

# ---------------------------------------------------------------------------
# Load + smooth grand-medians
# ---------------------------------------------------------------------------
d = np.load(HEM)
tTask = d['tTask']
left, right = d['left_TC'], d['right_TC']


def smooth_median(mat):
    med = np.median(mat, axis=0)
    win = int(round(SMOOTH_WIN_S * FS))
    halfK = win // 2
    sm = np.convolve(med, np.ones(win) / win, mode='valid')
    t_sm = tTask[halfK: halfK + sm.size]
    return t_sm, sm


t_sm, sm_L = smooth_median(left)
_,    sm_R = smooth_median(right)


# ---------------------------------------------------------------------------
# Fitters
# ---------------------------------------------------------------------------
def _fit_given_tau(t_sm, sm, tau):
    """Linear A,B fit for a given tau; returns yhat, SSE, and diagnostics."""
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
    t_trough = tau + (1.0 / E) * Teff
    return dict(tau=tau, A=A, B=B, R2=r2, SSE=ss_res, t_trough=t_trough,
                yhat=yhat, t_fit=t_fit)


def fit_locked(t_sm, sm, tau=TAU_LOCKED):
    return _fit_given_tau(t_sm, sm, tau)


def fit_free_tau(t_sm, sm, bounds=(0.5, 20.0)):
    def obj(tau):
        return _fit_given_tau(t_sm, sm, tau)['SSE']
    res = minimize_scalar(obj, bounds=bounds, method='bounded',
                          options={'xatol': 1e-4})
    return _fit_given_tau(t_sm, sm, float(res.x))


# ---------------------------------------------------------------------------
# Fit both hemispheres both ways
# ---------------------------------------------------------------------------
L_lock = fit_locked(t_sm, sm_L)
R_lock = fit_locked(t_sm, sm_R)
L_free = fit_free_tau(t_sm, sm_L)
R_free = fit_free_tau(t_sm, sm_R)

print('Hemisphere |   tau (s)   A (uV)   B (uV)    R^2    t_trough (s)')
print('-' * 68)
for lbl, f in [('Left  lock ', L_lock), ('Left  free ', L_free),
               ('Right lock ', R_lock), ('Right free ', R_free)]:
    print(f'{lbl} | {f["tau"]:7.2f}  {f["A"]:+7.2f}  {f["B"]:+6.2f}  '
          f'{f["R2"]:6.3f}   {f["t_trough"]:6.2f}')

# ---------------------------------------------------------------------------
# Plot — locked vs. free tau, same panel style as Figure 2 row B
# ---------------------------------------------------------------------------
plt.rcParams.update({'font.size': 10,
                     'axes.spines.top': False,
                     'axes.spines.right': False})

col_L, col_R = '#1f77b4', '#d62728'

fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.2), sharex=True, sharey=True)

for ax, variant, label in [(axes[0], 'locked', rf'Locked $\tau = {TAU_LOCKED:.2f}$ s'),
                           (axes[1], 'free',   r'Free $\tau$ (per hemisphere)')]:

    # Raw grand-median traces (faded)
    ax.plot(t_sm, sm_L, color=col_L, lw=1.0, alpha=0.35, label='Left parietal')
    ax.plot(t_sm, sm_R, color=col_R, lw=1.0, alpha=0.35, label='Right parietal')

    fL = L_lock if variant == 'locked' else L_free
    fR = R_lock if variant == 'locked' else R_free

    # Fit overlays
    ax.plot(fL['t_fit'], fL['yhat'], color=col_L, lw=2.2,
            label=(rf'L fit: $\tau={fL["tau"]:.2f}$ s,  $A={fL["A"]:+.1f}\,\mu$V,  '
                   rf'$R^2={fL["R2"]:.2f}$'))
    ax.plot(fR['t_fit'], fR['yhat'], color=col_R, lw=2.2,
            label=(rf'R fit: $\tau={fR["tau"]:.2f}$ s,  $A={fR["A"]:+.1f}\,\mu$V,  '
                   rf'$R^2={fR["R2"]:.2f}$'))

    # Trough markers for each fit
    trough_pts = []
    for f, col in [(fL, col_L), (fR, col_R)]:
        y_tr = f['A'] * (1.0 / E) * np.exp(-1.0) + f['B']
        ax.plot(f['t_trough'], y_tr, 'o', mfc=col, mec='white',
                mew=1.2, ms=7, zorder=5)
        trough_pts.append((f['t_trough'], y_tr))

    # Tau-lag verticals (dotted) for each fit
    for f, col in [(fL, col_L), (fR, col_R)]:
        ax.axvline(f['tau'], color=col, ls=':', lw=0.9, alpha=0.6)

    # Trough-timing callout (single label when locked; one per hemisphere when free)
    if variant == 'locked':
        t_tr = fL['t_trough']  # identical for L and R
        y_mid = 0.5 * (trough_pts[0][1] + trough_pts[1][1])
        ax.plot([t_tr, t_tr], [trough_pts[0][1], trough_pts[1][1]],
                ls='--', color='0.4', lw=0.9, zorder=4)
        ax.text(t_tr, y_mid,
                rf'$t_{{\mathrm{{trough}}}} = {t_tr:.1f}$ s',
                ha='center', va='center', fontsize=8.5,
                bbox=dict(boxstyle='round,pad=0.25', fc='white',
                          ec='0.7', lw=0.6),
                zorder=6)
    else:
        for (t_tr, y_tr), col in zip(trough_pts, [col_L, col_R]):
            ax.annotate(rf'$t_{{\mathrm{{trough}}}} = {t_tr:.1f}$ s',
                        xy=(t_tr, y_tr), xytext=(6, -10),
                        textcoords='offset points',
                        fontsize=8.5, color=col)

    ax.axhline(0, color='k', lw=0.4, alpha=0.3)
    ax.set_xlim(t_sm.min(), t_sm.max())
    ax.set_xlabel('Time (s)')
    ax.set_title(label, pad=6)
    ax.legend(loc='upper center', frameon=False, fontsize=8.3, ncol=1,
              handlelength=1.8)

axes[0].set_ylabel(r'ROI-mean scalp potential ($\mu$V)')

fig.suptitle(r'Parietal hemispheres — locked vs. free $\tau$ variants '
             r'of $S(x) = A\,x\,e^{-e\,x} + B$',
             y=1.00, fontsize=11, fontweight='bold')
fig.tight_layout()
fig.savefig(OUT, dpi=300, bbox_inches='tight')
print(f'\nSaved: {OUT}')
