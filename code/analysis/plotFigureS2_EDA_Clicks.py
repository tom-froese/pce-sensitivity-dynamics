"""
Figure S2 — EDA and click response-time distribution fit independently
to the rate-sensitivity shape  x * exp(-e * x).

Left column: per-participant grand-average task EDA on the trial horizon
(x = (t-tau)/(T-tau), rate = e) with R(x) = A0 * exp(-e*x) + B overlaid.
The signed rate-sensitivity curve  dR/dlambda|_{lambda=e} = -A0*x*exp(-e*x)
is shown on a secondary axis (green, trough at x = 1/e, i.e.
t = tau + (T-tau)/e).

Right column: distribution of click times for the augmented moving-target
cohort (avatar + shadow + rescued 'none'-clicks, N = 948), together with
the reflection-corrected KDE produced by Matlab's ksdensity. The shape
A*x*exp(-e*x) is fit to the KDE with tau as a free parameter; the
amplitude A is obtained by projection.

Bottom row: a cohort-sensitivity panel showing how the click-side tau
estimate depends on the definition of the moving-target cohort:
    strict       N = 870  (OriginalTarget in {avatar, shadow} only)
    augmented    N = 948  (strict + 78 'none'-clicks whose nearest prior
                           contact was with a moving target — canonical)
    reassigned   N = 884  (every click re-routed to its nearest prior
                           contact; moving-target-only subset)

Inputs:
    data/preprocessed/EDA/EDA_Task_Preprocessed.csv
    data/preprocessed/ClickTimes/ClickResponseTimes.csv
    code/analysis/_scratch/reassigned_clicks.csv  (for the rescued and
        reassigned cohorts; produced by reassign_none.py)
    results/click_kde_matlab.csv  (Matlab KDE for each cohort; produced
        once via the MATLAB MCP call documented in FigureS2_methods.md)

Output:
    results/FigureS2_EDA_Clicks_IndependentLags.png
    results/FigureS2_EDA_Clicks_IndependentLags.pdf

See code/analysis/FigureS2_methods.md for the full methodology and the
audit of the 84 'none'-classified clicks.
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path('/sessions/great-gallant-dirac/mnt/pce-sensitivity-dynamics')
CLICK = ROOT / 'data/preprocessed/ClickTimes/ClickResponseTimes.csv'
EDA   = ROOT / 'data/preprocessed/EDA/EDA_Task_Preprocessed.csv'
REAS  = ROOT / 'code/analysis/_scratch/reassigned_clicks.csv'
MKDE  = ROOT / 'results/click_kde_matlab.csv'
OUT_PNG = ROOT / 'results/FigureS2_EDA_Clicks_IndependentLags.png'
OUT_PDF = OUT_PNG.with_suffix('.pdf')

E = np.e
T_TRIAL = 60.0
FS = 25
N_SAMPLES = int(T_TRIAL * FS)
TAU_GRID = np.arange(0.0, 15.01, 0.1)

# =========================================================================
# EDA fit (same as before)
# =========================================================================
df_eda = pd.read_csv(EDA)
df_eda['PID'] = df_eda.DyadID * 10 + df_eda.ParticipantID
part_matrix = np.full((df_eda.PID.nunique(), N_SAMPLES), np.nan)
for i, pid in enumerate(np.sort(df_eda.PID.unique())):
    pd_ = df_eda[df_eda.PID == pid]
    trial_means = []
    for _, tg in pd_.groupby('TrialNum'):
        if len(tg) == N_SAMPLES:
            trial_means.append(tg.sort_values('Time_s').EDA_uS.to_numpy())
    if trial_means:
        part_matrix[i] = np.mean(np.stack(trial_means), axis=0)
part_matrix = part_matrix[~np.isnan(part_matrix).any(axis=1)]
n_part = part_matrix.shape[0]
time_eda = np.arange(N_SAMPLES) / FS
y = part_matrix.mean(axis=0)
sem = part_matrix.std(axis=0, ddof=0) / np.sqrt(n_part)

def R_matlab(x_norm, A0, B):
    return A0 * np.exp(-E * x_norm) + B

rows = []
for tau in TAU_GRID:
    T_eff = T_TRIAL - tau
    if T_eff < 10:
        continue
    ts = int(round(tau * FS))
    yw = y[ts:]; tw = time_eda[ts:]
    xe = (tw - tau) / T_eff
    popt, _ = curve_fit(R_matlab, xe, yw, p0=[yw[0]-yw[-1], yw[-1]],
                        bounds=([0, -50], [50, 50]), maxfev=20000)
    yfit = R_matlab(xe, *popt)
    r2 = 1 - np.sum((yw - yfit)**2) / np.sum((yw - yw.mean())**2)
    rmse = np.sqrt(np.mean((yw - yfit)**2))
    rows.append(dict(tau=tau, A0=popt[0], B=popt[1], R2=r2, RMSE=rmse))
eda_sweep = pd.DataFrame(rows)
best_eda = eda_sweep.loc[eda_sweep.R2.idxmax()]
tau_EDA   = float(best_eda.tau)
T_eff_EDA = T_TRIAL - tau_EDA
A0_EDA    = float(best_eda.A0)
B_EDA     = float(best_eda.B)
R2_eda    = float(best_eda.R2)
RMSE_eda  = float(best_eda.RMSE)
peak_EDA  = tau_EDA + T_eff_EDA / E

ts = int(round(tau_EDA * FS))
t_fit = time_eda[ts:]
y_fit = y[ts:]
x_fit = (t_fit - tau_EDA) / T_eff_EDA
yhat  = R_matlab(x_fit, A0_EDA, B_EDA)
res_eda = y_fit - yhat
sens = -A0_EDA * x_fit * np.exp(-E * x_fit)
sens_trough = -A0_EDA / (E ** 2)

print(f'EDA: tau = {tau_EDA:.2f} s, T_eff = {T_eff_EDA:.2f} s, '
      f'R^2 = {R2_eda:.4f}, peak_sens = {peak_EDA:.2f} s')

# =========================================================================
# Click cohorts
# =========================================================================
kde_df = pd.read_csv(MKDE)
xi = kde_df.t_s.to_numpy()

df_click = pd.read_csv(CLICK)
df_reas  = pd.read_csv(REAS)

# Strict (OriginalTarget in avatar/shadow)
mask_strict = (df_click.Clicked == 1) & df_click.ClickTarget.isin(['avatar', 'shadow'])
clicks_strict = df_click.loc[mask_strict, 'ClickTime_s'].to_numpy()
clicks_strict = clicks_strict[np.isfinite(clicks_strict) & (clicks_strict >= 0) & (clicks_strict < T_TRIAL)]

# Augmented (strict + ex-none clicks whose nearest prior contact is moving)
mask_rec = ((df_reas.Clicked == 1) & (df_reas.OriginalTarget == 'none') &
            df_reas.AssignedTarget.isin(['avatar', 'shadow']))
clicks_rec = df_reas.loc[mask_rec, 'ClickTime_s'].to_numpy()
clicks_rec = clicks_rec[np.isfinite(clicks_rec) & (clicks_rec >= 0) & (clicks_rec < T_TRIAL)]
clicks_aug = np.concatenate([clicks_strict, clicks_rec])

# Reassigned (nearest-prior-contact, moving-target-only subset)
mask_reas = ((df_reas.Clicked == 1) & df_reas.AssignedTarget.isin(['avatar', 'shadow']))
clicks_reas = df_reas.loc[mask_reas, 'ClickTime_s'].to_numpy()
clicks_reas = clicks_reas[np.isfinite(clicks_reas) & (clicks_reas >= 0) & (clicks_reas < T_TRIAL)]

cohorts = {
    'strict':     dict(clicks=clicks_strict, f=kde_df.f_moving.to_numpy(),
                       label='strict (avatar + shadow)',      n=len(clicks_strict)),
    'augmented':  dict(clicks=clicks_aug,    f=kde_df.f_moving_augmented.to_numpy(),
                       label='augmented (+ 78 rescued none)', n=len(clicks_aug)),
    'reassigned': dict(clicks=clicks_reas,   f=kde_df.f_moving_reassigned.to_numpy(),
                       label='reassigned (nearest prior)',    n=len(clicks_reas)),
}
print()
for k, v in cohorts.items():
    print(f'  {k:10s}: N = {v["n"]}')

BIN = 2.0
edges = np.arange(0, T_TRIAL + BIN, BIN)
centres = edges[:-1] + BIN / 2

def shape(t, tau, T_eff):
    x = np.maximum((t - tau) / T_eff, 0.0)
    out = x * np.exp(-E * x)
    out[t < tau] = 0.0
    return out

def fit_amp(tau, T_eff, xi, f):
    s = shape(xi, tau, T_eff)
    A = np.dot(s, f) / np.dot(s, s)
    fit = A * s
    ss_res = np.sum((f - fit)**2)
    ss_tot = np.sum((f - f.mean())**2)
    return A, fit, 1 - ss_res/ss_tot, np.sqrt(ss_res/len(f))

def sweep_tau(f):
    rows = []
    for tau in TAU_GRID:
        T_eff = T_TRIAL - tau
        if T_eff < 10:
            continue
        A, fit, r2, rmse = fit_amp(tau, T_eff, xi, f)
        rows.append(dict(tau=tau, T_eff=T_eff, A=A, R2=r2, RMSE=rmse))
    return pd.DataFrame(rows)

for k, v in cohorts.items():
    sw = sweep_tau(v['f'])
    best = sw.loc[sw.R2.idxmax()]
    tau_c = float(best.tau); T_c = T_TRIAL - tau_c
    A_c, fit_c, r2_c, rmse_c = fit_amp(tau_c, T_c, xi, v['f'])
    v['sweep'] = sw
    v['tau']   = tau_c
    v['T_eff'] = T_c
    v['A']     = A_c
    v['fit']   = fit_c
    v['R2']    = r2_c
    v['RMSE']  = rmse_c
    v['peak']  = tau_c + T_c / E
    v['counts'], _ = np.histogram(v['clicks'], bins=edges)
    v['kde_scale'] = v['n'] * BIN
    print(f'  {k:10s}: tau = {tau_c:.2f} s, peak = {v["peak"]:.2f} s, '
          f'R^2 = {r2_c:.4f}')

# =========================================================================
# Plot
# =========================================================================
fig = plt.figure(figsize=(13.5, 9.0))
gs = fig.add_gridspec(3, 2, height_ratios=[2.2, 1.0, 1.1],
                      hspace=0.55, wspace=0.22)

# --- EDA panel (top left) ---
axL = fig.add_subplot(gs[0, 0])
axL.fill_between(time_eda, y - sem, y + sem, color='0.75', alpha=0.55,
                 linewidth=0, label=f'SEM across participants (N = {n_part})')
axL.plot(time_eda, y, color='0.15', lw=1.4,
         label='mean EDA (per-participant grand average)')
axL.plot(t_fit, yhat, color='#1f77b4', lw=2.2,
         label=(r'$R(x) = A_0\,e^{-e\,x} + B$'
                f'\n    $\\tau_{{\\mathrm{{EDA}}}}$ = {tau_EDA:.2f} s,  '
                f'$R^2$ = {R2_eda:.3f}'))
axL.axvline(tau_EDA, color='#1f77b4', lw=0.8, ls=':', alpha=0.7)
axL.set_ylabel(r'EDA ($\mu$S)')
axL.set_title(r'Task EDA on trial horizon ($x \in [0,1]$, rate $= e$)')
axL.spines[['top']].set_visible(False)

axLr = axL.twinx()
axLr.plot(t_fit, sens, color='#2ca02c', lw=2.0, ls='-.',
          label=(r'$\partial R/\partial\lambda|_{\lambda=e}'
                 r' = -A_0\,x\,e^{-e\,x}$'))
axLr.axhline(0, color='#2ca02c', lw=0.4, alpha=0.4)
axLr.set_ylabel(r'$\partial R/\partial\lambda$ at $\lambda=e$  ($\mu$S)',
                color='#2ca02c')
axLr.tick_params(axis='y', labelcolor='#2ca02c')
axLr.spines[['top']].set_visible(False)
axLr.spines['right'].set_color('#2ca02c')
axLr.set_ylim(top=abs(sens_trough) * 0.1, bottom=sens_trough * 1.25)
h1, l1 = axL.get_legend_handles_labels()
h2, l2 = axLr.get_legend_handles_labels()
axL.legend(h1 + h2, l1 + l2, loc='upper right', frameon=False, fontsize=9.0)

# --- Click panel (top right) — AUGMENTED cohort ---
axR = fig.add_subplot(gs[0, 1])
v = cohorts['augmented']
axR.bar(centres, v['counts'], width=BIN, color='#6A85B5', alpha=0.30,
        edgecolor='white', lw=0.5, label=f'click histogram (N = {v["n"]})')
axR.plot(xi, v['f'] * v['kde_scale'], color='#1f3a6a', lw=1.6,
         label='Matlab ksdensity (reflection)')
axR.plot(xi, v['fit'] * v['kde_scale'], color='#d62728', lw=2.2,
         label=(r'$A\cdot x\,e^{-e\,x}$ fit, $\tau$ free'
                f'\n    $\\tau_{{\\mathrm{{click}}}}$ = {v["tau"]:.2f} s,  '
                f'$R^2$ = {v["R2"]:.3f},  peak {v["peak"]:.1f} s'))
axR.axvline(v['tau'],  color='#d62728', lw=0.8, ls=':', alpha=0.7)
axR.axvline(v['peak'], color='#d62728', lw=0.8, ls=':', alpha=0.4)
axR.axvline(tau_EDA,   color='#1f77b4', lw=0.8, ls=':', alpha=0.7)
axR.set_ylabel('click count (2 s bins)')
axR.set_title(r'Click response times — augmented moving-target cohort'
              f' (avatar + shadow + rescued none, N = {v["n"]})')
axR.spines[['top', 'right']].set_visible(False)
axR.legend(loc='upper right', frameon=False, fontsize=9.0)

# --- Residuals (middle row) ---
axLb = fig.add_subplot(gs[1, 0])
axLb.axhline(0, color='k', lw=0.8, alpha=0.7)
axLb.plot(t_fit, res_eda, color='#d62728', lw=1.3,
          label=f'residual   (RMSE = {RMSE_eda:.3f} μS)')
axLb.axvline(tau_EDA, color='#1f77b4', lw=0.8, ls=':', alpha=0.7)
axLb.set_xlabel('trial time $t$ (s)')
axLb.set_ylabel(r'residual ($\mu$S)')
axLb.set_title(f'EDA residual (fit window $t \\geq {tau_EDA:.2f}$ s)')
axLb.legend(loc='upper right', frameon=False, fontsize=9.0)
axLb.spines[['top', 'right']].set_visible(False)

axRb = fig.add_subplot(gs[1, 1])
axRb.axhline(0, color='k', lw=0.8, alpha=0.7)
mask_fit = xi >= v['tau']
axRb.plot(xi[mask_fit], ((v['f'] - v['fit']) * v['kde_scale'])[mask_fit],
          color='#d62728', lw=1.3,
          label=f'residual   (RMSE = {v["RMSE"] * v["kde_scale"]:.2f} counts)')
axRb.axvline(v['tau'], color='#d62728', lw=0.8, ls=':', alpha=0.7)
axRb.axvline(tau_EDA,  color='#1f77b4', lw=0.8, ls=':', alpha=0.7)
axRb.set_xlabel('trial time $t$ (s)')
axRb.set_ylabel('residual (count units)')
axRb.set_title(r'Click residual (augmented cohort, fit window '
               f'$t \\geq {v["tau"]:.2f}$ s)')
axRb.legend(loc='upper right', frameon=False, fontsize=9.0)
axRb.spines[['top', 'right']].set_visible(False)

# --- Cohort sensitivity (bottom row, spans both columns) ---
axS1 = fig.add_subplot(gs[2, 0])
axS2 = fig.add_subplot(gs[2, 1])

# Bottom-left: R^2(tau) curves for all three cohorts
col = {'strict': '#888888', 'augmented': '#d62728', 'reassigned': '#1f77b4'}
for k, v in cohorts.items():
    axS1.plot(v['sweep'].tau, v['sweep'].R2, color=col[k], lw=1.6,
              label=f'{k} (N = {v["n"]},  $\\tau^*$ = {v["tau"]:.2f} s,  '
                    f'$R^2$ = {v["R2"]:.3f})')
    axS1.scatter([v['tau']], [v['R2']], color=col[k], s=25, zorder=5)
axS1.axvline(tau_EDA, color='#2ca02c', lw=1.0, ls='--', alpha=0.8,
             label=rf'$\tau_{{\mathrm{{EDA}}}}$ = {tau_EDA:.2f} s')
axS1.set_xlabel(r'lag $\tau$ (s)')
axS1.set_ylabel(r'$R^2$ of click-KDE fit')
axS1.set_title(r'Cohort sensitivity: $R^2(\tau)$ for three definitions '
               'of the moving-target set')
axS1.spines[['top', 'right']].set_visible(False)
axS1.legend(loc='lower center', frameon=False, fontsize=8.5)
axS1.set_xlim(0, 10)

# Bottom-right: the three KDEs overlaid to show how similar they are
for k, v in cohorts.items():
    axS2.plot(xi, v['f'], color=col[k], lw=1.5,
              label=f'{k} (N = {v["n"]})')
axS2.set_xlabel('trial time $t$ (s)')
axS2.set_ylabel('KDE (density, 1/s)')
axS2.set_title('Matlab KDE by cohort (reflection-corrected)')
axS2.spines[['top', 'right']].set_visible(False)
axS2.legend(loc='upper right', frameon=False, fontsize=8.5)
axS2.set_xlim(0, T_TRIAL)

for ax in (axL, axR, axLb, axRb):
    ax.set_xlim(0, T_TRIAL)

fig.suptitle(
    r'EDA and clicks fit independently — augmented moving-target cohort: '
    rf'$\tau_{{\mathrm{{EDA}}}} = {tau_EDA:.2f}$ s,  '
    rf'$\tau_{{\mathrm{{click}}}} = {cohorts["augmented"]["tau"]:.2f}$ s  '
    rf'($\Delta = {cohorts["augmented"]["tau"] - tau_EDA:+.2f}$ s)',
    fontsize=11.5)

fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(OUT_PNG, dpi=200, bbox_inches='tight')
fig.savefig(OUT_PDF, bbox_inches='tight')
print(f'\nWrote {OUT_PNG}')
print(f'Wrote {OUT_PDF}')
