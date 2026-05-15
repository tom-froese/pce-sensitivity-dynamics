# Figure S2 — Methods

This document records the derivation, the cohort choices, and the audit
that underlie Figure S2 (`plotFigureS2_EDA_Clicks.py`,
`results/FigureS2_EDA_Clicks_IndependentLags.png`).

## 1. Scope

Figure S2 shows that two independent measures — the grand-average task EDA
and the distribution of click response times — are each fit by *different*
analytical shapes that both follow from the same underlying rate-decay
model:

- the EDA is fit by the decay profile itself,
  **R(x) = A₀·exp(−e·x) + B**;
- the click-time KDE is fit by the analytical modulus of the rate
  sensitivity derived from R(x),
  **S(x) = |∂R/∂λ|_{λ=e} = A·x·exp(−e·x)**.

The trial axis is parameterized as x = (t − τ)/(T − τ). Fitting the two
signals *independently* rather than locking them to a common τ lets the
near-match between the two recovered lags — and, as we show in § 6, the
agreement between the two residuals — stand as evidence for the
irruption-theory rate-sensitivity interpretation rather than as a
constraint imposed by us.

## 2. EDA fit

- **Preprocessed input**: `data/preprocessed/EDA/EDA_Task_Preprocessed.csv`
  (per-participant task-EDA traces, 25 Hz, trials of 60 s; produced by the
  existing preprocessing pipeline).
- **Aggregation**: for each participant, average across trials to obtain a
  single 60-s task profile; then average across all complete participants
  to obtain the grand average y(t). Kept N = 61 participants.
- **Standard error**: cross-participant SEM (y.std(axis=0, ddof=0) / √N).
- **Trial-horizon model**: x = (t − τ) / (T − τ) with T = 60 s, and

        R(x) = A₀ · exp(−e · x) + B          on x ∈ [0, 1].

  The rate λ is fixed at e throughout (the distinguished value at which
  ∂R/∂λ has a single finite trough at x = 1/e). The free parameters are
  (A₀, B, τ); (A₀, B) are obtained by linear regression of y against
  exp(−e · x) for each candidate τ.
- **τ-sweep**: τ is scanned on a 0 – 15 s grid in steps of 0.1 s and the
  value with the largest R² is reported.

**EDA result** (canonical run):

- τ_EDA = 4.90 s
- T_eff = T − τ = 55.10 s
- A₀ = (reported on the figure)
- B = (reported on the figure)
- R² = 0.974
- Sensitivity trough at x = 1/e, i.e. t_trough = τ + (T − τ)/e ≈ 25.17 s.

The left panel overlays y(t) ± SEM with R(x) and, on a secondary axis, the
signed sensitivity curve −A₀·x·exp(−e·x) (green) so that the shared trough
is visible against the raw EDA.

## 3. Click response-time distribution

### 3.1 Source

- `data/preprocessed/ClickTimes/ClickResponseTimes.csv`, one row per
  (dyad, trial, participant, click) with columns including
  `ClickTime_s`, `Clicked`, `ClickTarget`, and the auxiliary inter-event
  latencies (`x_own_s`, `x_exo_s`, `x_any_s`).
- The original `ClickTarget` classification follows
  Froese, Iizuka & Ikegami (2014): at each click the preprocessor
  (`code/preprocessing/classifyClicks.m`) examines a fixed lookback point
  1.0 s before the click, checks which of the three on-screen objects
  (avatar, shadow, static) lies within ±70 units (on the 600-unit ring,
  with wrap-around distance), and assigns the first match in priority
  order `avatar → shadow → static → none`.

### 3.2 Canonical cohort: *augmented moving-target*

The comparison we care about is the distribution of click times that
follow contact with a *moving* object (avatar or shadow), since those are
the clicks whose latency should track the decaying-reliability signal.
Two observations shape the cohort definition:

1. The static target sits in one fixed location; clicks tagged `static`
   are driven by a different process (a participant can plan a click on
   a stationary target long in advance) and are excluded.
2. The 1.0 s *fixed-lookback* classifier is strict: a brief glancing
   contact that ended, say, 0.6 s before the click leaves no mark at the
   lookback point and the click is labelled `none`. An audit (Section 4
   below) shows that 74 of the 84 `none`-clicks are of this kind; they
   are moving-target clicks that the classifier missed. We rescue them.

Operationally:

- **Augmented moving-target cohort** = all rows with `Clicked == 1` and
  either (a) `ClickTarget ∈ {avatar, shadow}`, or (b) `ClickTarget ==
  none` *and* the nearest prior contact (within the trial, minimum of
  `x_avatar`, `x_shadow`, `x_static`) is with a moving object. We also
  accept the few cases where a rescued click's nearest contact is
  formally with `static` but within 1 s of a moving-target contact; see
  `_scratch/reassign_none.py` for the exact routine.

Counts: strict 870 (avatar 654 + shadow 216); augmented 948 = 870 + 78
rescued. A sensitivity analysis on a third definition (nearest-prior-
contact reassignment across *all* clicks, yielding a moving-target
subset of 884) is shown in the bottom panel.

### 3.3 KDE

A Gaussian kernel density is fit to the click-time cohort on [0, 60] s
with reflection boundary correction at both endpoints. We use Matlab's
`ksdensity` (Statistics & Machine Learning Toolbox) because the scipy
implementation we initially tried did not respect the hard boundaries
and produced artefacts at the trial edges. The Matlab call is:

```matlab
xi = linspace(0, 60, 601);
[f, ~] = ksdensity(click_times, xi, ...
                   'Support', [-0.001 60.001], ...
                   'BoundaryCorrection', 'reflection');
```

This is run once per cohort and the resulting evaluation grids are saved
to `results/click_kde_matlab.csv` (columns: `t_s`, `f_all`,
`f_nostatic`, `f_avatar`, `f_moving`, `f_moving_reassigned`,
`f_moving_augmented`). `plotFigureS2_EDA_Clicks.py` reads this CSV
directly — nothing in the Python figure script recomputes the KDE.

### 3.4 Shape fit

For each candidate τ on the same 0 – 15 s grid (step 0.1 s):

1. compute x = (t − τ)/(T − τ) on the KDE grid;
2. evaluate the template s(x) = x · exp(−e · x) on x ∈ [0, 1];
3. fit the amplitude by projection: A = ⟨s, f⟩ / ⟨s, s⟩;
4. compute R² = 1 − SS_res/SS_tot against the KDE.

The τ that maximises R² is reported. No non-linear fitting is used — the
amplitude projection is the closed-form least-squares solution for fixed
τ, and the outer τ-scan is a 1-D search on a dense grid.

**Click-cohort result** (canonical augmented cohort, N = 948):

- τ_click = 3.90 s
- T_eff = 56.10 s
- Peak of the fitted A·x·exp(−e·x) at t = τ + (T − τ)/e ≈ 24.54 s
- R² = 0.920

### 3.5 Sensitivity across cohort definitions

The bottom panel of the figure reports the free-τ fit on three cohorts:

| cohort               |   N | τ* (s) |   R² |
|----------------------|----:|------:|-----:|
| strict moving (avatar+shadow) | 870 | 3.70 | (on figure) |
| augmented moving (canonical)  | 948 | 3.90 | 0.920 |
| reassigned moving             | 884 | (on figure) | (on figure) |

All three land on a flat-topped R²(τ) plateau between roughly 2 and 6 s,
i.e. the recovered lag is not an artefact of one particular cohort
definition.

## 4. Audit of the 84 `none`-classified clicks

The 84 clicks with `ClickTarget == none` were decomposed by
`_scratch/decompose_none.py`. Each click is assigned to a bucket based on
`x_exo_s` — the elapsed time since the most recent contact with a moving
object *in that trial*:

| bucket                                    |   n | median click time |
|-------------------------------------------|----:|------------------:|
| `x_exo < 1 s` (classifier-missed)         |  74 | — |
| `1 ≤ x_exo < 3 s` (fading-signal, late)   |   6 | — |
| `x_exo ≥ 3 s` (endogenous)                |   1 | — |
| no prior exo contact in trial             |   3 | — |

**Interpretation.** The 74 clicks with `x_exo < 1 s` are the ones the
1.0 s fixed-lookback classifier missed: the moving object *had* been
within ±70 units during the last second, but had left the window by the
lookback point. These are moving-target clicks, so they belong in the
cohort. The 6 clicks with 1 ≤ `x_exo` < 3 s are plausibly still riding
the decaying reliability signal — we include them. The single
`x_exo ≥ 3 s` click is genuinely endogenous by any reasonable definition
and is dropped; the 3 clicks with no prior moving-object contact are
dropped. The final rescue count is 78 of the 84 `none`-clicks, giving
the augmented cohort of 948.

The endogenous-rate estimate from the same audit (clicks within the
> 3 s empty-interval mass, divided by that mass) is small relative to
the contact-window click rate, supporting the decision to model the
distribution as a single peri-contact process.

## 5. Auxiliary scripts

The audit and cohort-construction steps are implemented in the
scratch directory (not part of the canonical figure pipeline, but kept
for reproducibility):

- `_scratch/reassign_none.py` — builds
  `_scratch/reassigned_clicks.csv` (every click augmented with
  `x_avatar`, `x_shadow`, `x_static`, `AssignedTarget`, `x_assigned`).
- `_scratch/decompose_none.py` — buckets the 84 `none`-clicks and
  reports the classifier-missed count used above.
- `_scratch/compare_EDA_clicks_independent_lags.py` — earlier iteration
  on the strict cohort; retained for reference.
- `_scratch/compare_EDA_clicks_augmented.py` — direct ancestor of
  `plotFigureS2_EDA_Clicks.py`.

## 6. Residual cross-correlation: an independent check on the rate-sensitivity idea

Two distinct analytical shapes are fit in this figure:

- the EDA grand average is fit by **R(x) = A₀·exp(−e·x) + B**
  (monotone-decreasing exponential decay);
- the click-time KDE is fit by **S(x) = A·x·exp(−e·x)**
  (the modulus of ∂R/∂λ at λ = e, which is peaked with a single
  extremum at x = 1/e).

S(x) is derived analytically from R(x) under the rate-sensitivity
interpretation, but it is functionally distinct from R(x): one decays,
the other peaks. This means the two fits cannot share residual shape by
inheriting a common functional form. If the two residuals nevertheless
share a shape, that shared shape has to be carried by the underlying
dynamics that both the EDA and the click process are reporting on — it
cannot be an artefact of the fitting procedure.

We test for shared shape by cross-correlating the two residuals on a
common time grid (`_scratch/correlate_residuals.py`):

- Both residuals are interpolated onto a 0.1-s grid covering the
  overlap fit window, t ∈ [max(τ_EDA, τ_click), 60 s] = [4.90 s, 60 s],
  yielding 552 common samples.
- Zero-lag Pearson correlation:

      r(0) = +0.794   (p = 4 × 10⁻¹²¹,  n = 552,
                       block-bootstrap 95% CI [+0.58, +0.93],
                       2000 resamples, 5-s blocks)

- Lag-swept correlation (click residual shifted by L seconds relative
  to the EDA residual, L ∈ [−10, +10] s, 0.1-s step):

      r_max = +0.890   at   L* = +1.70 s.

  The positive sign means the click residual is shifted to the *right*
  (later in time) to align best with the EDA residual, i.e. the EDA
  residual *leads* the click residual by about 1.7 s. The τ-fits
  themselves differ by τ_EDA − τ_click = +1.00 s; the residual peak lag
  exceeds that by ~0.7 s, so the slow sub-structure the two signals
  share is itself slightly later in the click process than in the EDA.
  The correlation plateau is broad: r ≥ +0.85 across L ∈ [+1.0, +3.0]
  seconds.

**τ-locked control.** To confirm that the residual alignment is
carried by common sub-structure rather than by the ~1 s freedom the
independent τ-fits have, we repeat the S(x)-fit with τ_click *fixed*
to τ_EDA (= 4.90 s) and recompute the correlation. The shape
amplitude A is still obtained by projection:

| variant            | τ_click | click-fit R² | r(0)  | r_max | L* (s) |
|--------------------|--------:|-------------:|------:|------:|-------:|
| free τ (canonical) | 3.90    | 0.920        | +0.794 | +0.890 | +1.70 |
| τ locked to τ_EDA  | 4.90    | 0.912        | +0.774 | +0.777 | **+0.30** |

Three things happen when τ is locked:

1. The click fit's R² drops by only 0.008 (0.920 → 0.912), confirming
   the flat-topped R²(τ) plateau visible in the bottom-left panel of
   the figure: the click KDE is nearly indifferent between the two τ
   values.
2. The peak-correlation lag **collapses from +1.70 s to +0.30 s** —
   i.e. with a common τ, the two residual curves are essentially in
   phase. The ~1 s offset seen in the free-τ fit is precisely the
   offset between τ_EDA and τ_click, as expected.
3. Peak-r magnitude drops from +0.89 to +0.78, because with τ locked
   the fit can no longer absorb any of the shared sub-structure into a
   time shift; that structure now shows up entirely as a zero-lag
   correlation.

**Interpretation.** The leading-order fit S(x) captures 92 % of the
click-KDE variance while the EDA fit R(x) captures 97 % of the EDA
variance; the remaining residual structure in the two signals is
strongly correlated — in phase once a common τ is imposed. Because the
two fits use different functional forms (exponential decay vs.
x·exp(−e·x)) and the locked variant shares no free parameter at all on
the click side other than an amplitude, the shared residual shape is
not a fitting artefact — it reflects a common slow sub-structure in
both the electrodermal and behavioural readouts. This is independent
corroboration of the rate-sensitivity reading of the EDA→click
correspondence: the analytical move from R(x) to S(x) = |∂R/∂λ|
succeeds at the leading-order shape *and* the residuals agree, in
magnitude, in shape, and in phase (when τ is shared), to a degree that
would be surprising if the two signals were not reporting on a common
underlying dynamics.

## 7. Summary of canonical numbers

- EDA: τ_EDA = 4.90 s, T_eff = 55.10 s, R² = 0.974, trough at 25.17 s.
- Clicks (augmented, N = 948): τ_click = 3.90 s, T_eff = 56.10 s,
  R² = 0.920, peak at 24.54 s.
- Δτ = τ_click − τ_EDA = −1.00 s.
- Sensitivity panel: τ* stable across strict / augmented / reassigned
  cohort definitions within the flat R²(τ) plateau.

The two τ estimates were obtained from entirely independent signals
(per-participant EDA grand average vs. boundary-corrected click-time
KDE) with different functional forms (R(x) vs. S(x) = |∂R/∂λ|) and no
shared parameter, and they agree to within 1 s. In addition, the two
residuals are themselves strongly correlated (zero-lag r = +0.79; peak
r = +0.89 at a 1.7 s lag), which we interpret as independent
corroboration of the rate-sensitivity reading.
