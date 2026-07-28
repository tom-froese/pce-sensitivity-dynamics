# Sensitivity Dynamics in Perceptual Crossing

Code and data for:

> **Clock-free optimal stopping in decision-making: Stochastic sensitivity dynamics predict behavioral, neural, and perceptual transition timing**
>
> Tom Froese (2026). *In preparation* (aimed at *npj Complexity*).

## Overview

This repository contains all code needed to reproduce the data figures reported in the paper. We show that participants in a perceptual crossing experiment (PCE) exhibit sensitivity dynamics consistent with the derivative of an exponential reliability decay with rate parameter Λ = *e*:

    R(x) = exp(-e * x)                # reliability of the prepared state
    S(x) = dR/dΛ = -x * exp(-e * x)   # signed rate sensitivity (a trough)

The signed sensitivity `S(x)` reaches its extremum — a trough — at *x* = 1/*e*, the optimal stopping point of the 1/*e*-law. Each empirical readout is an **affine observation** of this geometric function, `A·S(x) + B`, with the readout's direction carried in the fitted gain `A`: a readout that *peaks* (click times, the 1/f aperiodic exponent) simply has `A < 0`, while the global scalp potential troughs with `A > 0`. All four landmark *x* = 1/*e* — across behavioral responses (click times), autonomic arousal (electrodermal activity), neural activity (global scalp potential), and perceptual awareness (PAS ratings).

## Quick Start

### One-command reproduction

In MATLAB:
```matlab
reproduce            % from preprocessed data (default)
reproduce('raw')     % from raw data (requires data/raw/)
```

### Prerequisites

- MATLAB R2022a or later with Statistics and Optimization toolboxes
- [EEGLAB](https://sccn.ucsd.edu/eeglab/) (for topographic scalp maps)
- Python 3.9+ for Figure 2 (panels B–D) and the aperiodic-exponent pipeline.
  Install the Python dependencies with `pip install -r requirements.txt`
  (key packages: `mne` + `fooof`, plus the usual
  `numpy`/`pandas`/`scipy`/`matplotlib`).

## Data

### Preprocessed data

The smaller derived files (click times, haptics, EDA, and the global-scalp-potential MAT files) are included in this repository under `data/preprocessed/`. The larger EEG inputs are archived on Zenodo — concept DOI [10.5281/zenodo.19425013](https://doi.org/10.5281/zenodo.19425013), which always resolves to the latest version:

| File on Zenodo | Unzip into | Contents |
|---|---|---|
| `preprocessed.zip` | `data/` | the preprocessed CSV/MAT files, including the 64-channel `allchannel_data.mat` |
| `eeg_task_raw_fif.zip` (~4.3 GB) | `data/preprocessed/EEG/` | the 62 cleaned, **preprocessed** continuous EEG recordings (one per participant), as per-dyad `pceXX/` folders, read by the aperiodic-exponent (1/f) pipeline |

> **Note on the `.fif` naming.** The per-participant files are named `pceXX_PY_task-raw.fif`. The `-raw.fif` suffix is MNE-Python's required label for a *continuous* recording (an `mne.Raw` object); it does **not** mean the data are unprocessed — they are bandpass-filtered 1–40 Hz, resampled to 250 Hz, bad-channel interpolated, and average-referenced.

### Raw data (on OSF, required only for re-running preprocessing)

The raw experimental data is archived on OSF as part of the [Perceptual Crossing Dataset](https://osf.io/47n3p). To reproduce the full pipeline from scratch, download the raw data and organize it into `data/raw/`:

| Directory | OSF source | Size |
|-----------|------------|------|
| `data/raw/Behavior/` | [osf.io/7hfec](https://osf.io/7hfec) | ~1.8 GB |
| `data/raw/EDA/` | [osf.io/47n3p](https://osf.io/47n3p) — Peripheral Physiological Data / EDA.zip | ~38 MB |
| `data/raw/EEG/` | [osf.io/47n3p](https://osf.io/47n3p) — Raw EEG Hyperscanning Data (32 per-dyad zips) | ~15 GB |

## Repository Structure

```
pce-sensitivity-dynamics/
├── reproduce.m                         # Master reproduction script
├── code/
│   ├── _config.py                       # Shared Python config (paths + EEG)
│   ├── analysis/                       # Figure generation and statistics
│   │   ├── plotFigure1_Behavioral.m    # Fig 1: clicks + haptics + EDA
│   │   ├── plotFigure2_Neural.py       # Fig 2: GSP + scalp maps + 1/f exponent
│   │   ├── computeAperiodicExponent.py # Per-(subj × bin) FOOOF exponent (Fig 2D)
│   │   ├── fitExponentSensitivity.py   # S(x) fit + bootstrap CI on exponent peak
│   │   ├── computePASCrossover.m       # PAS 4/3 logistic crossover + click-time/PAS Spearman
│   │   ├── plotFigure3_Perceptual.m    # Fig 3: sensitivity + PAS
│   │   ├── computeClickPAS.py          # Per-click click-time + PAS table (ClickPAS.csv)
│   │   ├── computeEDAFreeLambda.py     # EDA exponential-vs-linear test (free-lambda, dR2)
│   │   ├── computeWithinTrialComplexity.py  # Within-trial multichannel-LZc slope + spectral entropy
│   │   ├── computeDecouplingBF.py      # Neural<->phenomenal decoupling Bayes factors
│   │   ├── _eeg_io.py                  # Dataset + MNE loaders (Python pipeline)
│   │   └── _exponent_common.py         # Shared constants for the exponent code
│   └── preprocessing/                  # Data extraction pipelines
│       ├── preprocessClicks.m          # Click response times from raw trials
│       ├── preprocessHaptics.m         # Haptic feedback from raw trials
│       ├── preprocessEDA.m             # Electrodermal activity from raw EDA
│       ├── preprocessGSP.m            # Global scalp potential from raw EEG
│       ├── computeGSPStats.m           # Hierarchical sensitivity model fits
│       ├── extractAllChannels.m        # All 64 EEG channels from raw MAT
│       ├── extractParietalHemispheres.m  # L/R parietal cluster extraction
│       ├── computePerChannelFits.m     # Per-channel free-tau + locked-tau fits
│       ├── computePASProportions.m     # Unsmoothed disjoint-bin PAS proportions
│       └── preprocessEEGForExponent.py # Raw .mat → 250 Hz cleaned .fif (Python)
├── data/
│   ├── preprocessed/                   # Tracked in git (included in repo)
│   │   ├── ClickTimes/                 # Behavioral responses
│   │   ├── ClickPAS/                   # Per-click click time + PAS rating (derived; see computeClickPAS.py)
│   │   ├── EDA/                        # Electrodermal activity (task + rest)
│   │   ├── EEG/                        # Neural data (MAT files)
│   │   └── Haptics/                    # Haptic feedback time series
│   └── raw/                            # NOT tracked (download from OSF)
├── results/                            # Generated figures and intermediate CSVs
├── README.md
├── LICENSE                             # MIT
└── CITATION.cff
```

## Preprocessing Pipeline

The master script `reproduce.m` runs these steps in order:

| Step | Script | Output |
|------|--------|--------|
| 1a | `preprocessClicks.m` | `data/preprocessed/ClickTimes/` |
| 1b | `preprocessHaptics.m` | `data/preprocessed/Haptics/` |
| 1c | `preprocessEDA.m` | `data/preprocessed/EDA/` |
| 1d | `preprocessGSP.m` | `data/preprocessed/EEG/globalScalpPotential_data.mat` |
| 1e | `computeGSPStats.m` | `data/preprocessed/EEG/globalScalpPotential_stats.mat` |
| 1f | `extractAllChannels.m` | `data/preprocessed/EEG/allchannel_data.mat` |
| 1g | `extractParietalHemispheres.m` | `data/preprocessed/EEG/parietal_hemisphere_data.mat` |
| 1h | `computePerChannelFits.m` | `results/Figure2_perchannel_fits.csv` |
| 1i | `computePASProportions.m` | `results/Figure3_pas_proportions.csv` |
| 1j | `computeClickPAS.py` | `data/preprocessed/ClickPAS/ClickPAS.csv` (merges click times with the raw PAS questionnaires; needs `data/raw/Behavior/`) |
| 1k | `preprocessEEGForExponent.py` | `data/preprocessed/EEG/pceXX/pceXX_PY_task-raw.fif` (per dyad/participant; Python via MNE — minimal chain: bandpass 1–40 Hz [Nyquist alias safety + drift control] → resample 250 Hz → bad-channel LOF interp → average reference. The recording-level chain — 60 Hz notch, FCz common reference, 1000 Hz sample, 0.016–1000 Hz analog cutoff — is documented in Lerique et al. 2024.) |

### Derived statistics

| Step | Script | Output |
|------|--------|--------|
| 2a | `computePerChannelFits.m`     | `results/Figure2_perchannel_fits.csv` |
| 2b | `computeAperiodicExponent.py` | `results/aperiodic_exponent_per_participant{,_band-2-20}.csv` (per-(subj × bin) FOOOF aperiodic exponent on 2-40 Hz + 2-20 Hz EMG-reduced band) |
| 2c | `fitExponentSensitivity.py`   | `results/aperiodic_exponent_within_trial.{csv,png,json}` (cohort mean ± SEM + S(x) fit + bootstrap peak CI; powers Panel D of Figure 2) |
| 2d | `computeClickPAS.py`          | `data/preprocessed/ClickPAS/ClickPAS.csv` (per-click click time + PAS rating; merges `ClickResponseTimes.csv` with the raw PAS questionnaires — needs `data/raw/Behavior/`; bundled in the repo so the two PAS analyses below run without raw) |
| 2e | `computePASCrossover.m`       | `results/Figure3_crossover_stats.csv` (adds the within-participant click-time/PAS **Spearman** `ρ=−0.12`, `t(57)=−2.77`, `p=0.008`, `N=58` alongside the crossover stats) |
| 2f | `computeEDAFreeLambda.py`     | `results/eda_free_lambda.csv` (EDA exponential-vs-linear decisive test: free-`λ`≈2.62≈*e*, `ΔR²`=R²[exp@e]−R²[linear]=+0.07) |
| 2g | `computeWithinTrialComplexity.py` | `results/within_trial_complexity_pertrial.csv` (per-trial multichannel-LZc slope `mc_slope`) + `results/within_trial_windows.csv` (per-window global `specent_global`) — from the cleaned `.fif` EEG; the compute-heavy step (~10 min) |
| 2h | `computeDecouplingBF.py`      | `results/decoupling_bayes_factors.csv` (three neural↔phenomenal decoupling JZS Pearson Bayes factors: BF₁₀=0.17 [coupling, n=58], 0.25 [timing, n=25], 0.28 [magnitude, n=58] — all favor the null) |

### Manuscript figures

| Step | Script | Output |
|------|--------|--------|
| 3a | `plotFigure1_Behavioral.m` | `results/Figure1_Behavioral.{png,pdf}` |
| 3b | `plotFigure2_Neural.py`    | `results/Figure2_Neural.{png,pdf}` |
| 3c | `computePASCrossover.m`    | `results/Figure3_crossover_stats.csv` |
| 3d | `plotFigure3_Perceptual.m` | `results/Figure3_Perceptual.{png,pdf}` |

## Dataset

32 dyads (64 participants) performed the perceptual crossing experiment. Each session comprised 18 trials of 60 seconds and 4 rest periods of 180 seconds. Simultaneous recordings include 64-channel EEG (1000 Hz), electrodermal activity, respiration, haptic feedback, behavioral responses (button presses), and subjective perceptual awareness ratings (PAS, 1-4 scale).

## Related Repositories

For the full analysis pipeline including respiration fits, see [icdl2026-embodied-anticipation](https://github.com/tom-froese/icdl2026-embodied-anticipation).

## License

MIT. See [LICENSE](LICENSE).
