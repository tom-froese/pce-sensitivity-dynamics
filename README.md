# Sensitivity Dynamics in Perceptual Crossing

Code and data for:

> **Clock-free optimal stopping in social cognition: Stochastic sensitivity dynamics predict behavioral, neural, and perceptual transition timing**
>
> Tom Froese (2026). *Under review.*

## Overview

This repository contains all code needed to reproduce the three figures reported in the paper. We show that participants in a perceptual crossing experiment (PCE) exhibit sensitivity dynamics consistent with the derivative of an exponential reliability decay with rate parameter lambda = *e*:

    S(x) = |dR/dlambda| = x * exp(-e * x)

This sensitivity function peaks at *x* = 1/*e*, the optimal stopping point of the 1/*e*-law, and is reflected in behavioral responses (click times), autonomic arousal (electrodermal activity), neural activity (global scalp potential), and perceptual awareness (PAS ratings).

## Quick Start

From MATLAB (R2023b or later), run each figure script from `code/analysis/`:

```matlab
cd code/analysis
plotFigure1_Behavioral   % Figure 1: Click times + haptics + EDA
plotFigure2_Neural       % Figure 2: ROI sensitivity model fits
plotFigure3_Perceptual   % Figure 3: Sensitivity schematic + PAS clarity
```

All three scripts use relative paths and output 600 dpi PNG files to `results/`.

## Data

### Preprocessed data (included in this repository)

The figure scripts load these precomputed files directly:

| Directory | Contents | Size |
|-----------|----------|------|
| `data/preprocessed/ClickTimes/` | Click response times (CSV + JSON sidecar) | 32 KB |
| `data/preprocessed/Haptics/` | Haptic feedback time series (gzipped CSV) | 18 MB |
| `data/preprocessed/EDA/` | Electrodermal activity (CSV + JSON sidecars) | 50 MB |
| `data/preprocessed/PAS/` | PAS sensitivity and crossover analysis (MAT) | 40 KB |
| `data/preprocessed/EEG/` | Global scalp potential data and statistics (MAT) | 7 MB |

### Raw data (on Zenodo, required only for re-running preprocessing)

To reproduce the full pipeline from scratch, download the raw data archive from Zenodo, unzip it into `data/raw/`:

| Directory | Contents | Size |
|-----------|----------|------|
| `data/raw/Behavior/` | Raw trial CSVs and PAS questionnaires (32 dyads) | ~1.8 GB |
| `data/raw/EDA/` | Raw EDA .mat files (1000 Hz, task + rest) | ~107 MB |
| `data/raw/EEG/` | Raw EEG .mat files (64 ch, 1000 Hz, task + rest) | ~15 GB |

Zenodo DOI: [10.5281/zenodo.19425014](https://doi.org/10.5281/zenodo.19425014)

After downloading, your `data/raw/` directory should contain `Behavior/`, `EDA/`, and `EEG/` subfolders, each with `pce*` experiment folders inside.

## Toolbox Requirements

- Signal Processing Toolbox (filtering, spectral estimation)
- Statistics and Machine Learning Toolbox (logistic regression, permutation tests)
- Optimization Toolbox (curve fitting in Figure 1 Panels B and C)

No third-party packages are required.

## Repository Structure

```
pce-sensitivity-dynamics/
├── code/
│   ├── analysis/                    # Figure generation scripts
│   │   ├── plotFigure1_Behavioral.m # Fig 1: clicks + haptics + EDA
│   │   ├── plotFigure2_Neural.m     # Fig 2: ROI sensitivity fits
│   │   └── plotFigure3_Perceptual.m # Fig 3: sensitivity + PAS
│   └── preprocessing/               # Data extraction pipelines
│       ├── preprocessClicks.m       # Click response times from raw trials
│       ├── preprocessHaptics.m      # Haptic feedback from raw trials
│       ├── preprocessEDA.m          # Electrodermal activity from raw EDA
│       ├── preprocessPAS.m          # PAS ratings from questionnaires
│       ├── preprocessGSP.m         # Global scalp potential from raw EEG
│       ├── computeGSPStats.m        # Hierarchical sensitivity model fits
│       ├── computePASSensitivity.m  # PAS early vs. late analysis
│       └── computePASCrossover.m    # PAS 4/3 crossover + ROI mapping
├── data/
│   ├── preprocessed/                # Tracked in git (included in repo)
│   │   ├── ClickTimes/              # Behavioral responses
│   │   ├── EDA/                     # Electrodermal activity (task + rest)
│   │   ├── EEG/                     # Neural data (MAT files)
│   │   ├── Haptics/                 # Haptic feedback time series
│   │   └── PAS/                     # Perceptual awareness analysis
│   └── raw/                         # NOT tracked (download from Zenodo)
│       ├── Behavior/                # Raw trial CSVs + questionnaires
│       ├── EDA/                     # Raw EDA .mat files
│       └── EEG/                     # Raw EEG .mat files
├── results/                         # Publication figures (600 dpi PNG)
├── README.md
├── LICENSE                          # MIT
└── CITATION.cff
```

## Preprocessing

The figure scripts load precomputed data files included in `data/preprocessed/`. The preprocessing scripts document how these files were generated from raw experimental data.

To reproduce from scratch, download the raw data from Zenodo into `data/raw/`, then run from `code/preprocessing/`:

```matlab
cd code/preprocessing

% Step 1: Extract behavioral data (independent, from raw trial CSVs)
preprocessClicks       % -> data/preprocessed/ClickTimes/
preprocessHaptics      % -> data/preprocessed/Haptics/

% Step 2: Preprocess physiological signals
preprocessEDA          % -> data/preprocessed/EDA/

% Step 3: Extract perceptual awareness ratings
preprocessPAS          % -> data/preprocessed/PAS/

% Step 4: Extract global scalp potential from raw EEG
preprocessGSP          % -> data/preprocessed/EEG/globalScalpPotential_data.mat

% Step 5: Compute sensitivity model statistics
computeGSPStats        % -> data/preprocessed/EEG/globalScalpPotential_stats.mat

% Step 6: PAS analysis (requires steps 1, 2, 5)
computePASSensitivity  % -> data/preprocessed/PAS/gsp_sensitivity_pas.mat
computePASCrossover    % -> data/preprocessed/PAS/gsp_pas_roi_crossover.mat
```

All preprocessing scripts use `mfilename('fullpath')` to resolve paths, so they work regardless of MATLAB's current directory.

## Figures

**Figure 1 -- Bodily Evidence for Reliability Decay Dynamics.** Three panels: (A) Click response-time distribution fitted by the sensitivity function S(x); (B) Haptic contact proportion saturating at 1/e; (C) Electrodermal activity following the reliability function R(x) = exp(-e*x).

**Figure 2 -- Neural Evidence: Sensitivity Model Fits by Brain Region.** Seven ROI time courses with continuous anterior-to-posterior color gradient, showing the spatial distribution of sensitivity dynamics in the global scalp potential.

**Figure 3 -- Perceptual Evidence: Sensitivity Dynamics and PAS.** Two panels: (A) Theoretical sensitivity of R(x) to rate perturbations, with rejection and selection phases; (B) Moving-window PAS ratings with bootstrap crossover CI and ROI trough markers.

## Dataset

32 dyads (64 participants) performed the perceptual crossing experiment. Each session comprised 18 trials of 60 seconds and 4 rest periods of 180 seconds. Simultaneous recordings include 64-channel EEG (1000 Hz), electrodermal activity, respiration, haptic feedback, behavioral responses (button presses), and subjective perceptual awareness ratings (PAS, 1-4 scale).

## Related Repositories

For the full analysis pipeline including respiration fits, see [icdl2026-embodied-anticipation](https://github.com/tom-froese/icdl2026-embodied-anticipation).

## License

MIT. See [LICENSE](LICENSE).
