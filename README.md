# Sensitivity Dynamics in Perceptual Crossing

Code and data for:

> **Clock-free optimal stopping in decision-making: Stochastic sensitivity dynamics predict behavioral, neural, and perceptual transition timing**
>
> Tom Froese (2026). *Under review.*

## Overview

This repository contains all code needed to reproduce the three figures reported in the paper. We show that participants in a perceptual crossing experiment (PCE) exhibit sensitivity dynamics consistent with the derivative of an exponential reliability decay with rate parameter lambda = *e*:

    S(x) = |dR/dlambda| = x * exp(-e * x)

This sensitivity function peaks at *x* = 1/*e*, the optimal stopping point of the 1/*e*-law, and is reflected in behavioral responses (click times), autonomic arousal (electrodermal activity), neural activity (global scalp potential), and perceptual awareness (PAS ratings).

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

## Data

### Preprocessed data (included in this repository; also archived on Zenodo)

| Directory | Contents | Size |
|-----------|----------|------|
| `data/preprocessed/ClickTimes/` | Click response times (CSV + JSON sidecar) | 32 KB |
| `data/preprocessed/Haptics/` | Haptic feedback time series (gzipped CSV) | 18 MB |
| `data/preprocessed/EDA/` | Electrodermal activity (CSV + JSON sidecars) | 50 MB |
| `data/preprocessed/EEG/` | Global scalp potential, per-channel, and parietal hemisphere data (MAT files) | 7 MB |

These files are also archived on Zenodo: [10.5281/zenodo.19425014](https://doi.org/10.5281/zenodo.19425014). To set up from the Zenodo archive, download `preprocessed.zip` and unzip it into `data/`.

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
│   ├── analysis/                       # Figure generation and statistics
│   │   ├── plotFigure1_Behavioral.m    # Fig 1: clicks + haptics + EDA
│   │   ├── plotFigure2_Neural.m        # Fig 2: EEG sensitivity fits
│   │   ├── computePASCrossover.m       # PAS 4/3 logistic crossover stats
│   │   └── plotFigure3_Perceptual.m    # Fig 3: sensitivity + PAS
│   └── preprocessing/                  # Data extraction pipelines
│       ├── preprocessClicks.m          # Click response times from raw trials
│       ├── preprocessHaptics.m         # Haptic feedback from raw trials
│       ├── preprocessEDA.m             # Electrodermal activity from raw EDA
│       ├── preprocessGSP.m            # Global scalp potential from raw EEG
│       ├── computeGSPStats.m           # Hierarchical sensitivity model fits
│       ├── extractAllChannels.m        # All 64 EEG channels from raw MAT
│       ├── extractParietalHemispheres.m  # L/R parietal cluster extraction
│       ├── computePerChannelFits.m     # Per-channel free-tau + locked-tau fits
│       └── computePASProportions.m     # Unsmoothed disjoint-bin PAS proportions
├── data/
│   ├── preprocessed/                   # Tracked in git (included in repo)
│   │   ├── ClickTimes/                 # Behavioral responses
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

### Analysis and Figure Generation

| Step | Script | Output |
|------|--------|--------|
| 2a | `plotFigure1_Behavioral.m` | `results/Figure1_Behavioral.{png,pdf}` |
| 2b | `plotFigure2_Neural.py` | `results/Figure2_Neural.{png,pdf}` (updated prototype: retains R² scalp in Panel C; adds Panel D for 1/f aperiodic exponent tracking S(x) with R²≈0.85, aligned under B) |
| 2c | `computePASCrossover.m` | `results/Figure3_crossover_stats.csv` |
| 2d | `plotFigure3_Perceptual.m` | `results/Figure3_Perceptual.{png,pdf}` |

## Dataset

32 dyads (64 participants) performed the perceptual crossing experiment. Each session comprised 18 trials of 60 seconds and 4 rest periods of 180 seconds. Simultaneous recordings include 64-channel EEG (1000 Hz), electrodermal activity, respiration, haptic feedback, behavioral responses (button presses), and subjective perceptual awareness ratings (PAS, 1-4 scale).

## Related Repositories

For the full analysis pipeline including respiration fits, see [icdl2026-embodied-anticipation](https://github.com/tom-froese/icdl2026-embodied-anticipation).

## License

MIT. See [LICENSE](LICENSE).
