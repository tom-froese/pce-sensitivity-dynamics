# Sensitivity Dynamics in Perceptual Crossing

Code and data for:

> **Clock-free optimal stopping in social cognition: Stochastic sensitivity dynamics predict behavioral, neural, and perceptual transition timing**
>
> Tom Froese (2026). *Under review.*

## Overview

This repository contains all code needed to reproduce the three figures reported in the paper. We show that participants in a perceptual crossing experiment (PCE) exhibit sensitivity dynamics consistent with the derivative of an exponential decay with rate parameter lambda = *e*:

    S(x) = |dP(0)/dlambda| = x * exp(-e * x)

This sensitivity function peaks at *x* = 1/*e*, the optimal stopping point of the 1/*e*-law, and is reflected in behavioral responses (click times), neural activity (global scalp potential), and perceptual awareness (PAS ratings).

## Quick Start

From MATLAB (R2023b or later), run each figure script from `code/analysis/`:

```matlab
cd code/analysis
plotFigure1_Behavioral   % Figure 1: Click times + haptic contact
plotFigure2_Neural       % Figure 2: Global scalp potential + ROI analysis
plotFigure3_Perceptual   % Figure 3: Sensitivity schematic + PAS clarity
```

All three scripts use relative paths and output 600 dpi PNG files to `results/`.

## Data

Most data files are included in this repository. The per-channel EEG time series (~1.6 GB) are archived separately on Zenodo due to their size.

**Included in this repository:**

| Directory | Contents | Size |
|-----------|----------|------|
| `data/ClickTimes/` | Click response times (CSV + JSON sidecar) | 32 KB |
| `data/Haptics/` | Haptic feedback time series (gzipped CSV) | 18 MB |
| `data/PAS/` | PAS sensitivity and crossover analysis (MAT) | 40 KB |
| `data/EEG/` | Global scalp potential data and statistics (MAT) | 7 MB |

**On Zenodo (required only for re-running preprocessing):**

| Contents | Size |
|----------|------|
| 62 per-channel EEG task CSVs (64 ch, 10 Hz, 60 s/trial) | 975 MB |
| 62 per-channel EEG rest CSVs (64 ch, 10 Hz, 180 s/rest) | 651 MB |

Zenodo DOI: *[to be added upon upload]*

To use the Zenodo data, download and extract into `data/EEG/` so that the per-channel CSV files sit alongside the existing MAT files.

## Toolbox Requirements

- Signal Processing Toolbox (filtering, spectral estimation)
- Statistics and Machine Learning Toolbox (logistic regression, permutation tests)
- Optimization Toolbox (curve fitting in Figure 1 Panel B)

No third-party packages are required.

## Repository Structure

```
pce-sensitivity-dynamics/
├── code/
│   ├── analysis/                    # Figure generation scripts
│   │   ├── plotFigure1_Behavioral.m # Fig 1: clicks + haptics
│   │   ├── plotFigure2_Neural.m     # Fig 2: GSP + ROI + gradient
│   │   └── plotFigure3_Perceptual.m # Fig 3: sensitivity + PAS
│   └── preprocessing/               # Data extraction pipelines
│       ├── preprocessClicks.m       # Click response times from raw trials
│       ├── preprocessHaptics.m      # Haptic feedback from raw trials
│       ├── preprocessPAS.m          # PAS ratings from questionnaires
│       ├── preprocessGSP.m          # Global scalp potential from raw EEG
│       ├── computeGSPStats.m        # Hierarchical sensitivity model fits
│       ├── computePASSensitivity.m  # PAS early vs. late analysis
│       └── computePASCrossover.m    # PAS 4/3 crossover + ROI mapping
├── data/
│   ├── ClickTimes/                  # Behavioral responses
│   ├── EEG/                         # Neural data (MAT + Zenodo CSVs)
│   ├── Haptics/                     # Haptic feedback time series
│   └── PAS/                         # Perceptual awareness analysis
├── results/                         # Publication figures (600 dpi PNG)
├── README.md
├── LICENSE                          # MIT
└── CITATION.cff
```

## Preprocessing

The figure scripts load precomputed data files that are included in this repository. The preprocessing scripts document how these files were generated from raw experimental data and are provided for transparency. To re-run preprocessing from scratch, you need access to the raw PCE dataset (available from the corresponding author upon request).

Execution order:

1. `preprocessClicks.m` and `preprocessHaptics.m` (independent, from raw trial CSVs)
2. `preprocessPAS.m` (from raw questionnaire CSVs)
3. `preprocessGSP.m` (from raw EEG MAT files; produces `globalScalpPotential_data.mat`)
4. `computeGSPStats.m` (from step 3 output; produces `globalScalpPotential_stats.mat`)
5. `computePASSensitivity.m` and `computePASCrossover.m` (from steps 1, 2, 4)

## Dataset

32 dyads (64 participants) performed the perceptual crossing experiment. Each session comprised 18 trials of 60 seconds and 4 rest periods of 180 seconds. Simultaneous recordings include 64-channel EEG (1000 Hz), electrodermal activity, respiration, haptic feedback, behavioral responses (button presses), and subjective perceptual awareness ratings (PAS, 1-4 scale).

## Related Repositories

For the full analysis pipeline including EDA and respiration fits, see [icdl2026-embodied-anticipation](https://github.com/tom-froese/icdl2026-embodied-anticipation).

## License

MIT. See [LICENSE](LICENSE).
