#!/bin/bash
# reproduce.sh — Master script to reproduce all results from Froese (2026)
# "Clock-free optimal stopping in decision-making"
#
# Usage:
#   ./reproduce.sh              # from preprocessed data (default)
#   ./reproduce.sh --from-raw   # from raw data (requires data/raw/)
#
# Prerequisites:
#   MATLAB (R2022a or later) with Statistics and Optimization toolboxes
#   Python 3.9+ with MNE-Python, matplotlib, numpy, pandas, scipy
#
# Data:
#   Preprocessed: https://doi.org/10.5281/zenodo.19425014
#   Raw:          https://osf.io/47n3p

set -euo pipefail
cd "$(dirname "$0")"
ROOT="$(pwd)"

FROM_RAW=false
if [[ "${1:-}" == "--from-raw" ]]; then
    FROM_RAW=true
fi

echo "============================================="
echo "  Reproducing: Froese (2026) PNAS Nexus Research Report"
echo "  Mode: $(if $FROM_RAW; then echo 'from raw data'; else echo 'from preprocessed data'; fi)"
echo "============================================="
echo ""

# ------------------------------------------------------------------
# STEP 1: Preprocessing (only needed if --from-raw)
# ------------------------------------------------------------------
if $FROM_RAW; then
    echo "--- Step 1: Preprocessing raw data ---"

    if [ ! -d "data/raw" ]; then
        echo "ERROR: data/raw/ not found. Download from https://osf.io/47n3p"
        exit 1
    fi

    echo "  [1a] Preprocessing click response times..."
    matlab -batch "cd('code/preprocessing'); preprocessClicks"

    echo "  [1b] Preprocessing haptic feedback..."
    matlab -batch "cd('code/preprocessing'); preprocessHaptics"

    echo "  [1c] Preprocessing EDA..."
    matlab -batch "cd('code/preprocessing'); preprocessEDA"

    echo "  [1d] Preprocessing global scalp potential..."
    matlab -batch "cd('code/preprocessing'); preprocessGSP"

    echo "  [1e] Computing GSP statistics..."
    matlab -batch "cd('code/preprocessing'); computeGSPStats"

    echo "  [1f] Extracting all 64 EEG channels..."
    matlab -batch "cd('code/preprocessing'); extractAllChannels"

    echo "  [1g] Computing per-channel sensitivity fits..."
    matlab -batch "cd('code/preprocessing'); computePerChannelFits"

    echo "  [1i] Computing PAS proportions (unsmoothed)..."
    matlab -batch "cd('code/preprocessing'); computePASProportions"

    echo "  [1j] Computing PAS crossover statistics..."
    matlab -batch "cd('code/analysis'); computePASCrossover"

    echo "  [1k] Preprocessing 250 Hz cleaned EEG for the aperiodic-exponent analysis..."
    python3 code/preprocessing/preprocessEEGForExponent.py --all

    echo "  Step 1 complete."
    echo ""
else
    echo "--- Step 1: Skipped (using preprocessed data) ---"
    echo "  Checking preprocessed data exists..."

    MISSING=false
    for f in \
        data/preprocessed/ClickTimes/ClickResponseTimes.csv \
        data/preprocessed/Haptics/HapticFeedback.csv.gz \
        data/preprocessed/EDA/EDA_Task_Preprocessed.csv \
        data/preprocessed/EEG/globalScalpPotential_data.mat \
        data/preprocessed/EEG/globalScalpPotential_stats.mat \
        data/preprocessed/EEG/allchannel_data.mat \
        data/preprocessed/EEG/pce01/pce01_P1_task-raw.fif \
        results/Figure2_perchannel_fits.csv \
        results/Figure3_pas_proportions.csv \
        results/Figure3_crossover_stats.csv; do
        if [ ! -f "$f" ]; then
            echo "  MISSING: $f"
            MISSING=true
        fi
    done

    if $MISSING; then
        echo ""
        echo "  Some preprocessed files are missing. Either:"
        echo "    1. Download preprocessed data from https://doi.org/10.5281/zenodo.19425014"
        echo "    2. Run ./reproduce.sh --from-raw (requires raw data)"
        exit 1
    fi
    echo "  All preprocessed data present."
    echo ""
fi

# ------------------------------------------------------------------
# STEP 2: Generate figures
# ------------------------------------------------------------------
# This script generates the four code-derived figures. Their internal file
# names (Figure1/2/3 + FigConcept) map to the manuscript as:
#   Figure1_Behavioral   -> manuscript Fig 3 (behavioral/bodily)
#   FigConcept_Sensitivity -> manuscript Fig 2 (R(x), S(x) concept curves)
#   Figure2_Neural       -> manuscript Fig 4 (neural: GSP, topography, 1/f)
#   Figure3_Perceptual   -> manuscript Fig 5 (perceptual: PAS)
# Manuscript Fig 1 (the apparatus/task/protocol) is a composite of curated
# images, not generated here.
echo "--- Step 2: Generating manuscript figures ---"

echo "  [2a] Figure 1: Behavioral and bodily evidence..."
matlab -batch "cd('code/analysis'); plotFigure1_Behavioral"

echo "  [2b] Aperiodic exponent per (participant x within-trial bin) -- FOOOF..."
python3 code/analysis/computeAperiodicExponent.py

echo "  [2c] Aperiodic exponent S(x) fit + bootstrap peak CI..."
python3 code/analysis/fitExponentSensitivity.py

echo "  [2d] Figure 2: Neural evidence (Python/MNE)..."
python3 code/analysis/plotFigure2_Neural.py

echo "  [2e] Figure 3: Perceptual evidence..."
matlab -batch "cd('code/analysis'); plotFigure3_Perceptual"

echo "  [2f] Sensitivity-framework concept figure (R(x), S(x) curves)..."
python3 code/analysis/plotFigConcept_Sensitivity.py

echo ""
echo "============================================="
echo "  Done. Figures saved to results/"
echo "============================================="
ls -lh results/Figure*.pdf results/Figure*.png results/FigConcept*.pdf results/FigConcept*.png 2>/dev/null
