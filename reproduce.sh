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
echo "  Reproducing: Froese (2026) PNAS Brief Report"
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

    echo "  [1g] Extracting parietal hemisphere data..."
    matlab -batch "cd('code/preprocessing'); extractParietalHemispheres"

    echo "  [1h] Computing per-channel sensitivity fits..."
    matlab -batch "cd('code/preprocessing'); computePerChannelFits"

    echo "  [1i] Computing PAS proportions (unsmoothed)..."
    matlab -batch "cd('code/preprocessing'); computePASProportions"

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
        data/preprocessed/EEG/parietal_hemisphere_data.mat \
        results/FigureS2_GSP_TopoMap_FreeTau_perchannel.csv \
        results/pas_unsmoothed_proportions.csv; do
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
echo "--- Step 2: Generating manuscript figures ---"

echo "  [2a] Figure 1: Behavioral and bodily evidence..."
matlab -batch "cd('code/analysis'); plotFigure1_Behavioral"

echo "  [2b] Figure 2: Neural evidence (Python/MNE)..."
python3 code/analysis/plotFigure2_Neural.py

echo "  [2c] PAS crossover statistics..."
matlab -batch "cd('code/analysis'); computePASCrossover"

echo "  [2d] Figure 3: Perceptual evidence..."
matlab -batch "cd('code/analysis'); plotFigure3_Perceptual"

echo ""
echo "============================================="
echo "  Done. Figures saved to results/"
echo "============================================="
ls -lh results/Figure*.pdf results/Figure*.png 2>/dev/null
