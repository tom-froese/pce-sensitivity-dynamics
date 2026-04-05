%% globalScalpPotential_stats.m
% =========================================================================
% Global Scalp Potential: Per-Participant P(1) Fits & Group Statistics
% =========================================================================
%
% PURPOSE:
%   Loads the participant-level time courses from globalScalpPotential.m
%   and performs a hierarchical statistical analysis:
%
%   (A) PRIMARY: Per-participant P(1) model fitting
%       - Smooth each participant's time course (5 s moving average)
%       - Invert (negate) the signal
%       - Fit P(1)(x) = A * lambda * x * exp(-lambda * x) + B
%         with lambda = e fixed, tau (onset lag) by grid search, A & B by OLS
%       - Extract per-participant: R², tau, peak time, amplitude A
%       - Group-level tests:
%         * One-sample t-test & Wilcoxon: R² > 0
%         * One-sample t-test & Wilcoxon: amplitude A ≠ 0
%         * Consistency of peak time (mean ± SD, CI)
%         * Cohen's d effect sizes
%
%   (B) ROBUSTNESS: Permutation test on grand-average fit
%       - Compute grand-average time course (median across participants)
%       - Fit P(1) → observed R²
%       - Permutation null: for each of nPerm iterations, circularly shift
%         each participant's time course by a random amount, recompute
%         grand average, refit P(1) → null R² distribution
%       - p-value = proportion of null R² ≥ observed R²
%
%   (C) ROI ANALYSIS: Repeat (A) for each of 7 anterior-posterior ROIs
%
%   (D) REST CONTROL: Fit P(1) to rest data to confirm no systematic
%       temporal structure (R² should be near zero)
%
% INPUT:
%   EEG/globalScalpPotential_data.mat (from globalScalpPotential.m)
%
% OUTPUT:
%   EEG/globalScalpPotential_stats.mat
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   April 2026
% =========================================================================

clearvars; close all; clc;

%% ========================================================================
%  1. LOAD DATA
%  ========================================================================

dataFile = fullfile(pwd, 'EEG', 'globalScalpPotential_data.mat');
if ~exist(dataFile, 'file')
    error('Data file not found: %s\nRun globalScalpPotential.m first.', dataFile);
end

fprintf('Loading %s ...\n', dataFile);
D = load(dataFile);

taskTC   = D.taskParticipantTC;   % [nPart x nSampTask]
restTC   = D.restParticipantTC;   % [nPartRest x nSampRest]
roiTC    = D.roiParticipantTC;    % {nROI x 1} cell
tTask    = D.tTask;
tRest    = D.tRest;
cfg      = D.cfg;
roi      = D.roi;
nROI     = D.nROI;
nPart    = D.nPart;

fprintf('  Participants: %d\n', nPart);
fprintf('  Task samples: %d (%.0f s at %d Hz)\n', ...
    D.nSampTask, cfg.trialDur_s, cfg.fsOut);
fprintf('\n');

%% ========================================================================
%  2. SMOOTHING CONFIGURATION
%  ========================================================================

smoothWin_s = 5.0;          % Smoothing window (seconds)
smoothK     = round(smoothWin_s * cfg.fsOut);  % kernel width in samples
kernel      = ones(1, smoothK) / smoothK;

% After conv(..., 'valid'), edges are truncated
halfK   = floor(smoothK / 2);

% Compute truncated time vector
taskTC_sm_test = conv(taskTC(1,:), kernel, 'valid');
nSampSm = length(taskTC_sm_test);
tTaskSm = tTask(halfK+1 : halfK+nSampSm);

lambda = exp(1);  % Fixed Poisson rate parameter

fprintf('Smoothing: %.1f s window (%d samples)\n', smoothWin_s, smoothK);
fprintf('Fitting window: %.1f – %.1f s (%d samples)\n', ...
    tTaskSm(1), tTaskSm(end), nSampSm);

%% ========================================================================
%  3. PER-PARTICIPANT P(1) FITS (GLOBAL)
%  ========================================================================

fprintf('\n--- Per-Participant P(1) Fits (Global Scalp Potential) ---\n');

partR2   = nan(nPart, 1);
partTau  = nan(nPart, 1);
partA    = nan(nPart, 1);
partB    = nan(nPart, 1);
partPeak = nan(nPart, 1);
partFits = nan(nPart, nSampSm);    % fitted curves
partSmTC = nan(nPart, nSampSm);    % smoothed inverted time courses

for pi = 1:nPart
    % Smooth
    sm = conv(taskTC(pi,:), kernel, 'valid');
    % Invert (negate) — we expect a negative-going potential that,
    % when inverted, shows the P(1) rise-peak-fall shape
    sm_inv = -sm;

    partSmTC(pi,:) = sm_inv;

    % Fit P(1)
    [yFit, R2, tau, A, B] = fitP1_gsp(tTaskSm, sm_inv, lambda);

    partR2(pi)      = R2;
    partTau(pi)     = tau;
    partA(pi)       = A;
    partB(pi)       = B;
    partFits(pi,:)  = yFit;

    % Peak time of the fitted curve
    [~, pkIdx] = max(yFit);
    partPeak(pi) = tTaskSm(pkIdx);
end

fprintf('  Per-participant R²: mean=%.3f, median=%.3f, SD=%.3f\n', ...
    mean(partR2), median(partR2), std(partR2));
fprintf('  Per-participant peak: mean=%.1f s, median=%.1f s, SD=%.1f s\n', ...
    mean(partPeak), median(partPeak), std(partPeak));

%% ========================================================================
%  4. GROUP-LEVEL STATISTICS (GLOBAL)
%  ========================================================================

fprintf('\n--- Group-Level Statistics ---\n');

df = nPart - 1;

% --- R² > 0 ---
[~, pR2_t, ~, sR2] = ttest(partR2, 0, 'Tail', 'right');
dz_R2 = mean(partR2) / std(partR2);
pR2_w = signrank(partR2, 0, 'tail', 'right');

fprintf('\n  R² > 0:\n');
fprintf('    t(%d) = %.3f, p = %.2e, d_z = %.3f\n', df, sR2.tstat, pR2_t, dz_R2);
fprintf('    Wilcoxon signed-rank: p = %.2e\n', pR2_w);
fprintf('    Mean R² = %.4f, 95%% CI = [%.4f, %.4f]\n', ...
    mean(partR2), ...
    mean(partR2) - 1.96*std(partR2)/sqrt(nPart), ...
    mean(partR2) + 1.96*std(partR2)/sqrt(nPart));

% --- Amplitude A ≠ 0 ---
[~, pA_t, ~, sA] = ttest(partA);
dz_A = mean(partA) / std(partA);
pA_w = signrank(partA);

fprintf('\n  Amplitude A ≠ 0:\n');
fprintf('    t(%d) = %.3f, p = %.2e, d_z = %.3f\n', df, sA.tstat, pA_t, dz_A);
fprintf('    Wilcoxon signed-rank: p = %.2e\n', pA_w);
fprintf('    Mean A = %.4f, SD = %.4f\n', mean(partA), std(partA));

% --- Peak time consistency ---
fprintf('\n  Peak time:\n');
fprintf('    Mean = %.2f s, SD = %.2f s\n', mean(partPeak), std(partPeak));
fprintf('    Median = %.2f s\n', median(partPeak));
fprintf('    95%% CI = [%.2f, %.2f] s\n', ...
    mean(partPeak) - 1.96*std(partPeak)/sqrt(nPart), ...
    mean(partPeak) + 1.96*std(partPeak)/sqrt(nPart));
fprintf('    IQR = [%.2f, %.2f] s\n', ...
    prctile(partPeak, 25), prctile(partPeak, 75));

% --- Tau (onset lag) ---
fprintf('\n  Onset lag (tau):\n');
fprintf('    Mean = %.2f s, SD = %.2f s\n', mean(partTau), std(partTau));
fprintf('    Median = %.2f s\n', median(partTau));

%% ========================================================================
%  5. GRAND-AVERAGE FIT (for comparison with existing analysis)
%  ========================================================================

fprintf('\n--- Grand-Average P(1) Fit ---\n');

grandMedian     = median(taskTC, 1, 'omitnan');
grandMedian_sm  = conv(grandMedian, kernel, 'valid');
grandMedian_inv = -grandMedian_sm;

[grandFit, grandR2, grandTau, grandA, grandB] = ...
    fitP1_gsp(tTaskSm, grandMedian_inv, lambda);

[~, grandPkIdx] = max(grandFit);
grandPeak = tTaskSm(grandPkIdx);

fprintf('  Grand-average R² = %.4f\n', grandR2);
fprintf('  Grand-average tau = %.2f s\n', grandTau);
fprintf('  Grand-average peak = %.2f s\n', grandPeak);
fprintf('  Grand-average A = %.4f, B = %.4f\n', grandA, grandB);

%% ========================================================================
%  6. PERMUTATION TEST ON GRAND-AVERAGE R²
%  ========================================================================

fprintf('\n--- Permutation Test (circular shift) ---\n');

nPerm = 5000;
nullR2 = nan(nPerm, 1);

rng(42, 'twister');  % reproducibility

for pi_perm = 1:nPerm
    % Circular-shift each participant's time course by a random amount
    shifted = nan(nPart, D.nSampTask);
    for si = 1:nPart
        shift_amt = randi(D.nSampTask);
        shifted(si,:) = circshift(taskTC(si,:), shift_amt);
    end

    % Recompute grand average
    nullGrand     = median(shifted, 1, 'omitnan');
    nullGrand_sm  = conv(nullGrand, kernel, 'valid');
    nullGrand_inv = -nullGrand_sm;

    % Fit P(1)
    [~, nullR2(pi_perm)] = fitP1_gsp(tTaskSm, nullGrand_inv, lambda);

    if mod(pi_perm, 1000) == 0
        fprintf('  Permutation %d/%d done\n', pi_perm, nPerm);
    end
end

pPerm = mean(nullR2 >= grandR2);
fprintf('  Observed grand-average R² = %.4f\n', grandR2);
fprintf('  Permutation p-value = %.4f (based on %d permutations)\n', pPerm, nPerm);
fprintf('  Null R² distribution: mean=%.4f, max=%.4f, 95th=%.4f\n', ...
    mean(nullR2), max(nullR2), prctile(nullR2, 95));

%% ========================================================================
%  7. REST CONTROL
%  ========================================================================

if ~isempty(restTC)
    fprintf('\n--- Rest Control ---\n');

    nPartRest = size(restTC, 1);

    % For rest, use same smoothing approach
    smoothK_rest = round(smoothWin_s * cfg.fsOut);
    kernel_rest  = ones(1, smoothK_rest) / smoothK_rest;
    halfK_rest   = floor(smoothK_rest / 2);

    restSm_test  = conv(restTC(1,:), kernel_rest, 'valid');
    nSampRest_sm = length(restSm_test);
    tRestSm      = tRest(halfK_rest+1 : halfK_rest+nSampRest_sm);

    restR2 = nan(nPartRest, 1);
    for pi = 1:nPartRest
        sm = conv(restTC(pi,:), kernel_rest, 'valid');
        sm_inv = -sm;
        [~, restR2(pi)] = fitP1_gsp(tRestSm, sm_inv, lambda);
    end

    fprintf('  Rest R²: mean=%.4f, median=%.4f, SD=%.4f\n', ...
        mean(restR2), median(restR2), std(restR2));

    % Compare task vs rest R²
    [~, pTaskRest, ~, sTaskRest] = ttest2(partR2, restR2);
    fprintf('  Task vs Rest R²: t = %.3f, p = %.2e\n', sTaskRest.tstat, pTaskRest);

    % Grand rest average
    grandRest     = median(restTC, 1, 'omitnan');
    grandRest_sm  = conv(grandRest, kernel_rest, 'valid');
    grandRest_inv = -grandRest_sm;
    [grandRestFit, grandRestR2] = fitP1_gsp(tRestSm, grandRest_inv, lambda);

    fprintf('  Grand rest R² = %.4f (should be near zero / much < task)\n', grandRestR2);
else
    fprintf('\n  No rest data available.\n');
    restR2 = []; grandRestR2 = NaN; tRestSm = []; grandRest_inv = []; grandRestFit = [];
end

%% ========================================================================
%  8. ROI ANALYSIS
%  ========================================================================

fprintf('\n--- ROI-Level Per-Participant P(1) Fits ---\n');
fprintf('%-18s  Mean_R2   SD_R2    t(%d)     p_t        d_z     Mean_Peak  SD_Peak\n', ...
    'ROI', df);
fprintf('%s\n', repmat('-', 1, 95));

roiPartR2   = cell(nROI, 1);
roiPartPeak = cell(nROI, 1);
roiPartA    = cell(nROI, 1);
roiPartTau  = cell(nROI, 1);
roiPartFits = cell(nROI, 1);
roiPartSmTC = cell(nROI, 1);

roiGrandR2   = nan(nROI, 1);
roiGrandTau  = nan(nROI, 1);
roiGrandPeak = nan(nROI, 1);
roiGrandA    = nan(nROI, 1);
roiGrandB    = nan(nROI, 1);
roiGrandFit  = nan(nROI, nSampSm);
roiGrandInv  = nan(nROI, nSampSm);

for ri = 1:nROI
    rTC = roiTC{ri};      % [nPart x nSampTask]
    nP  = size(rTC, 1);

    rR2   = nan(nP, 1);
    rPeak = nan(nP, 1);
    rA    = nan(nP, 1);
    rTau  = nan(nP, 1);
    rFits = nan(nP, nSampSm);
    rSmTC = nan(nP, nSampSm);

    for pi = 1:nP
        sm = conv(rTC(pi,:), kernel, 'valid');
        sm_inv = -sm;
        rSmTC(pi,:) = sm_inv;

        [yFit, R2, tau, A] = fitP1_gsp(tTaskSm, sm_inv, lambda);
        rR2(pi)     = R2;
        rTau(pi)    = tau;
        rA(pi)      = A;
        rFits(pi,:) = yFit;

        [~, pkIdx] = max(yFit);
        rPeak(pi)   = tTaskSm(pkIdx);
    end

    roiPartR2{ri}   = rR2;
    roiPartPeak{ri} = rPeak;
    roiPartA{ri}    = rA;
    roiPartTau{ri}  = rTau;
    roiPartFits{ri} = rFits;
    roiPartSmTC{ri} = rSmTC;

    % Group stats for this ROI
    [~, pval, ~, s] = ttest(rR2, 0, 'Tail', 'right');
    dz = mean(rR2) / std(rR2);

    fprintf('%-18s  %.4f    %.4f   %+6.3f    %-10s  %+.3f   %.1f s      %.1f s\n', ...
        roi.names{ri}, mean(rR2), std(rR2), s.tstat, ...
        formatP(pval), dz, mean(rPeak), std(rPeak));

    % Grand average for ROI
    rGrand     = median(rTC, 1, 'omitnan');
    rGrand_sm  = conv(rGrand, kernel, 'valid');
    rGrand_inv = -rGrand_sm;
    [rGFit, rGR2, rGTau, rGA, rGB] = fitP1_gsp(tTaskSm, rGrand_inv, lambda);

    roiGrandR2(ri)    = rGR2;
    roiGrandTau(ri)   = rGTau;
    roiGrandA(ri)     = rGA;
    roiGrandB(ri)     = rGB;
    roiGrandFit(ri,:) = rGFit;
    roiGrandInv(ri,:) = rGrand_inv;

    [~, pkIdx] = max(rGFit);
    roiGrandPeak(ri) = tTaskSm(pkIdx);
end

fprintf('\n--- Grand-Average ROI Fits ---\n');
fprintf('%-18s  R2_grand  Tau(s)  Peak(s)  A         B\n', 'ROI');
fprintf('%s\n', repmat('-', 1, 70));
for ri = 1:nROI
    fprintf('%-18s  %.4f    %.1f     %.1f      %+.4f    %+.4f\n', ...
        roi.names{ri}, roiGrandR2(ri), roiGrandTau(ri), ...
        roiGrandPeak(ri), roiGrandA(ri), roiGrandB(ri));
end

%% ========================================================================
%  9. SAVE ALL RESULTS
%  ========================================================================

outputFile = fullfile(pwd, 'EEG', 'globalScalpPotential_stats.mat');
save(outputFile, ...
    'partR2', 'partTau', 'partA', 'partB', 'partPeak', ...
    'partFits', 'partSmTC', ...
    'grandMedian_inv', 'grandFit', 'grandR2', 'grandTau', 'grandPeak', 'grandA', 'grandB', ...
    'nullR2', 'pPerm', 'nPerm', ...
    'restR2', 'grandRestR2', 'tRestSm', 'grandRest_inv', 'grandRestFit', ...
    'roiPartR2', 'roiPartPeak', 'roiPartA', 'roiPartTau', ...
    'roiPartFits', 'roiPartSmTC', ...
    'roiGrandR2', 'roiGrandTau', 'roiGrandPeak', 'roiGrandA', 'roiGrandB', ...
    'roiGrandFit', 'roiGrandInv', ...
    'tTaskSm', 'smoothWin_s', 'lambda', ...
    'cfg', 'roi', 'nROI', 'nPart', ...
    '-v7.3');

fprintf('\nSaved: %s\n', outputFile);

%% ========================================================================
%  10. FINAL SUMMARY
%  ========================================================================

fprintf('\n==========================================================\n');
fprintf('  GLOBAL SCALP POTENTIAL — STATISTICAL SUMMARY\n');
fprintf('==========================================================\n');
fprintf('  N = %d participants\n', nPart);
fprintf('\n  GLOBAL (per-participant P(1) fits):\n');
fprintf('    R²:        M=%.3f, SD=%.3f, t(%d)=%.2f, p=%.2e, d_z=%.3f\n', ...
    mean(partR2), std(partR2), df, sR2.tstat, pR2_t, dz_R2);
fprintf('    Amplitude: M=%.4f, SD=%.4f, t(%d)=%.2f, p=%.2e\n', ...
    mean(partA), std(partA), df, sA.tstat, pA_t);
fprintf('    Peak time: M=%.1f s, SD=%.1f s\n', mean(partPeak), std(partPeak));
fprintf('    Onset lag: M=%.1f s, SD=%.1f s\n', mean(partTau), std(partTau));
fprintf('\n  GRAND-AVERAGE:\n');
fprintf('    R² = %.4f, tau = %.1f s, peak = %.1f s\n', grandR2, grandTau, grandPeak);
fprintf('    Permutation p = %.4f (%d perms)\n', pPerm, nPerm);
if ~isempty(restR2)
    fprintf('\n  REST CONTROL:\n');
    fprintf('    Grand rest R² = %.4f\n', grandRestR2);
    fprintf('    Per-participant rest R²: M=%.3f, SD=%.3f\n', mean(restR2), std(restR2));
end
fprintf('\n  ROI SUMMARY (grand-average R²):\n');
for ri = 1:nROI
    fprintf('    %-18s  R² = %.4f, peak = %.1f s\n', ...
        roi.names{ri}, roiGrandR2(ri), roiGrandPeak(ri));
end
fprintf('==========================================================\n');
fprintf('Done.\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function [yFit, bestR2, bestTau, bestA, bestB] = fitP1_gsp(tVec, yData, lambda)
%FITP1_GSP  Fit Poisson first-order term P(1) to a time series.
%
%   Model: y(t) = A * lambda * x * exp(-lambda * x) + B
%          where x = (t_rel - tau) / (T_rel - tau),
%          t_rel = tVec - tVec(1), T_rel = tVec(end) - tVec(1).
%   lambda is fixed (= e). A, B estimated by OLS at each candidate tau.
%   tau is found by grid search maximizing R².

    t0   = tVec(1);
    tRel = tVec(:) - t0;
    Trel = tRel(end);

    tauGrid = 0:0.1:min(15, Trel/4);
    nTau    = length(tauGrid);

    y     = yData(:);
    SStot = sum((y - mean(y)).^2);
    if SStot == 0
        yFit = zeros(size(yData));
        bestR2 = 0; bestTau = 0; bestA = 0; bestB = 0;
        return;
    end

    bestR2 = -Inf;
    yFit   = zeros(size(yData));
    bestTau = 0; bestA = 0; bestB = 0;

    for k = 1:nTau
        tau  = tauGrid(k);
        Teff = Trel - tau;
        if Teff <= 0, continue; end

        x = (tRel - tau) / Teff;
        x(x < 0) = 0;

        shape = lambda .* x .* exp(-lambda .* x);
        X     = [shape, ones(size(shape))];
        beta  = X \ y;
        yhat  = X * beta;
        SSres = sum((y - yhat).^2);
        R2    = 1 - SSres / SStot;

        if R2 > bestR2
            bestR2  = R2;
            bestTau = tau;
            bestA   = beta(1);
            bestB   = beta(2);
            yFit    = yhat(:)';
        end
    end
end


function s = formatP(p)
%FORMATP  Format a p-value for console display.
    if p < 0.001
        s = sprintf('%.2e', p);
    else
        s = sprintf('%.4f', p);
    end
end
