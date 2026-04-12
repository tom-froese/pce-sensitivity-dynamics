%% globalScalpPotential_stats.m
% =========================================================================
% Global Scalp Potential: Per-Participant Sensitivity Fits & Group Stats
% =========================================================================
%
% PURPOSE:
%   Loads the participant-level time courses from globalScalpPotential.m
%   and performs a hierarchical statistical analysis:
%
%   (A) PRIMARY: Per-participant sensitivity model fitting
%       - Smooth each participant's time course (5 s moving average)
%       - Fit S(x) = A * x * exp(-e * x) + B  directly (no inversion)
%         with lambda = e fixed, tau (onset lag) by grid search, A & B by OLS
%       - A is expected to be negative (the GSP has a natural trough)
%       - Extract per-participant: R², tau, trough time, amplitude A
%       - Group-level tests:
%         * One-sample t-test & Wilcoxon: R² > 0
%         * One-sample t-test & Wilcoxon: amplitude A ≠ 0
%         * Consistency of trough time (mean ± SD, CI)
%         * Cohen's d effect sizes
%
%   (B) ROBUSTNESS: Permutation test on grand-average fit
%       - Compute grand-average time course (median across participants)
%       - Fit S(x) → observed R²
%       - Permutation null: for each of nPerm iterations, circularly shift
%         each participant's time course by a random amount, recompute
%         grand average, refit S(x) → null R² distribution
%       - p-value = proportion of null R² ≥ observed R²
%
%   (C) ROI ANALYSIS: Repeat (A) for each of 7 anterior-posterior ROIs
%
%   (D) REST CONTROL: Fit S(x) to rest data to confirm no systematic
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

scriptDir = fileparts(mfilename('fullpath'));
eegDir   = fullfile(scriptDir, '..', '..', 'data', 'preprocessed', 'EEG');
dataFile = fullfile(eegDir, 'globalScalpPotential_data.mat');
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

lambda = exp(1);  % Fixed rate parameter (sensitivity framework)

fprintf('Smoothing: %.1f s window (%d samples)\n', smoothWin_s, smoothK);
fprintf('Fitting window: %.1f – %.1f s (%d samples)\n', ...
    tTaskSm(1), tTaskSm(end), nSampSm);

%% ========================================================================
%  3. PER-PARTICIPANT P(1) FITS (GLOBAL)
%  ========================================================================

fprintf('\n--- Per-Participant Sensitivity Fits (Global Scalp Potential) ---\n');

partR2     = nan(nPart, 1);
partTau    = nan(nPart, 1);
partA      = nan(nPart, 1);
partB      = nan(nPart, 1);
partTrough = nan(nPart, 1);
partFits   = nan(nPart, nSampSm);    % fitted curves
partSmTC   = nan(nPart, nSampSm);    % smoothed time courses (raw, no inversion)

for pi = 1:nPart
    % Smooth (no inversion — fit sensitivity directly to raw GSP)
    sm = conv(taskTC(pi,:), kernel, 'valid');
    partSmTC(pi,:) = sm;

    % Fit S(x) = A*x*exp(-e*x) + B  (A expected negative → trough)
    [yFit, R2, tau, A, B] = fitSensitivity_gsp(tTaskSm, sm, lambda);

    partR2(pi)      = R2;
    partTau(pi)     = tau;
    partA(pi)       = A;
    partB(pi)       = B;
    partFits(pi,:)  = yFit;

    % Trough time of the fitted curve (sensitivity minimum)
    [~, trIdx] = min(yFit);
    partTrough(pi) = tTaskSm(trIdx);
end

fprintf('  Per-participant R²: mean=%.3f, median=%.3f, SD=%.3f\n', ...
    mean(partR2), median(partR2), std(partR2));
fprintf('  Per-participant trough: mean=%.1f s, median=%.1f s, SD=%.1f s\n', ...
    mean(partTrough), median(partTrough), std(partTrough));

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

% --- Trough time consistency ---
fprintf('\n  Trough time:\n');
fprintf('    Mean = %.2f s, SD = %.2f s\n', mean(partTrough), std(partTrough));
fprintf('    Median = %.2f s\n', median(partTrough));
fprintf('    95%% CI = [%.2f, %.2f] s\n', ...
    mean(partTrough) - 1.96*std(partTrough)/sqrt(nPart), ...
    mean(partTrough) + 1.96*std(partTrough)/sqrt(nPart));
fprintf('    IQR = [%.2f, %.2f] s\n', ...
    prctile(partTrough, 25), prctile(partTrough, 75));

% --- Tau (onset lag) ---
fprintf('\n  Onset lag (tau):\n');
fprintf('    Mean = %.2f s, SD = %.2f s\n', mean(partTau), std(partTau));
fprintf('    Median = %.2f s\n', median(partTau));

%% ========================================================================
%  5. GRAND-AVERAGE FIT (for comparison with existing analysis)
%  ========================================================================

fprintf('\n--- Grand-Average Sensitivity Fit ---\n');

grandMedian    = median(taskTC, 1, 'omitnan');
grandMedianSm  = conv(grandMedian, kernel, 'valid');

[grandFit, grandR2, grandTau, grandA, grandB] = ...
    fitSensitivity_gsp(tTaskSm, grandMedianSm, lambda);

[~, grandTrIdx] = min(grandFit);
grandTrough = tTaskSm(grandTrIdx);

fprintf('  Grand-average R² = %.4f\n', grandR2);
fprintf('  Grand-average tau = %.2f s\n', grandTau);
fprintf('  Grand-average trough = %.2f s\n', grandTrough);
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
    nullGrand    = median(shifted, 1, 'omitnan');
    nullGrand_sm = conv(nullGrand, kernel, 'valid');

    % Fit S(x)
    [~, nullR2(pi_perm)] = fitSensitivity_gsp(tTaskSm, nullGrand_sm, lambda);

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
        [~, restR2(pi)] = fitSensitivity_gsp(tRestSm, sm, lambda);
    end

    fprintf('  Rest R²: mean=%.4f, median=%.4f, SD=%.4f\n', ...
        mean(restR2), median(restR2), std(restR2));

    % Compare task vs rest R²
    [~, pTaskRest, ~, sTaskRest] = ttest2(partR2, restR2);
    fprintf('  Task vs Rest R²: t = %.3f, p = %.2e\n', sTaskRest.tstat, pTaskRest);

    % Grand rest average
    grandRest    = median(restTC, 1, 'omitnan');
    grandRestSm  = conv(grandRest, kernel_rest, 'valid');
    [grandRestFit, grandRestR2] = fitSensitivity_gsp(tRestSm, grandRestSm, lambda);

    fprintf('  Grand rest R² = %.4f (should be near zero / much < task)\n', grandRestR2);
else
    fprintf('\n  No rest data available.\n');
    restR2 = []; grandRestR2 = NaN; tRestSm = []; grandRestSm = []; grandRestFit = [];
end

%% ========================================================================
%  8. ROI ANALYSIS
%  ========================================================================

fprintf('\n--- ROI-Level Per-Participant Sensitivity Fits ---\n');
fprintf('%-18s  Mean_R2   SD_R2    t(%d)     p_t        d_z     Mean_Tr    SD_Tr\n', ...
    'ROI', df);
fprintf('%s\n', repmat('-', 1, 95));

roiPartR2     = cell(nROI, 1);
roiPartTrough = cell(nROI, 1);
roiPartA      = cell(nROI, 1);
roiPartTau    = cell(nROI, 1);
roiPartFits   = cell(nROI, 1);
roiPartSmTC   = cell(nROI, 1);

roiGrandR2     = nan(nROI, 1);
roiGrandTau    = nan(nROI, 1);
roiGrandTrough = nan(nROI, 1);
roiGrandA      = nan(nROI, 1);
roiGrandB      = nan(nROI, 1);
roiGrandFit    = nan(nROI, nSampSm);
roiGrandSm     = nan(nROI, nSampSm);

for ri = 1:nROI
    rTC = roiTC{ri};      % [nPart x nSampTask]
    nP  = size(rTC, 1);

    rR2     = nan(nP, 1);
    rTrough = nan(nP, 1);
    rA      = nan(nP, 1);
    rTau    = nan(nP, 1);
    rFits   = nan(nP, nSampSm);
    rSmTC   = nan(nP, nSampSm);

    for pi = 1:nP
        sm = conv(rTC(pi,:), kernel, 'valid');
        rSmTC(pi,:) = sm;    % raw smoothed (no inversion)

        [yFit, R2, tau, A] = fitSensitivity_gsp(tTaskSm, sm, lambda);
        rR2(pi)     = R2;
        rTau(pi)    = tau;
        rA(pi)      = A;
        rFits(pi,:) = yFit;

        [~, trIdx] = min(yFit);   % trough (A < 0)
        rTrough(pi) = tTaskSm(trIdx);
    end

    roiPartR2{ri}     = rR2;
    roiPartTrough{ri} = rTrough;
    roiPartA{ri}      = rA;
    roiPartTau{ri}    = rTau;
    roiPartFits{ri}   = rFits;
    roiPartSmTC{ri}   = rSmTC;

    % Group stats for this ROI
    [~, pval, ~, s] = ttest(rR2, 0, 'Tail', 'right');
    dz = mean(rR2) / std(rR2);

    fprintf('%-18s  %.4f    %.4f   %+6.3f    %-10s  %+.3f   %.1f s      %.1f s\n', ...
        roi.names{ri}, mean(rR2), std(rR2), s.tstat, ...
        formatP(pval), dz, mean(rTrough), std(rTrough));

    % Grand average for ROI (no inversion)
    rGrand    = median(rTC, 1, 'omitnan');
    rGrand_sm = conv(rGrand, kernel, 'valid');
    [rGFit, rGR2, rGTau, rGA, rGB] = fitSensitivity_gsp(tTaskSm, rGrand_sm, lambda);

    roiGrandR2(ri)      = rGR2;
    roiGrandTau(ri)     = rGTau;
    roiGrandA(ri)       = rGA;
    roiGrandB(ri)       = rGB;
    roiGrandFit(ri,:)   = rGFit;
    roiGrandSm(ri,:)    = rGrand_sm;

    [~, trIdx] = min(rGFit);   % trough (A < 0)
    roiGrandTrough(ri) = tTaskSm(trIdx);
end

fprintf('\n--- Grand-Average ROI Fits ---\n');
fprintf('%-18s  R2_grand  Tau(s)  Trough(s)  A         B\n', 'ROI');
fprintf('%s\n', repmat('-', 1, 70));
for ri = 1:nROI
    fprintf('%-18s  %.4f    %.1f     %.1f        %+.4f    %+.4f\n', ...
        roi.names{ri}, roiGrandR2(ri), roiGrandTau(ri), ...
        roiGrandTrough(ri), roiGrandA(ri), roiGrandB(ri));
end

%% ========================================================================
%  9. SAVE ALL RESULTS
%  ========================================================================

outputFile = fullfile(eegDir, 'globalScalpPotential_stats.mat');
save(outputFile, ...
    'partR2', 'partTau', 'partA', 'partB', 'partTrough', ...
    'partFits', 'partSmTC', ...
    'grandMedianSm', 'grandFit', 'grandR2', 'grandTau', 'grandTrough', 'grandA', 'grandB', ...
    'nullR2', 'pPerm', 'nPerm', ...
    'restR2', 'grandRestR2', 'tRestSm', 'grandRestSm', 'grandRestFit', ...
    'roiPartR2', 'roiPartTrough', 'roiPartA', 'roiPartTau', ...
    'roiPartFits', 'roiPartSmTC', ...
    'roiGrandR2', 'roiGrandTau', 'roiGrandTrough', 'roiGrandA', 'roiGrandB', ...
    'roiGrandFit', 'roiGrandSm', ...
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
fprintf('\n  GLOBAL (per-participant sensitivity fits):\n');
fprintf('    R²:        M=%.3f, SD=%.3f, t(%d)=%.2f, p=%.2e, d_z=%.3f\n', ...
    mean(partR2), std(partR2), df, sR2.tstat, pR2_t, dz_R2);
fprintf('    Amplitude: M=%.4f, SD=%.4f, t(%d)=%.2f, p=%.2e\n', ...
    mean(partA), std(partA), df, sA.tstat, pA_t);
fprintf('    Trough:    M=%.1f s, SD=%.1f s\n', mean(partTrough), std(partTrough));
fprintf('    Onset lag: M=%.1f s, SD=%.1f s\n', mean(partTau), std(partTau));
fprintf('\n  GRAND-AVERAGE:\n');
fprintf('    R² = %.4f, tau = %.1f s, trough = %.1f s\n', grandR2, grandTau, grandTrough);
fprintf('    Permutation p = %.4f (%d perms)\n', pPerm, nPerm);
if ~isempty(restR2)
    fprintf('\n  REST CONTROL:\n');
    fprintf('    Grand rest R² = %.4f\n', grandRestR2);
    fprintf('    Per-participant rest R²: M=%.3f, SD=%.3f\n', mean(restR2), std(restR2));
end
fprintf('\n  ROI SUMMARY (grand-average R²):\n');
for ri = 1:nROI
    fprintf('    %-18s  R² = %.4f, trough = %.1f s\n', ...
        roi.names{ri}, roiGrandR2(ri), roiGrandTrough(ri));
end
fprintf('==========================================================\n');
fprintf('Done.\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function [yFit, bestR2, bestTau, bestA, bestB] = fitSensitivity_gsp(tVec, yData, lambda)
%FITSENSITIVITY_GSP  Fit sensitivity function S(x) to a GSP time series.
%
%   Model: y(t) = A * lambda * x * exp(-lambda * x) + B
%          where x = (t_rel - tau) / (T_rel - tau),
%          t_rel = tVec - tVec(1), T_rel = tVec(end) - tVec(1).
%   lambda is fixed (= e). A, B estimated by OLS at each candidate tau.
%   A is expected negative for raw GSP (trough at x = 1/lambda).
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
