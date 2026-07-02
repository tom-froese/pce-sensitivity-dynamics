%% computePASCrossover.m
% =========================================================================
% PAS 4 vs PAS 3 Logistic Crossover Analysis
% =========================================================================
%
% Computes the time at which PAS 4 ("clear experience") ceases to
% dominate over PAS 3 ("almost clear") using binomial logistic
% regression on click times.
%
% Methods:
%   1. Trial-level logistic regression (all PAS 3 & 4 trials pooled)
%   2. Within-participant logistic regression (per participant,
%      then median/mean of individual crossover times)
%   3. Early vs late comparison (split at predicted sensitivity peak)
%
% The crossover time is tau-independent (purely data-driven).
% The early/late split uses locked tau = 3.90 s (Figures 1-2).
%
% INPUT:
%   data/preprocessed/ClickTimes/ClickResponseTimes.csv
%   data/raw/Behavior/pce*/questionnaires/pair_*_PAS_confidence_absence.csv
%
% OUTPUT:
%   results/Figure3_crossover_stats.csv
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   May 2026
% =========================================================================

scriptDir = fileparts(mfilename('fullpath'));
ROOT      = fullfile(scriptDir, '..', '..');
dataDir   = fullfile(ROOT, 'data');
resDir    = fullfile(ROOT, 'results');

TAU_LOCKED = 3.90;
T          = 60;
lambda     = exp(1);
tPeak      = TAU_LOCKED + (T - TAU_LOCKED) / lambda;

fprintf('Locked tau = %.1f s,  tPeak = %.1f s\n\n', TAU_LOCKED, tPeak);

%% ========================================================================
%  1. LOAD CLICK TIMES
%  ========================================================================

CT = readtable(fullfile(dataDir, 'preprocessed', 'ClickTimes', ...
    'ClickResponseTimes.csv'));
CT.DyadNum = floor(CT.DyadID / 1e6);
CT = CT(CT.DyadNum ~= 31, :);

%% ========================================================================
%  2. LOAD PAS FROM RAW QUESTIONNAIRES
%  ========================================================================

rawBehDir = fullfile(dataDir, 'raw', 'Behavior');
allDyads  = dir(rawBehDir);
allDyads  = allDyads([allDyads.isdir]);

pasRows = [];
for i = 1:length(allDyads)
    tok = regexp(allDyads(i).name, '^pce(\d{2})\d{6}$', 'tokens');
    if isempty(tok), continue; end
    dN = str2double(tok{1}{1});
    if dN == 31, continue; end
    for p = [1 2]
        fn = sprintf('pair_%02d_P%d_PAS_confidence_absence.csv', dN, p);
        fp = fullfile(rawBehDir, allDyads(i).name, 'questionnaires', fn);
        if ~exist(fp, 'file'), continue; end
        tbl = readtable(fp, 'TextType', 'string');
        pasMask = tbl.question_id == "click_presence";
        if sum(pasMask) == 0, continue; end
        trials_0 = tbl.trial_id(pasMask);
        pasVals  = tbl.answer(pasMask);
        for ti = 1:length(trials_0)
            pasRows = [pasRows; dN, p, trials_0(ti)+1, pasVals(ti)]; %#ok<AGROW>
        end
    end
end
PAS = array2table(pasRows, 'VariableNames', ...
    {'DyadNum','ParticipantID','TrialNum','PAS'});

fprintf('Loaded %d click entries, %d PAS entries\n', height(CT), height(PAS));

%% ========================================================================
%  3. MERGE
%  ========================================================================

CT.Key  = CT.DyadNum * 1000 + CT.ParticipantID * 100 + CT.TrialNum;
PAS.Key = PAS.DyadNum * 1000 + PAS.ParticipantID * 100 + PAS.TrialNum;
[~, iCT, iPAS] = intersect(CT.Key, PAS.Key);

merged = table();
merged.DyadNum       = CT.DyadNum(iCT);
merged.ParticipantID = CT.ParticipantID(iCT);
merged.TrialNum      = CT.TrialNum(iCT);
merged.ClickTime     = CT.ClickTime_s(iCT);
merged.Clicked       = CT.Clicked(iCT);
merged.PAS           = PAS.PAS(iPAS);
merged = merged(merged.Clicked == 1 & ~isnan(merged.ClickTime), :);
merged.PartID = merged.DyadNum * 10 + merged.ParticipantID;

fprintf('Merged: %d click-PAS trials\n\n', height(merged));

%% ========================================================================
%  4. TRIAL-LEVEL LOGISTIC REGRESSION
%  ========================================================================

m34 = merged(merged.PAS == 3 | merged.PAS == 4, :);
binaryPAS = double(m34.PAS == 4);

mdl = fitglm(m34.ClickTime, binaryPAS, 'Distribution', 'binomial');
b0 = mdl.Coefficients.Estimate(1);
b1 = mdl.Coefficients.Estimate(2);
p_b1 = mdl.Coefficients.pValue(2);
tcross_trial = -b0 / b1;

fprintf('=== TRIAL-LEVEL LOGISTIC (n = %d) ===\n', height(m34));
fprintf('  beta0 = %.4f, beta1 = %.4f\n', b0, b1);
fprintf('  Crossover = %.1f s\n', tcross_trial);
fprintf('  p(beta1) = %.6f\n\n', p_b1);

%% ========================================================================
%  5. WITHIN-PARTICIPANT LOGISTIC REGRESSION
%  ========================================================================

uParts = unique(m34.PartID);
nP     = length(uParts);

partCrossover = nan(nP, 1);
partBeta1     = nan(nP, 1);

for pi = 1:nP
    mask = m34.PartID == uParts(pi);
    ct   = m34.ClickTime(mask);
    pas  = double(m34.PAS(mask) == 4);

    if length(ct) < 5, continue; end
    if all(pas == 0) || all(pas == 1), continue; end

    try
        mdl_p = fitglm(ct, pas, 'Distribution', 'binomial');
        b0p = mdl_p.Coefficients.Estimate(1);
        b1p = mdl_p.Coefficients.Estimate(2);
        partBeta1(pi) = b1p;

        if b1p < 0
            xover = -b0p / b1p;
            if xover > 0 && xover < T
                partCrossover(pi) = xover;
            end
        end
    catch
    end
end

validXover = ~isnan(partCrossover);
nValid     = sum(validXover);
xovers     = partCrossover(validXover);

cross_mean   = mean(xovers);
cross_median = median(xovers);
cross_sem    = std(xovers) / sqrt(nValid);
ci_lo        = cross_mean - 1.96 * cross_sem;
ci_hi        = cross_mean + 1.96 * cross_sem;

validBeta1 = ~isnan(partBeta1);
pWilcox    = signrank(partBeta1(validBeta1));

fprintf('=== WITHIN-PARTICIPANT LOGISTIC (N = %d valid crossovers) ===\n', nValid);
fprintf('  Mean:   %.1f s (SEM = %.1f)\n', cross_mean, cross_sem);
fprintf('  Median: %.1f s\n', cross_median);
fprintf('  95%% CI: [%.1f, %.1f] s\n', ci_lo, ci_hi);
fprintf('  Wilcoxon p(beta1 != 0) = %.4f\n\n', pWilcox);

%% ========================================================================
%  6. EARLY VS LATE (LOCKED TAU)
%  ========================================================================

uPartsAll = unique(merged.PartID);
nPall     = length(uPartsAll);

propEarly = nan(nPall, 1);
propLate  = nan(nPall, 1);

for pi = 1:nPall
    mask = merged.PartID == uPartsAll(pi);
    ct   = merged.ClickTime(mask);
    pas  = merged.PAS(mask);

    earlyM = ct <= tPeak;
    lateM  = ct > tPeak;

    if sum(earlyM) >= 2 && sum(lateM) >= 2
        propEarly(pi) = mean(pas(earlyM) == 4);
        propLate(pi)  = mean(pas(lateM) == 4);
    end
end

valid_el  = ~isnan(propEarly);
pE = propEarly(valid_el);
pL = propLate(valid_el);
nEL = sum(valid_el);
diffEL = pE - pL;
[~, p_el, ~, stats_el] = ttest(diffEL);
dz = mean(diffEL) / std(diffEL);

fprintf('=== EARLY vs LATE (tPeak = %.1f s, N = %d) ===\n', tPeak, nEL);
fprintf('  Early PAS=4: %.3f,  Late: %.3f\n', mean(pE), mean(pL));
fprintf('  t(%d) = %.2f, p = %.4f, d_z = %.3f\n\n', nEL-1, stats_el.tstat, p_el, dz);

%% ========================================================================
%  7. WITHIN-PARTICIPANT CLICK-TIME / PAS SPEARMAN CORRELATION
%  ========================================================================
%  Per participant (>= 4 clicks), the Spearman correlation between click
%  time and PAS rating (all ratings, not just 3/4); then a one-sample
%  t-test of the per-participant rho's against 0. Negative rho = subjective
%  clarity declines as the trial unfolds. Paper: rho = -0.12, t(57) = -2.77,
%  p = 0.008, N = 58.

uPartsSp = unique(merged.PartID);
partRho  = nan(length(uPartsSp), 1);

for pi = 1:length(uPartsSp)
    mask = merged.PartID == uPartsSp(pi);
    ct   = merged.ClickTime(mask);
    pas  = merged.PAS(mask);
    if length(ct) < 4, continue; end          % require >= 4 clicks
    if numel(unique(pas)) < 2 || numel(unique(ct)) < 2
        continue;                              % rho undefined for a constant input
    end
    partRho(pi) = corr(ct, pas, 'Type', 'Spearman');
end

validRho = ~isnan(partRho);
rhoVals  = partRho(validRho);
nRho     = numel(rhoVals);

rho_mean = mean(rhoVals);
[~, p_rho, ~, stats_rho] = ttest(rhoVals);     % one-sample t-test vs 0
t_rho    = stats_rho.tstat;

fprintf('=== WITHIN-PARTICIPANT CLICK-TIME / PAS SPEARMAN (N = %d) ===\n', nRho);
fprintf('  mean rho = %.4f\n', rho_mean);
fprintf('  t(%d) = %.4f, p = %.4f\n\n', nRho-1, t_rho, p_rho);

%% ========================================================================
%  8. SAVE RESULTS
%  ========================================================================

results = table( ...
    TAU_LOCKED, tPeak, height(m34), b0, b1, tcross_trial, p_b1, ...
    nValid, cross_mean, cross_median, cross_sem, ci_lo, ci_hi, pWilcox, ...
    nEL, mean(pE), mean(pL), stats_el.tstat, p_el, dz, ...
    nRho, rho_mean, t_rho, p_rho, ...
    'VariableNames', { ...
    'tau_locked', 'tPeak', 'n_trials_34', ...
    'trial_beta0', 'trial_beta1', 'trial_crossover_s', 'trial_beta1_p', ...
    'n_valid_crossovers', 'within_part_mean_s', 'within_part_median_s', ...
    'within_part_sem_s', 'within_part_ci_lo_s', 'within_part_ci_hi_s', ...
    'wilcoxon_beta1_p', ...
    'n_early_late', 'early_pas4_prop', 'late_pas4_prop', ...
    'early_late_t', 'early_late_p', 'early_late_dz', ...
    'n_spearman', 'spearman_rho_mean', 'spearman_t', 'spearman_p'});

outFile = fullfile(resDir, 'Figure3_crossover_stats.csv');
writetable(results, outFile);
fprintf('Saved: %s\n', outFile);
