%% gsp_pas_roi_crossover.m
% =========================================================================
% PAS 4 vs PAS 3 Crossover: Within-Participant Analysis + ROI Mapping
% =========================================================================
%
% PURPOSE:
%   1. Replicate the PAS 4 vs PAS 3 logistic crossover analysis using
%      within-participant methods (avoiding pseudoreplication)
%   2. Compare the crossover time with per-ROI P(1) peak times from the
%      global scalp potential analysis
%   3. Create a publication figure linking neural and perceptual timescales
%
% METHODS:
%   - Per-participant logistic regression: fit PAS4 vs PAS3 binomial GLM
%     to each participant separately, extract crossover time
%   - Grand crossover: median of per-participant crossovers
%   - Also: within-participant moving-window PAS 4 proportion (hierarchical)
%   - ROI comparison: overlay ROI peak times on the PAS time course
%
% INPUT:
%   EEG/globalScalpPotential_stats.mat
%   ../PCE optimal waiting analysis/Repository/data/ClickTimes/ClickResponseTimes.csv
%   Behavior/pce*/questionnaires/pair_*_PAS_confidence_absence.csv
%
% OUTPUT:
%   EEG/gsp_fig6_pas_roi_crossover.png
%   EEG/gsp_pas_roi_crossover.mat
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   April 2026
% =========================================================================

clearvars; close all; clc;

%% ========================================================================
%  1. LOAD GSP STATS (ROI peak times)
%  ========================================================================

S = load(fullfile(pwd, 'EEG', 'globalScalpPotential_stats.mat'));
roiNames = S.roi.names;
roiPeaks = S.roiGrandPeak;
roiR2    = S.roiGrandR2;
nROI     = S.nROI;
tTaskSm  = S.tTaskSm;

globalPeak = S.grandPeak;
globalTau  = S.grandTau;
T = 60;
lambda = exp(1);

fprintf('Global GSP peak: %.1f s (tau = %.1f s)\n', globalPeak, globalTau);
for r = 1:nROI
    fprintf('  %-18s peak = %5.1f s  (R2 = %.3f)\n', roiNames{r}, roiPeaks(r), roiR2(r));
end

%% ========================================================================
%  2. LOAD AND MERGE CLICK TIMES + PAS
%  ========================================================================

clickFile = fullfile(pwd, '..', 'PCE optimal waiting analysis', ...
    'Repository', 'data', 'ClickTimes', 'ClickResponseTimes.csv');
CT = readtable(clickFile);
CT.DyadNum = floor(CT.DyadID / 1e6);
CT = CT(CT.DyadNum ~= 31, :);

behavDir = fullfile(pwd, 'Behavior');
allDyads = dir(behavDir);
allDyads = allDyads([allDyads.isdir]);

pasRows = [];
for i = 1:length(allDyads)
    tok = regexp(allDyads(i).name, '^pce(\d{2})\d{6}$', 'tokens');
    if isempty(tok), continue; end
    dN = str2double(tok{1}{1});
    if dN == 31, continue; end
    for p = [1 2]
        fn = sprintf('pair_%02d_P%d_PAS_confidence_absence.csv', dN, p);
        fp = fullfile(behavDir, allDyads(i).name, 'questionnaires', fn);
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
PAS = array2table(pasRows, 'VariableNames', {'DyadNum','ParticipantID','TrialNum','PAS'});

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

fprintf('\nMerged: %d click-PAS trials\n', height(merged));

%% ========================================================================
%  3. TRIAL-LEVEL LOGISTIC REGRESSION (for reference, as in paper)
%  ========================================================================

% PAS 3 and 4 only
m34 = merged(merged.PAS == 3 | merged.PAS == 4, :);
binaryPAS = double(m34.PAS == 4);

mdl_trial = fitglm(m34.ClickTime, binaryPAS, 'Distribution', 'binomial');
b0_trial = mdl_trial.Coefficients.Estimate(1);
b1_trial = mdl_trial.Coefficients.Estimate(2);
tcross_trial = -b0_trial / b1_trial;

fprintf('\n=== TRIAL-LEVEL LOGISTIC (reference) ===\n');
fprintf('  n = %d (PAS 3 & 4)\n', height(m34));
fprintf('  beta0 = %.4f, beta1 = %.4f\n', b0_trial, b1_trial);
fprintf('  Crossover = %.1f s\n', tcross_trial);

%% ========================================================================
%  4. WITHIN-PARTICIPANT LOGISTIC REGRESSION
%  ========================================================================

uParts = unique(m34.PartID);
nP = length(uParts);

partCrossover = nan(nP, 1);
partBeta0     = nan(nP, 1);
partBeta1     = nan(nP, 1);
partN34       = nan(nP, 1);

for pi = 1:nP
    mask = m34.PartID == uParts(pi);
    ct  = m34.ClickTime(mask);
    pas = double(m34.PAS(mask) == 4);

    % Need at least 5 trials with both PAS 3 and PAS 4
    if length(ct) < 5, continue; end
    if all(pas == 0) || all(pas == 1), continue; end

    try
        mdl_p = fitglm(ct, pas, 'Distribution', 'binomial');
        b0p = mdl_p.Coefficients.Estimate(1);
        b1p = mdl_p.Coefficients.Estimate(2);

        % Crossover only meaningful if slope is negative
        if b1p < 0
            xover = -b0p / b1p;
            if xover > 0 && xover < T
                partCrossover(pi) = xover;
            end
        end
        partBeta0(pi) = b0p;
        partBeta1(pi) = b1p;
        partN34(pi)   = sum(mask);
    catch
        % Convergence failure — skip
    end
end

validXover = ~isnan(partCrossover);
nValid = sum(validXover);

grandCrossover_mean   = mean(partCrossover(validXover));
grandCrossover_median = median(partCrossover(validXover));
grandCrossover_SEM    = std(partCrossover(validXover)) / sqrt(nValid);

% Test whether slope is consistently negative
validBeta1 = ~isnan(partBeta1);
[~, pBeta1, ~, statsBeta1] = ttest(partBeta1(validBeta1));

fprintf('\n=== WITHIN-PARTICIPANT LOGISTIC (N = %d with valid crossover) ===\n', nValid);
fprintf('  Mean crossover:   %.1f s (SEM = %.1f)\n', grandCrossover_mean, grandCrossover_SEM);
fprintf('  Median crossover: %.1f s\n', grandCrossover_median);
fprintf('  Range: [%.1f, %.1f] s\n', min(partCrossover(validXover)), max(partCrossover(validXover)));
fprintf('  Beta1 across participants: mean = %.4f\n', mean(partBeta1(validBeta1)));
fprintf('  t(%d) = %.2f, p = %.4f\n', sum(validBeta1)-1, statsBeta1.tstat, pBeta1);

% Wilcoxon signed-rank as robustness check
pWilcox_beta1 = signrank(partBeta1(validBeta1));
fprintf('  Wilcoxon p(beta1 != 0) = %.4f\n', pWilcox_beta1);

%% ========================================================================
%  5. HIERARCHICAL MOVING-WINDOW PAS PROPORTIONS
%  ========================================================================

winSize    = 8;
winStep    = 1;
winCenters = (winSize/2):winStep:(T - winSize/2);
nWin       = length(winCenters);

uPartsAll = unique(merged.PartID);
nPall     = length(uPartsAll);

% Per-participant, per-window proportions for each PAS level
winProp = cell(4, 1);
for lv = 1:4
    winProp{lv} = nan(nPall, nWin);
end

for pi = 1:nPall
    mask = merged.PartID == uPartsAll(pi);
    ct   = merged.ClickTime(mask);
    pas  = merged.PAS(mask);

    for wi = 1:nWin
        wMask = ct >= (winCenters(wi) - winSize/2) & ...
                ct <  (winCenters(wi) + winSize/2);
        nw = sum(wMask);
        if nw >= 2
            for lv = 1:4
                winProp{lv}(pi, wi) = mean(pas(wMask) == lv);
            end
        end
    end
end

% Grand means across participants
grandProp = nan(4, nWin);
grandSEM  = nan(4, nWin);
for lv = 1:4
    grandProp(lv, :) = mean(winProp{lv}, 1, 'omitnan');
    nValid_w = sum(~isnan(winProp{lv}), 1);
    grandSEM(lv, :)  = std(winProp{lv}, 0, 1, 'omitnan') ./ sqrt(nValid_w);
end

% Convert to percentages for Panel A
grandPropPct = grandProp * 100;
grandSEMPct  = grandSEM * 100;

%% ========================================================================
%  6. FIND HIERARCHICAL PAS4-PAS3 CROSSOVER FROM MOVING WINDOW
%  ========================================================================

diff43 = grandProp(4, :) - grandProp(3, :);
% Find where diff43 crosses zero (PAS4 drops below PAS3)
crossIdx = find(diff43(1:end-1) > 0 & diff43(2:end) <= 0, 1);
if ~isempty(crossIdx)
    % Linear interpolation
    x1 = winCenters(crossIdx); x2 = winCenters(crossIdx+1);
    y1 = diff43(crossIdx); y2 = diff43(crossIdx+1);
    crossover_mw = x1 + (0 - y1) * (x2 - x1) / (y2 - y1);
    fprintf('\nMoving-window PAS4-PAS3 crossover: %.1f s\n', crossover_mw);
else
    crossover_mw = NaN;
    fprintf('\nNo PAS4-PAS3 crossover found in moving window\n');
end

%% ========================================================================
%  7. ROI-TO-PAS COMPARISON
%  ========================================================================

fprintf('\n=== ROI PEAKS vs PAS CROSSOVER ===\n');
fprintf('PAS4-PAS3 crossover (trial-level logistic): %.1f s\n', tcross_trial);
fprintf('PAS4-PAS3 crossover (within-part median):   %.1f s\n', grandCrossover_median);
fprintf('PAS4-PAS3 crossover (moving window):        %.1f s\n', crossover_mw);
fprintf('\n');
fprintf('%-18s  peak(s)  |diff from crossover|\n', 'ROI');
fprintf('%-18s  ------   ----\n', '---');
for r = 1:nROI
    fprintf('%-18s  %5.1f    %.1f s\n', roiNames{r}, roiPeaks(r), ...
        abs(roiPeaks(r) - tcross_trial));
end
fprintf('\nClosest ROI to PAS crossover: ');
[~, bestROI] = min(abs(roiPeaks - tcross_trial));
fprintf('%s (peak = %.1f s, diff = %.1f s)\n', roiNames{bestROI}, ...
    roiPeaks(bestROI), abs(roiPeaks(bestROI) - tcross_trial));

%% ========================================================================
%  8. SAVE RESULTS
%  ========================================================================

save(fullfile(pwd, 'EEG', 'gsp_pas_roi_crossover.mat'), ...
    'tcross_trial', 'b0_trial', 'b1_trial', ...
    'partCrossover', 'partBeta0', 'partBeta1', 'partN34', ...
    'grandCrossover_mean', 'grandCrossover_median', 'grandCrossover_SEM', ...
    'nValid', 'statsBeta1', 'pBeta1', 'pWilcox_beta1', ...
    'winCenters', 'winProp', 'grandProp', 'grandSEM', ...
    'crossover_mw', 'roiPeaks', 'roiNames', 'roiR2', ...
    'globalPeak', 'merged');

fprintf('\nSaved: EEG/gsp_pas_roi_crossover.mat\n');

%% ========================================================================
%  9. FIGURE
%  ========================================================================

% --- Styling ---
col_pas = [0.65 0.65 0.65;   % PAS 1 grey
           0.93 0.69 0.13;   % PAS 2 amber
           0.85 0.33 0.10;   % PAS 3 orange
           0.00 0.00 0.00];  % PAS 4 black
lw_pas  = [1.5 1.5 2.0 2.5];

col_sens   = [0.85 0.20 0.10];
col_global = [0.15 0.45 0.75];

roi_colors = [
    0.90 0.29 0.21;   % Prefrontal - red
    0.95 0.50 0.25;   % Frontal - orange
    0.20 0.65 0.32;   % Fronto-Central - green
    0.00 0.53 0.53;   % Central - teal
    0.24 0.33 0.63;   % Centro-Parietal - blue
    0.52 0.40 0.70;   % Parietal - purple
    0.49 0.38 0.28;   % Occipital - brown
];

font_main  = 11;
font_label = 12;
font_title = 13;
font_panel = 18;

fig = figure('Units','inches', 'Position',[0.5 0.5 8 11], 'Color','w', ...
    'PaperUnits','inches', 'PaperSize',[8 11], 'PaperPosition',[0 0 8 11]);

% ---- Panel A: PAS proportions over time (hierarchical) ----
ax_a = subplot(3,1,1);
hold on;
for lv = 1:4
    validW = ~isnan(grandPropPct(lv,:));
    xw = winCenters(validW);
    yw = grandPropPct(lv, validW);
    ew = grandSEMPct(lv, validW);
    fill([xw, fliplr(xw)], [yw+ew, fliplr(yw-ew)], ...
        col_pas(lv,:), 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
    plot(xw, yw, '-', 'Color', col_pas(lv,:), 'LineWidth', lw_pas(lv), ...
        'DisplayName', sprintf('PAS %d', lv));
end

% Mark within-participant logistic crossover (median)
xline(grandCrossover_median, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');
text(grandCrossover_median + 1, 53, ...
    sprintf('PAS crossover\n%.0f s', grandCrossover_median), ...
    'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.3 0.3 0.3]);

% Mark fronto-central ROI peak
xline(roiPeaks(3), ':', 'Color', roi_colors(3,:), 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');
text(roiPeaks(3) - 1, 53, ...
    sprintf('FC peak\n%.1f s', roiPeaks(3)), ...
    'FontSize', 9, 'FontWeight', 'bold', 'Color', roi_colors(3,:), ...
    'HorizontalAlignment', 'right');

hold off;
xlim([0 T]); ylim([0 58]);
ylabel('Proportion of clicks (%)', 'FontSize', font_label);
title(sprintf('(A)  PAS ratings over trial time  (%d s window, within-participant)', winSize), ...
    'FontSize', font_title, 'FontWeight', 'bold');
lgd = legend('Location', 'east', 'Box', 'on', 'FontSize', 10);
set(lgd, 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);
title(lgd, 'PAS Rating');
set(gca, 'TickDir', 'out', 'FontSize', font_main, 'Box', 'off');

% ---- Panel B: Per-participant crossover distribution + ROI peaks ----
ax_b = subplot(3,1,2);
hold on;

% Histogram of per-participant crossovers
validCross = partCrossover(validXover);
histogram(validCross, 15, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'w', ...
    'FaceAlpha', 0.6, 'HandleVisibility', 'off');

% Grand median line
yl = [0 max(histcounts(validCross, 15)) * 1.3];
ylim(yl);

plot([grandCrossover_median grandCrossover_median], yl, '-', ...
    'Color', 'k', 'LineWidth', 2, ...
    'DisplayName', sprintf('Median = %.1f s', grandCrossover_median));

% ROI peak markers at top — stagger vertically to avoid overlap
yLevels = yl(2) * [0.95, 0.87, 0.95, 0.87, 0.79, 0.87, 0.95];
% Short ROI labels
roiShort = {'PF','F','FC','C','CP','P','O'};
for r = 1:nROI
    plot(roiPeaks(r), yLevels(r), 'v', 'MarkerSize', 9, ...
        'MarkerFaceColor', roi_colors(r,:), 'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.6, 'HandleVisibility', 'off');
    text(roiPeaks(r), yLevels(r) + yl(2)*0.05, ...
        sprintf('%s\n%.1f', roiShort{r}, roiPeaks(r)), ...
        'FontSize', 7.5, 'HorizontalAlignment', 'center', ...
        'Color', roi_colors(r,:), 'FontWeight', 'bold');
end

% Trial-level crossover for reference
plot([tcross_trial tcross_trial], yl, '--', 'Color', [0.5 0.5 0.5], ...
    'LineWidth', 1.2, ...
    'DisplayName', sprintf('Trial-level = %.1f s', tcross_trial));

hold off;
xlim([0 T]);
xlabel('PAS 4 vs PAS 3 crossover time (s)', 'FontSize', font_label);
ylabel('Count', 'FontSize', font_label);
title(sprintf('(B)  Per-participant crossover times  (N = %d)  with ROI peak times', nValid), ...
    'FontSize', font_title, 'FontWeight', 'bold');
legend('Location', 'northwest', 'Box', 'off', 'FontSize', 10);
set(gca, 'TickDir', 'out', 'FontSize', font_main, 'Box', 'off');

% ---- Panel C: ROI peak times vs anterior-posterior position ----
ax_c = subplot(3,1,3);
hold on;

% ROI positions: 1=Prefrontal (most anterior) ... 7=Occipital (most posterior)
roiPos = 1:nROI;

% Bar chart of ROI peak times
for r = 1:nROI
    bar(r, roiPeaks(r), 0.6, 'FaceColor', roi_colors(r,:), ...
        'EdgeColor', 'k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
end

% Global GSP peak line
plot([0.3 nROI+0.7], [globalPeak globalPeak], '-', ...
    'Color', col_global, 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Global GSP peak = %.1f s', globalPeak));

% PAS crossover lines
plot([0.3 nROI+0.7], [tcross_trial tcross_trial], '--', ...
    'Color', 'k', 'LineWidth', 2, ...
    'DisplayName', sprintf('PAS 4/3 crossover = %.1f s', tcross_trial));

% CI band for PAS crossover
fill([0.3 nROI+0.7 nROI+0.7 0.3], [21.5 21.5 33.4 33.4], ...
    'k', 'FaceAlpha', 0.06, 'EdgeColor', 'none', ...
    'DisplayName', '95% CI [21.5, 33.4] s');

% R2 annotation on each bar
for r = 1:nROI
    text(r, roiPeaks(r) + 0.5, sprintf('.%03.0f', roiR2(r)*1000), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
end

hold off;
set(gca, 'XTick', 1:nROI, 'XTickLabel', roiNames, ...
    'XTickLabelRotation', 30);
xlim([0.3 nROI+0.7]);
ylim([20 32]);
xlabel('ROI (anterior \rightarrow posterior)', 'FontSize', font_label);
ylabel('P(1) peak time (s)', 'FontSize', font_label);
title('(C)  P(1) peak times by brain region with PAS crossover', ...
    'FontSize', font_title, 'FontWeight', 'bold');
legend('Location', 'southwest', 'Box', 'off', 'FontSize', 10);
set(gca, 'TickDir', 'out', 'FontSize', font_main, 'Box', 'off');

% No annotation arrows needed — axis label already says anterior -> posterior

%% ========================================================================
%  10. SAVE FIGURE
%  ========================================================================

outFile = fullfile(pwd, 'EEG', 'gsp_fig6_pas_roi_crossover.png');
exportgraphics(fig, outFile, 'Resolution', 300);
fprintf('\nSaved figure: %s\n', outFile);
close(fig);

fprintf('Done.\n');
