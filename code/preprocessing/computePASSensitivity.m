%% gsp_sensitivity_pas_figure.m
% =========================================================================
% Sensitivity Analysis of P(0) and PAS Connection: Publication Figure
% =========================================================================
%
% PURPOSE:
%   Creates a 3-panel publication figure showing:
%     (A) Mathematical relationship: P(0) decay and its rate sensitivity
%         |dP(0)/dlambda| with peak at x = 1/lambda
%     (B) PAS clarity (prop PAS=4) for early vs late clicks relative to
%         the sensitivity peak (within-participant paired comparison)
%     (C) Moving-window PAS=4 proportion over trial time overlaid with
%         the sensitivity curve
%
% INPUT:
%   EEG/globalScalpPotential_stats.mat  (for tau, grand fit parameters)
%   Behavior/pceXXYYMMDD/questionnaires/pair_XX_P{1,2}_PAS_confidence_absence.csv
%   ../PCE optimal waiting analysis/Repository/data/ClickTimes/ClickResponseTimes.csv
%
% OUTPUT:
%   EEG/gsp_fig5_sensitivity_pas.png
%   EEG/gsp_sensitivity_pas.mat  (all statistical results)
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   April 2026
% =========================================================================

clearvars; close all; clc;

%% ========================================================================
%  1. LOAD GSP STATS (for tau)
%  ========================================================================

statsFile = fullfile(pwd, 'EEG', 'globalScalpPotential_stats.mat');
S = load(statsFile);

tau    = S.grandTau;       % onset lag from grand P(1) fit
T      = 60;               % trial duration (s)
lambda = exp(1);           % fixed rate parameter
tPeak  = tau + (T - tau) / lambda;  % sensitivity peak in real time

fprintf('Sensitivity peak: %.1f s (tau = %.1f s)\n', tPeak, tau);

%% ========================================================================
%  2. LOAD CLICK TIMES
%  ========================================================================

clickFile = fullfile(pwd, '..', 'PCE optimal waiting analysis', ...
    'Repository', 'data', 'ClickTimes', 'ClickResponseTimes.csv');
CT = readtable(clickFile);

% Extract dyad number from DyadID (format XXYYYYYY)
CT.DyadNum = floor(CT.DyadID / 1e6);

% Exclude dyad 31
CT = CT(CT.DyadNum ~= 31, :);

fprintf('Loaded %d click entries from %d unique dyads\n', ...
    height(CT), length(unique(CT.DyadNum)));

%% ========================================================================
%  3. LOAD PAS DATA (from Behavior questionnaires)
%  ========================================================================

behavDir = fullfile(pwd, 'Behavior');
allDyads = dir(behavDir);
allDyads = allDyads([allDyads.isdir]);

% Build PAS table: DyadNum, ParticipantID, TrialNum (1-indexed), PAS
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

        % Extract PAS values (question_id == 'click_presence')
        pasMask = tbl.question_id == "click_presence";
        if sum(pasMask) == 0, continue; end

        trials_0idx = tbl.trial_id(pasMask);
        pasVals     = tbl.answer(pasMask);

        for ti = 1:length(trials_0idx)
            pasRows = [pasRows; dN, p, trials_0idx(ti)+1, pasVals(ti)]; %#ok<AGROW>
        end
    end
end

PAS = array2table(pasRows, 'VariableNames', {'DyadNum','ParticipantID','TrialNum','PAS'});
fprintf('Loaded %d PAS entries\n', height(PAS));

%% ========================================================================
%  4. MERGE CLICK TIMES WITH PAS
%  ========================================================================

% Create matching keys
CT.Key  = CT.DyadNum * 1000 + CT.ParticipantID * 100 + CT.TrialNum;
PAS.Key = PAS.DyadNum * 1000 + PAS.ParticipantID * 100 + PAS.TrialNum;

% Inner join
[~, iCT, iPAS] = intersect(CT.Key, PAS.Key);
merged = table();
merged.DyadNum       = CT.DyadNum(iCT);
merged.ParticipantID = CT.ParticipantID(iCT);
merged.TrialNum      = CT.TrialNum(iCT);
merged.ClickTime     = CT.ClickTime_s(iCT);
merged.Clicked       = CT.Clicked(iCT);
merged.PAS           = PAS.PAS(iPAS);

% Keep only trials with clicks
merged = merged(merged.Clicked == 1 & ~isnan(merged.ClickTime), :);
fprintf('Merged dataset: %d trials with clicks and PAS\n', height(merged));

%% ========================================================================
%  5. WITHIN-PARTICIPANT ANALYSIS: EARLY vs LATE
%  ========================================================================

merged.PartID = merged.DyadNum * 10 + merged.ParticipantID;
uParts = unique(merged.PartID);
nP = length(uParts);

propPAS4_early = nan(nP, 1);
propPAS4_late  = nan(nP, 1);
meanPAS_early  = nan(nP, 1);
meanPAS_late   = nan(nP, 1);
nEarly         = nan(nP, 1);
nLate          = nan(nP, 1);

for pi = 1:nP
    mask = merged.PartID == uParts(pi);
    ct   = merged.ClickTime(mask);
    pas  = merged.PAS(mask);

    earlyMask = ct <= tPeak;
    lateMask  = ct > tPeak;

    if sum(earlyMask) >= 2 && sum(lateMask) >= 2
        propPAS4_early(pi) = mean(pas(earlyMask) == 4);
        propPAS4_late(pi)  = mean(pas(lateMask) == 4);
        meanPAS_early(pi)  = mean(pas(earlyMask));
        meanPAS_late(pi)   = mean(pas(lateMask));
        nEarly(pi)         = sum(earlyMask);
        nLate(pi)          = sum(lateMask);
    end
end

% Remove participants without both early and late
valid = ~isnan(propPAS4_early);
propPAS4_early = propPAS4_early(valid);
propPAS4_late  = propPAS4_late(valid);
meanPAS_early  = meanPAS_early(valid);
meanPAS_late   = meanPAS_late(valid);
nPAS = sum(valid);

fprintf('\n=== EARLY vs LATE (N = %d) ===\n', nPAS);

% Paired t-test: prop PAS=4
diff_prop = propPAS4_early - propPAS4_late;
[~, pval_prop, ~, stats_prop] = ttest(diff_prop);
dz_prop = mean(diff_prop) / std(diff_prop);
fprintf('Prop(PAS=4): early=%.3f, late=%.3f\n', ...
    mean(propPAS4_early), mean(propPAS4_late));
fprintf('  t(%d) = %.2f, p = %.4f, d_z = %.3f\n', ...
    nPAS-1, stats_prop.tstat, pval_prop, dz_prop);

% Paired t-test: mean PAS
diff_mean = meanPAS_early - meanPAS_late;
[~, pval_mean, ~, stats_mean] = ttest(diff_mean);
dz_mean = mean(diff_mean) / std(diff_mean);
fprintf('Mean PAS: early=%.3f, late=%.3f\n', ...
    mean(meanPAS_early), mean(meanPAS_late));
fprintf('  t(%d) = %.2f, p = %.4f, d_z = %.3f\n', ...
    nPAS-1, stats_mean.tstat, pval_mean, dz_mean);

% Within-participant Spearman: click time vs PAS
spearRho = nan(nP, 1);
for pi = 1:nP
    mask = merged.PartID == uParts(pi);
    ct = merged.ClickTime(mask);
    pas = merged.PAS(mask);
    if length(ct) >= 5
        spearRho(pi) = corr(ct, pas, 'Type', 'Spearman');
    end
end
validRho = ~isnan(spearRho);
[~, pval_rho, ~, stats_rho] = ttest(spearRho(validRho));
nRho = sum(validRho);
dz_rho = mean(spearRho(validRho)) / std(spearRho(validRho));
fprintf('Spearman (click time x PAS): mean rho = %.3f\n', ...
    mean(spearRho(validRho)));
fprintf('  t(%d) = %.2f, p = %.4f, d_z = %.3f\n', ...
    nRho-1, stats_rho.tstat, pval_rho, dz_rho);

%% ========================================================================
%  6. MOVING-WINDOW PAS OVER TRIAL TIME
%  ========================================================================

winSize    = 10;   % seconds
winStep    = 1;    % seconds
winCenters = (winSize/2):winStep:(T - winSize/2);
nWin       = length(winCenters);

winPropPAS4 = nan(nP, nWin);
for pi = 1:nP
    mask = merged.PartID == uParts(pi);
    ct   = merged.ClickTime(mask);
    pas  = merged.PAS(mask);

    for wi = 1:nWin
        wMask = ct >= (winCenters(wi) - winSize/2) & ...
                ct <  (winCenters(wi) + winSize/2);
        if sum(wMask) >= 2
            winPropPAS4(pi, wi) = mean(pas(wMask) == 4);
        end
    end
end

grandPropPAS4 = mean(winPropPAS4, 1, 'omitnan');
grandSEM      = std(winPropPAS4, 0, 1, 'omitnan') ./ ...
                sqrt(sum(~isnan(winPropPAS4), 1));

%% ========================================================================
%  7. SAVE RESULTS
%  ========================================================================

save(fullfile(pwd, 'EEG', 'gsp_sensitivity_pas.mat'), ...
    'propPAS4_early', 'propPAS4_late', 'meanPAS_early', 'meanPAS_late', ...
    'nPAS', 'tPeak', 'tau', ...
    'stats_prop', 'pval_prop', 'dz_prop', ...
    'stats_mean', 'pval_mean', 'dz_mean', ...
    'spearRho', 'stats_rho', 'pval_rho', 'dz_rho', 'nRho', ...
    'winCenters', 'winPropPAS4', 'grandPropPAS4', 'grandSEM', ...
    'merged', 'uParts');

fprintf('\nSaved: EEG/gsp_sensitivity_pas.mat\n');

%% ========================================================================
%  8. STYLING
%  ========================================================================

col_decay   = [0.50 0.50 0.50];   % Grey for P(0)
col_sens    = [0.85 0.20 0.10];   % Red for sensitivity curve
col_early   = [0.15 0.45 0.75];   % Blue for early clicks
col_late    = [0.85 0.55 0.10];   % Orange for late clicks
col_pas4    = [0.20 0.65 0.32];   % Green for PAS=4

font_main   = 11;
font_label  = 12;
font_title  = 13;

%% ========================================================================
%  9. FIGURE
%  ========================================================================

fig = figure('Units','inches', 'Position',[0.5 0.5 7.5 9.5], 'Color','w', ...
    'PaperUnits','inches', 'PaperSize',[7.5 9.5], 'PaperPosition',[0 0 7.5 9.5]);

% ---- Panel A: Mathematical Derivation ----
ax_a = subplot(3,1,1);
hold on;

x = linspace(0, 1.5, 300);
P0   = exp(-lambda * x);
sens = x .* exp(-lambda * x);
sens_norm = sens / max(sens);

yyaxis left;
plot(x, P0, '-', 'Color', col_decay, 'LineWidth', 2.5);
ylabel('P(0) = e^{-\lambda x}', 'FontSize', font_label, 'Color', col_decay);
set(gca, 'YColor', col_decay);
ylim([-0.05 1.05]);

yyaxis right;
plot(x, sens_norm, '-', 'Color', col_sens, 'LineWidth', 2.5);
ylabel('|\partial P(0) / \partial\lambda|  (normalized)', ...
    'FontSize', font_label, 'Color', col_sens);
set(gca, 'YColor', col_sens);
ylim([-0.05 1.25]);

% Mark peak
xPeak = 1/lambda;
plot(xPeak, 1.0, 'v', 'MarkerSize', 10, 'MarkerFaceColor', col_sens, ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
text(xPeak + 0.04, 1.08, sprintf('x = 1/\\lambda = 1/e', xPeak), ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', col_sens);

% Shade early/late regions
yl = ylim;
fill([0 xPeak xPeak 0], [yl(1) yl(1) yl(2) yl(2)], ...
    col_early, 'FaceAlpha', 0.06, 'EdgeColor', 'none');
fill([xPeak x(end) x(end) xPeak], [yl(1) yl(1) yl(2) yl(2)], ...
    col_late, 'FaceAlpha', 0.06, 'EdgeColor', 'none');

% Region labels
text(0.15, 0.85, {'System robust', '(low sensitivity)'}, ...
    'FontSize', 9, 'Color', col_early, 'FontAngle', 'italic', ...
    'Units', 'normalized');
text(0.72, 0.85, {'System sensitive', '(high sensitivity)'}, ...
    'FontSize', 9, 'Color', col_late, 'FontAngle', 'italic', ...
    'Units', 'normalized');

hold off;
xlabel('Normalized time  x = (t - \tau) / (T - \tau)', 'FontSize', font_label);
title('(A)  Sensitivity of P(0) to rate perturbations', ...
    'FontSize', font_title, 'FontWeight', 'bold');
set(gca, 'TickDir', 'out', 'FontSize', font_main, 'Box', 'off');
xlim([0 1.2]);
lg = legend('P(0):  probability of no event', ...
            '|\partialP(0)/\partial\lambda|:  rate sensitivity', ...
            'Location', 'east');
lg.FontSize = 9; lg.Box = 'off';

% ---- Panel B: PAS Early vs Late ----
ax_b = subplot(3,1,2);
hold on;

barData = [mean(propPAS4_early), mean(propPAS4_late)];
barSEM  = [std(propPAS4_early)/sqrt(nPAS), std(propPAS4_late)/sqrt(nPAS)];

b = bar([1 2], barData, 0.55);
b.FaceColor = 'flat';
b.CData(1,:) = col_early;
b.CData(2,:) = col_late;
b.EdgeColor  = 'k';
b.LineWidth  = 0.8;

errorbar([1 2], barData, barSEM, 'k', 'LineStyle', 'none', ...
    'LineWidth', 1.5, 'CapSize', 8);

% Individual participants (jittered)
rng(42);
jE = 0.15 * (rand(nPAS,1) - 0.5);
jL = 0.15 * (rand(nPAS,1) - 0.5);
scatter(1 + jE, propPAS4_early, 12, col_early, 'filled', 'MarkerFaceAlpha', 0.3);
scatter(2 + jL, propPAS4_late,  12, col_late,  'filled', 'MarkerFaceAlpha', 0.3);

% Paired lines
for i = 1:nPAS
    plot([1+jE(i), 2+jL(i)], [propPAS4_early(i), propPAS4_late(i)], ...
        '-', 'Color', [0.5 0.5 0.5 0.12], 'LineWidth', 0.4);
end

% Significance bracket
yBracket = max(barData) + max(barSEM) + 0.04;
plot([1 1 2 2], [yBracket-0.005 yBracket yBracket yBracket-0.005], ...
    'k-', 'LineWidth', 1);
if pval_prop < 0.001
    sigStr = '***';
elseif pval_prop < 0.01
    sigStr = '**';
elseif pval_prop < 0.05
    sigStr = '*';
else
    sigStr = 'n.s.';
end
text(1.5, yBracket + 0.005, sigStr, 'HorizontalAlignment', 'center', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Stats text
text(1.5, yBracket + 0.07, ...
    sprintf('t(%d) = %.2f, p = %.4f, d_z = %.2f', ...
    nPAS-1, stats_prop.tstat, pval_prop, dz_prop), ...
    'HorizontalAlignment', 'center', 'FontSize', 9);

hold off;
set(gca, 'XTick', [1 2], 'XTickLabel', ...
    {sprintf('Early (before %.0f s)', tPeak), ...
     sprintf('Late (after %.0f s)', tPeak)});
ylabel('Proportion PAS = 4 (clear experience)', 'FontSize', font_label);
title(sprintf('(B)  Perceptual clarity: early vs. late clicks  (N = %d)', nPAS), ...
    'FontSize', font_title, 'FontWeight', 'bold');
ylim([0 yBracket + 0.14]);
set(gca, 'TickDir', 'out', 'FontSize', font_main, 'Box', 'off');

% ---- Panel C: Moving-window PAS + sensitivity curve ----
ax_c = subplot(3,1,3);
hold on;

% Plot PAS=4 proportion over time (left axis)
yyaxis left;
validW = ~isnan(grandPropPAS4) & ~isnan(grandSEM);
xw = winCenters(validW);
yw = grandPropPAS4(validW);
ew = grandSEM(validW);
fill([xw, fliplr(xw)], [yw+ew, fliplr(yw-ew)], ...
    col_pas4, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot(xw, yw, '-', 'Color', col_pas4, 'LineWidth', 2);
ylabel('Proportion PAS = 4', 'FontSize', font_label, 'Color', [0 0.4 0.15]);
set(gca, 'YColor', [0 0.4 0.15]);

% Sensitivity curve on right axis
yyaxis right;
tSm = S.tTaskSm;
tNorm = (tSm - tau) / (T - tau);
sensReal = tNorm .* exp(-lambda * tNorm);
sensReal(tNorm < 0) = 0;
sensReal = sensReal / max(sensReal);
plot(tSm, sensReal, '--', 'Color', col_sens, 'LineWidth', 2.5);
ylabel('Rate sensitivity (normalized)', 'FontSize', font_label, 'Color', col_sens);
set(gca, 'YColor', col_sens);
ylim([-0.05 1.3]);

% Mark peak
plot(tPeak, 1.0, 'v', 'MarkerSize', 9, 'MarkerFaceColor', col_sens, ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
xline(tPeak, ':', 'Color', [col_sens 0.5], 'LineWidth', 1);
text(tPeak + 1, 1.05, sprintf('%.1f s', tPeak), ...
    'FontSize', 9, 'FontWeight', 'bold', 'Color', col_sens);

hold off;
xlabel('Trial time (s)', 'FontSize', font_label);
title('(C)  Perceptual clarity tracks sensitivity dynamics', ...
    'FontSize', font_title, 'FontWeight', 'bold');
xlim([0 T]);
set(gca, 'TickDir', 'out', 'FontSize', font_main, 'Box', 'off');
lg2 = legend('Prop(PAS=4) \pm SEM', '', 'Rate sensitivity', ...
    'Location', 'northeast');
lg2.FontSize = 9; lg2.Box = 'off';

%% ========================================================================
%  10. SAVE FIGURE
%  ========================================================================

outFile = fullfile(pwd, 'EEG', 'gsp_fig5_sensitivity_pas.png');
exportgraphics(fig, outFile, 'Resolution', 300);
fprintf('\nSaved figure: %s\n', outFile);
close(fig);

fprintf('Done.\n');
