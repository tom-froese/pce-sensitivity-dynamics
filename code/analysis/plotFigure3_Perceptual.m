%% plotFigure3_Perceptual.m
% =========================================================================
% PNAS Figure 3: Perceptual Evidence — Sensitivity Dynamics and PAS
% =========================================================================
%
% Two-panel figure linking sensitivity dynamics to perceptual awareness:
%
%   Panel A — Sensitivity of R(x) to Rate Perturbations (Schematic)
%     The reliability function R(x) = exp(-lambda*x) and its sensitivity
%     dR/dlambda (normalized), which reaches its trough at x = 1/lambda.
%     At lambda = e, this trough occurs at x = 1/e — the optimal stopping
%     point. The left half is the "rejection phase" (low sensitivity,
%     system robust); the right half is the "selection phase" (high
%     sensitivity, system fragile).
%
%   Panel B — PAS Ratings Over Trial Time + Bootstrap Crossover CI
%     Moving-window PAS proportions (all 4 levels) across trial time,
%     showing how PAS 4 ("clear experience") declines from its initial
%     dominance. The PAS 4->3 crossover time is estimated via bootstrap
%     resampling of participants (median and 95% CI). ROI sensitivity
%     trough times are overlaid as markers.
%
% INPUT:
%   ../../data/preprocessed/EEG/globalScalpPotential_stats.mat
%   ../../data/preprocessed/PAS/gsp_pas_roi_crossover.mat
%
% OUTPUT:
%   ../../results/Figure3_Perceptual.png  (600 dpi)
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   April 2026
% =========================================================================

clear; clc; close all;

%% ========================================================================
%  1. LOAD PRECOMPUTED RESULTS
%  ========================================================================

statsFile     = '../../data/preprocessed/EEG/globalScalpPotential_stats.mat';
crossoverFile = '../../data/preprocessed/PAS/gsp_pas_roi_crossover.mat';

for f = {statsFile, crossoverFile}
    if ~isfile(f{1})
        error('Missing: %s\nSee README for instructions.', f{1});
    end
end

S    = load(statsFile);
XOVR = load(crossoverFile);

tau    = S.grandTau;
T      = 60;
lambda = exp(1);
tPeak  = tau + (T - tau) / lambda;

fprintf('Loaded precomputed results:\n');
fprintf('  GSP: N=%d, grand R^2=%.3f\n', S.nPart, S.grandR2);

%% ========================================================================
%  2. STYLING
%  ========================================================================

col_decay = [0.50 0.50 0.50];    % Grey for R(x)
col_sens  = [0.85 0.20 0.10];    % Red for sensitivity curve
col_early = [0.15 0.45 0.75];    % Blue
col_late  = [0.85 0.55 0.10];    % Orange

col_pas = [0.65 0.65 0.65;       % PAS 1 grey
           0.93 0.69 0.13;       % PAS 2 amber
           0.85 0.33 0.10;       % PAS 3 orange
           0.00 0.00 0.00];      % PAS 4 black
lw_pas  = [1.5 1.5 2.0 2.5];

% ROI gradient (same as Figure 2)
nROI = S.nROI;
baseColor  = [0.10 0.35 0.65];
lightColor = [0.70 0.85 0.95];
roi_colors = zeros(nROI, 3);
for r = 1:nROI
    frac = (r - 1) / (nROI - 1);
    roi_colors(r,:) = baseColor * (1 - frac) + lightColor * frac;
end
roiShort = {'PF','F','FC','C','CP','P','O'};

font_main  = 11;
font_label = 12;
font_title = 13;

%% ========================================================================
%  3. BOOTSTRAP CI FOR PAS 4/3 CROSSOVER
%  ========================================================================

winCenters = XOVR.winCenters;
nWin       = length(winCenters);
winProp3   = XOVR.winProp{3};    % nPart x nWin
winProp4   = XOVR.winProp{4};    % nPart x nWin
nPall      = size(winProp3, 1);

tcross_obs = XOVR.tcross_trial;
smKernel   = ones(1,3)/3;

nBoot = 10000;
rng(42);
bootCross = nan(nBoot, 1);
minCrossTime = 15;
minIdx = find(winCenters >= minCrossTime, 1);

for bi = 1:nBoot
    idx = randi(nPall, nPall, 1);
    gP3 = mean(winProp3(idx,:), 1, 'omitnan');
    gP4 = mean(winProp4(idx,:), 1, 'omitnan');
    diff43 = gP4 - gP3;
    diff43sm = conv(diff43, smKernel, 'same');
    for ci = minIdx:nWin-1
        if diff43sm(ci) > 0 && diff43sm(ci+1) <= 0
            x1 = winCenters(ci); x2 = winCenters(ci+1);
            y1 = diff43sm(ci); y2 = diff43sm(ci+1);
            bootCross(bi) = x1 + (0 - y1) * (x2 - x1) / (y2 - y1);
            break;
        end
    end
end

validBoot  = ~isnan(bootCross);
bootCI     = prctile(bootCross(validBoot), [2.5 97.5]);
bootMedian = median(bootCross(validBoot));

fprintf('\nBootstrap crossover (N=%d valid / %d):\n', sum(validBoot), nBoot);
fprintf('  Median = %.1f s,  95%% CI = [%.1f, %.1f] s\n', ...
    bootMedian, bootCI(1), bootCI(2));

%% ========================================================================
%  4. CREATE FIGURE — 2 rows
%  ========================================================================

fprintf('Creating Figure 3 ...\n');

fig = figure('Units', 'inches', 'Position', [0.5 0.5 7.5 6], 'Color', 'w', ...
    'PaperUnits', 'inches', 'PaperSize', [7.5 6], 'PaperPosition', [0 0 7.5 6]);

%% ---- Panel A: Reliability R(x) and its sensitivity ---------------------

ax_a = subplot(2,1,1);
hold on;

x = linspace(0, 1.0, 300);
Rx         = exp(-lambda * x);             % Reliability function
dRdlam     = -x .* exp(-lambda * x);       % dR/dlambda (negative)
dRdlam_norm = dRdlam / abs(min(dRdlam));   % Normalize so trough = -1
xTrough_norm = 1/lambda;

% Shade early/late FIRST so curves draw on top
yyaxis right;
ylim([-1.25 0.05]);
yl = ylim;
fill([0 xTrough_norm xTrough_norm 0], [yl(1) yl(1) yl(2) yl(2)], ...
    col_early, 'FaceAlpha', 0.06, 'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
fill([xTrough_norm x(end) x(end) xTrough_norm], [yl(1) yl(1) yl(2) yl(2)], ...
    col_late, 'FaceAlpha', 0.06, 'EdgeColor', 'none', ...
    'HandleVisibility', 'off');

text(0.08, 0.85, 'Rejection phase', ...
    'FontSize', 9, 'Color', col_early, 'FontAngle', 'italic', ...
    'Units', 'normalized');
text(0.58, 0.85, 'Selection phase', ...
    'FontSize', 9, 'Color', col_late, 'FontAngle', 'italic', ...
    'Units', 'normalized');

% Plot curves on top of shading
yyaxis left;
plot(x, Rx, '-', 'Color', col_decay, 'LineWidth', 2.5);
ylabel('R(x) = e^{-\lambda x}', 'FontSize', font_label, 'Color', col_decay);
set(gca, 'YColor', col_decay);
ylim([-0.05 1.05]);

yyaxis right;
plot(x, dRdlam_norm, '-', 'Color', col_sens, 'LineWidth', 3);
ylabel({'\partial R / \partial\lambda', '(normalized)'}, ...
    'FontSize', font_label, 'Color', col_sens);
set(gca, 'YColor', col_sens);
ylim([-1.25 0.05]);
set(gca, 'Color', 'none');
ax_a.Layer = 'top';

% Trough marker
plot(xTrough_norm, -1.0, '^', 'MarkerSize', 10, 'MarkerFaceColor', col_sens, ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
text(xTrough_norm + 0.03, -1.08, 'x = 1/\lambda = 1/e', ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', col_sens);

% Re-draw sensitivity curve on top
plot(x, dRdlam_norm, '-', 'Color', col_sens, 'LineWidth', 2.5, ...
    'HandleVisibility', 'off');

hold off;
xlabel('Normalized time  x = (t - \tau) / (T - \tau)', 'FontSize', font_label);
title('(A)  Sensitivity of reliability R(x) to rate perturbations', ...
    'FontSize', font_title, 'FontWeight', 'bold');
set(gca, 'TickDir', 'out', 'FontSize', font_main, 'Box', 'off');
xlim([0 1.0]);

lg = legend('R(x): reliability of metastable state', ...
            '\partialR/\partial\lambda: rate sensitivity', ...
            'Location', 'southeast');
lg.FontSize = 9; lg.Box = 'on';
lg.Color = 'w';
lg.EdgeColor = [0.7 0.7 0.7];

%% ---- Panel B: Moving-window PAS + bootstrap CI + ROI trough markers ----

ax_b = subplot(2,1,2);
hold on;

% PAS levels
pasHandles = gobjects(4,1);
for lv = 1:4
    validW = ~isnan(XOVR.grandProp(lv,:));
    xw = winCenters(validW);
    yw = XOVR.grandProp(lv, validW) * 100;
    ew = XOVR.grandSEM(lv, validW) * 100;
    fill([xw, fliplr(xw)], [yw+ew, fliplr(yw-ew)], ...
        col_pas(lv,:), 'FaceAlpha', 0.10, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
    pasHandles(lv) = plot(xw, yw, '-', 'Color', col_pas(lv,:), ...
        'LineWidth', lw_pas(lv), ...
        'DisplayName', sprintf('PAS %d', lv));
end

% Bootstrap CI shading
ciPatch = fill([bootCI(1) bootCI(2) bootCI(2) bootCI(1)], ...
    [0 0 60 60], ...
    [0.7 0.3 0.3], 'FaceAlpha', 0.10, 'EdgeColor', 'none', ...
    'DisplayName', sprintf('Crossover 95%% CI [%.1f, %.1f] s', bootCI(1), bootCI(2)));

% Crossover median line
crossLineH = plot([bootMedian bootMedian], [0 60], '--', ...
    'Color', [0.5 0.1 0.1], 'LineWidth', 2, ...
    'DisplayName', sprintf('PAS 4/3 crossover = %.1f s', bootMedian));
text(bootMedian + 1, 55, sprintf('%.1f s', bootMedian), ...
    'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.5 0.1 0.1]);

% ROI trough markers at top
roiTroughs = XOVR.roiTroughs;
yTop     = 57;
yLevels  = [yTop, yTop-3, yTop, yTop-3, yTop-6, yTop-3, yTop];

for r = 1:nROI
    plot(roiTroughs(r), yLevels(r), 'v', 'MarkerSize', 8, ...
        'MarkerFaceColor', roi_colors(r,:), 'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.6, 'HandleVisibility', 'off');
    text(roiTroughs(r), yLevels(r) + 2.0, roiShort{r}, ...
        'FontSize', 7, 'HorizontalAlignment', 'center', ...
        'Color', roi_colors(r,:), 'FontWeight', 'bold');
end

% Dummy for ROI legend entry
roiDummy = plot(NaN, NaN, 'kv', 'MarkerSize', 7, ...
    'MarkerFaceColor', [0.4 0.6 0.8], ...
    'DisplayName', 'ROI sensitivity troughs');

hold off;
xlim([0 T]); ylim([0 62]);
xlabel('Trial time (s)', 'FontSize', font_label);
ylabel('Proportion of clicks (%)', 'FontSize', font_label);
title('(B)  PAS ratings over trial time with crossover CI and ROI troughs', ...
    'FontSize', font_title, 'FontWeight', 'bold');

lgd = legend([pasHandles; crossLineH; ciPatch; roiDummy], ...
    'Location', 'southeast', 'Box', 'on', 'FontSize', 8, ...
    'NumColumns', 2);
set(lgd, 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);
set(gca, 'TickDir', 'out', 'FontSize', font_main, 'Box', 'off');

%% ========================================================================
%  5. SAVE
%  ========================================================================

outFile = '../../results/Figure3_Perceptual.png';
exportgraphics(fig, outFile, 'Resolution', 600);
fprintf('  Saved: %s (600 dpi)\n', outFile);

%% ========================================================================
%  6. SUMMARY
%  ========================================================================

fprintf('\n==========================================================\n');
fprintf('  FIGURE 3 SUMMARY (Perceptual Evidence)\n');
fprintf('==========================================================\n');
fprintf('  Panel A — Sensitivity Schematic\n');
fprintf('    R(x) = exp(-lambda*x), trough at x = 1/e\n');
fprintf('    lambda = e (fixed, not fitted)\n');
fprintf('  --\n');
fprintf('  Panel B — PAS Dynamics + Bootstrap Crossover CI\n');
fprintf('    Bootstrap median crossover: %.1f s\n', bootMedian);
fprintf('    95%% CI: [%.1f, %.1f] s\n', bootCI(1), bootCI(2));
fprintf('    Valid bootstrap samples: %d / %d\n', sum(validBoot), nBoot);
fprintf('==========================================================\n');
fprintf('Done.\n');
