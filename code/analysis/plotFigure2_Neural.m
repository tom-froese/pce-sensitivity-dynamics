%% plotFigure2_Neural.m
% =========================================================================
% PNAS Figure 2: Neural Evidence — Sensitivity Model Fits by Brain Region
% =========================================================================
%
% Single-panel figure showing 7-ROI sensitivity model fits overlaid with
% a continuous anterior-to-posterior color gradient (dark blue → light blue).
%
% Each ROI's grand-median GSP time course is shown (faint) alongside its
% sensitivity model fit S(x) = A*x*exp(-e*x) + B (solid). Trough markers
% (triangles) indicate the time of maximum sensitivity deflection for each
% region, revealing an anterior-to-posterior gradient in timing.
%
% Cross-figure reference lines are overlaid:
%   - Behavioral S(x) peak from click response times (Figure 1)
%   - PAS 4->3 perceptual crossover time (Figure 3)
%
% INPUT:
%   ../../data/preprocessed/EEG/globalScalpPotential_stats.mat  (from computeGSPStats.m)
%   ../../data/preprocessed/ClickTimes/ClickResponseTimes.csv   (for Fig 1 peak)
%   ../../data/preprocessed/PAS/gsp_pas_roi_crossover.mat       (for Fig 3 crossover)
%
% OUTPUT:
%   ../../results/Figure2_Neural.png  (600 dpi)
%
% DEPENDENCIES: Statistics and Machine Learning Toolbox
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   April 2026
% =========================================================================

clear; clc; close all;

%% ========================================================================
%  1. LOAD PRECOMPUTED DATA
%  ========================================================================

statsFile = '../../data/preprocessed/EEG/globalScalpPotential_stats.mat';

if ~isfile(statsFile)
    error(['Missing: %s\n' ...
           'Run computeGSPStats.m first, or download from Zenodo.\n' ...
           'See README for instructions.'], statsFile);
end

S = load(statsFile);

tTaskSm = S.tTaskSm;
nPart   = S.nPart;
roi     = S.roi;
nROI    = S.nROI;

fprintf('Loaded GSP data: N=%d participants, %d ROIs\n', nPart, nROI);

%% ========================================================================
%  1b. CROSS-FIGURE REFERENCE VALUES (Figs 1 & 3)
%  ========================================================================

% --- Behavioral S(x) peak from Figure 1 (click response times) ---
clickFile = '../../data/preprocessed/ClickTimes/ClickResponseTimes.csv';
if isfile(clickFile)
    data_click = readtable(clickFile);
    if ismember('ClickTime_s', data_click.Properties.VariableNames)
        clickTimes = data_click.ClickTime_s(data_click.Clicked == 1);
    else
        clickTimes = data_click.ResponseTime_s(data_click.Clicked == 1);
    end
    clickTimes = clickTimes(~isnan(clickTimes) & clickTimes >= 0 & clickTimes < 60);

    T_trial = 60;  lambda_e = exp(1);
    xi_click = linspace(0, 60, 500);
    [f_click, xi_click] = ksdensity(clickTimes, xi_click, ...
        'Support', [-0.001 60.001], 'BoundaryCorrection', 'reflection');
    offsets = 0:0.1:15;
    R2_click = zeros(1, length(offsets));
    for k = 1:length(offsets)
        res = fit_sensitivity_kde(offsets(k), xi_click, f_click, T_trial, lambda_e);
        R2_click(k) = res.R2;
    end
    [~, best_idx] = max(R2_click);
    bestClick = fit_sensitivity_kde(offsets(best_idx), xi_click, f_click, T_trial, lambda_e);
    clickPeakTime = bestClick.peak_time;
    fprintf('  Behavioral S(x) peak (Fig 1): %.1f s\n', clickPeakTime);
else
    warning('Click data not found — skipping behavioral peak line.');
    clickPeakTime = NaN;
end

% --- PAS 4->3 crossover from Figure 3 (bootstrap median) ---
crossoverFile = '../../data/preprocessed/PAS/gsp_pas_roi_crossover.mat';
if isfile(crossoverFile)
    XOVR = load(crossoverFile);
    winCenters_xovr = XOVR.winCenters;
    nWin_xovr       = length(winCenters_xovr);
    winProp3 = XOVR.winProp{3};
    winProp4 = XOVR.winProp{4};
    nPall    = size(winProp3, 1);

    smKernel = ones(1,3)/3;
    nBoot = 10000;  rng(42);
    bootCross = nan(nBoot, 1);
    minCrossTime = 15;
    minIdx_xovr  = find(winCenters_xovr >= minCrossTime, 1);
    for bi = 1:nBoot
        idx = randi(nPall, nPall, 1);
        gP3 = mean(winProp3(idx,:), 1, 'omitnan');
        gP4 = mean(winProp4(idx,:), 1, 'omitnan');
        diff43   = gP4 - gP3;
        diff43sm = conv(diff43, smKernel, 'same');
        for ci = minIdx_xovr:nWin_xovr-1
            if diff43sm(ci) > 0 && diff43sm(ci+1) <= 0
                x1 = winCenters_xovr(ci);  x2 = winCenters_xovr(ci+1);
                y1 = diff43sm(ci);          y2 = diff43sm(ci+1);
                bootCross(bi) = x1 + (0 - y1) * (x2 - x1) / (y2 - y1);
                break;
            end
        end
    end
    pasCrossoverTime = median(bootCross(~isnan(bootCross)));
    fprintf('  PAS 4->3 crossover (Fig 3): %.1f s\n', pasCrossoverTime);
else
    warning('PAS crossover data not found — skipping crossover line.');
    pasCrossoverTime = NaN;
end

%% ========================================================================
%  2. STYLING — continuous anterior-to-posterior gradient
%  ========================================================================

baseColor  = [0.10 0.35 0.65];   % Dark blue (anterior)
lightColor = [0.70 0.85 0.95];   % Light blue (posterior)
roi_colors = zeros(nROI, 3);
for r = 1:nROI
    frac = (r - 1) / (nROI - 1);
    roi_colors(r,:) = baseColor * (1 - frac) + lightColor * frac;
end

font_main  = 11;
font_label = 12;
font_title = 13;

%% ========================================================================
%  3. CREATE FIGURE — single panel
%  ========================================================================

fprintf('Creating Figure 2 ...\n');

fig_w = 7.5;  fig_h = 3.5;
fig = figure('Units', 'inches', 'Position', [0.5 0.5 fig_w fig_h], 'Color', 'w', ...
    'PaperUnits', 'inches', 'PaperSize', [fig_w fig_h], 'PaperPosition', [0 0 fig_w fig_h]);

% Explicit axes: leave 40% of figure width for legend + annotation on right
ax = axes(fig, 'Position', [0.08 0.15 0.52 0.78]);
hold on;

lineHandles = gobjects(nROI, 1);
for ri = 1:nROI
    col = roi_colors(ri,:);

    % Smoothed grand median (faint)
    plot(tTaskSm, S.roiGrandSm(ri,:), '-', 'Color', [col 0.3], 'LineWidth', 0.8, ...
        'HandleVisibility', 'off');

    % Sensitivity fit (solid, labeled with trough as T_eff/e)
    lineHandles(ri) = plot(tTaskSm, S.roiGrandFit(ri,:), '-', 'Color', col, ...
        'LineWidth', 2.2, ...
        'DisplayName', sprintf('%.1f s   %s', ...
            S.roiGrandTrough(ri), roi.names{ri}));

    % Trough marker (minimum of fit, A < 0)
    [trVal, trIdx] = min(S.roiGrandFit(ri,:));
    plot(tTaskSm(trIdx), trVal, 'v', 'MarkerSize', 6, ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, ...
        'HandleVisibility', 'off');
end

% --- Cross-figure reference lines (Figs 1 & 3) ---
yl = ylim;
yLabel = yl(2) - 0.02 * (yl(2) - yl(1));
if ~isnan(clickPeakTime)
    xline(clickPeakTime, ':', 'Color', [0.80 0.15 0.15], ...
        'LineWidth', 1.5, 'Alpha', 0.7, 'HandleVisibility', 'off');
    text(clickPeakTime - 0.5, yLabel, {'Click peak', sprintf('%.1f s', clickPeakTime)}, ...
        'FontSize', 7, 'Color', [0.80 0.15 0.15], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
end
if ~isnan(pasCrossoverTime)
    xline(pasCrossoverTime, '--', 'Color', [0.50 0.10 0.10], ...
        'LineWidth', 1.5, 'Alpha', 0.7, 'HandleVisibility', 'off');
    text(pasCrossoverTime + 0.5, yLabel, {'PAS 4\rightarrow3', sprintf('%.1f s', pasCrossoverTime)}, ...
        'FontSize', 7, 'Color', [0.50 0.10 0.10], ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end

yline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, 'HandleVisibility', 'off');
hold off;

xlabel('Time (s)', 'FontSize', font_label);
ylabel(['Amplitude (' char(181) 'V)'], 'FontSize', font_label);
title(sprintf('Sensitivity model fits by brain region  (N = %d)', nPart), ...
    'FontSize', font_title, 'FontWeight', 'bold');
xlim([tTaskSm(1) tTaskSm(end)]);
box off;
set(gca, 'TickDir', 'out', 'FontSize', font_main);

% --- Legend (right panel, upper portion) ---
lgd = legend(lineHandles, 'Box', 'off', 'FontSize', 8);
title(lgd, 'T_{eff}/e   (Anterior \rightarrow Posterior)');
lgd.Units = 'normalized';
lgd.Position = [0.62, 0.48, 0.36, 0.48];

% --- Model annotation (right panel, below legend) ---
tauRange = [min(S.roiGrandTau), max(S.roiGrandTau)];
R2range  = [min(S.roiGrandR2), max(S.roiGrandR2)];
annotation(fig, 'textbox', [0.62, 0.04, 0.36, 0.42], ...
    'String', {
        '\bfModel:\rm  S(x) = A \cdot x \cdot exp(-ex) + B', ...
        '', ...
        '\lambda = e  fixed;  A, B  fitted (OLS)', ...
        '\tau  by grid search  (onset lag)', ...
        'Trough at x = 1/e,  i.e.  t = \tau + T_{eff}/e', ...
        '', ...
        sprintf('\\tau range:   %.1f – %.1f s', tauRange(1), tauRange(2)), ...
        sprintf('R^2 range:  %.2f – %.2f', R2range(1), R2range(2))
    }, ...
    'FontSize', 8, 'EdgeColor', [0.7 0.7 0.7], ...
    'BackgroundColor', 'w', 'Margin', 5, ...
    'VerticalAlignment', 'top', 'FitBoxToText', 'off', ...
    'Interpreter', 'tex');

%% ========================================================================
%  4. SAVE
%  ========================================================================

outFile = '../../results/Figure2_Neural.png';
exportgraphics(fig, outFile, 'Resolution', 600);
fprintf('  Saved: %s (600 dpi)\n', outFile);

%% ========================================================================
%  5. SUMMARY
%  ========================================================================

fprintf('\n==========================================================\n');
fprintf('  FIGURE 2 SUMMARY (Neural Evidence)\n');
fprintf('==========================================================\n');
fprintf('  N participants: %d\n', nPart);
fprintf('  ROI sensitivity fits (anterior -> posterior):\n');
for ri = 1:nROI
    fprintf('    %-18s  R^2 = %.3f,  trough = %.1f s\n', ...
        roi.names{ri}, S.roiGrandR2(ri), S.roiGrandTrough(ri));
end
fprintf('  Cross-figure references:\n');
if ~isnan(clickPeakTime)
    fprintf('    Click S(x) peak (Fig 1):    %.1f s\n', clickPeakTime);
end
if ~isnan(pasCrossoverTime)
    fprintf('    PAS 4->3 crossover (Fig 3): %.1f s\n', pasCrossoverTime);
end
fprintf('==========================================================\n');
fprintf('Done.\n');

%% ========================================================================
%  LOCAL FUNCTION
%  ========================================================================

function res = fit_sensitivity_kde(offset, xi_kde, f_kde, T, lambda)
% Fit sensitivity function S(x) = A * x * exp(-lambda * x) to KDE density.
% x = (t - offset) / (T - offset),  lambda fixed at e.

    T_eff = T - offset;
    mask  = xi_kde > offset;
    t_fit = xi_kde(mask);
    f_fit = f_kde(mask);

    if sum(mask) < 10
        res.offset = offset; res.A = 0; res.R2 = 0;
        res.RMSE = Inf; res.T_eff = T_eff; res.peak_time = NaN;
        return;
    end

    x = (t_fit - offset) / T_eff;
    x = max(x, eps);
    y_shape = x .* exp(-lambda .* x);

    A = dot(y_shape(:), f_fit(:)) / dot(y_shape(:), y_shape(:));
    y_fitted = A .* y_shape;

    SS_res = sum((f_fit(:) - y_fitted(:)).^2);
    SS_tot = sum((f_fit(:) - mean(f_fit)).^2);

    res.offset    = offset;
    res.A         = A;
    res.R2        = max(1 - SS_res / SS_tot, 0);
    res.RMSE      = sqrt(mean((f_fit(:) - y_fitted(:)).^2));
    res.T_eff     = T_eff;
    res.peak_time = offset + T_eff / lambda;
end
