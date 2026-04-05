%% plotFigure1_Behavioral.m
% =========================================================================
% PNAS Figure 1: Behavioral Evidence for Sensitivity Dynamics
% =========================================================================
%
% Two-panel figure combining behavioral evidence:
%
%   Panel A — Click Response-Time Distribution
%     Click-response times follow the sensitivity function:
%       S(x) = A * x * exp(-e * x)
%     where x = (t - tau) / (T - tau) is normalised trial time,
%     lambda = e (Euler's number) is fixed, and tau is onset lag.
%     This is |dP(0)/dlambda| evaluated at lambda = e: the sensitivity
%     of the exponential decay P(0) = exp(-lambda*x) to rate perturbations.
%     A is estimated by scalar projection; tau by grid search.
%
%   Panel B — Haptic Contact Proportion
%     The proportion of participant-trials with active haptic feedback
%     saturates at ~1/e, confirming that the interaction statistics
%     remain stable throughout the trial. The neural and perceptual
%     shifts reported in Figures 2-3 therefore reflect changes in
%     sensitivity to perturbations, not changes in the perturbations
%     themselves.
%
% INPUT:
%   ../../data/ClickTimes/ClickResponseTimes.csv
%   ../../data/Haptics/HapticFeedback.csv.gz  (auto-decompressed)
%
% OUTPUT:
%   ../../results/Figure1_Behavioral.png  (600 dpi)
%
% DEPENDENCIES: Signal Processing Toolbox, Optimization Toolbox
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   April 2026
% =========================================================================

clear; clc; close all;

%% ========================================================================
%  0. CHECK DATA AVAILABILITY
%  ========================================================================

clickFile   = '../../data/ClickTimes/ClickResponseTimes.csv';
hapticFile  = '../../data/Haptics/HapticFeedback.csv';
hapticGz    = '../../data/Haptics/HapticFeedback.csv.gz';

if ~isfile(clickFile)
    error('Missing: %s\nSee README for data download instructions.', clickFile);
end

% Auto-decompress haptic data if only .gz is present
if ~isfile(hapticFile) && isfile(hapticGz)
    fprintf('Decompressing %s ...\n', hapticGz);
    gunzip(hapticGz, fileparts(hapticFile));
end
if ~isfile(hapticFile)
    error('Missing: %s\nSee README for data download instructions.', hapticFile);
end

%% ========================================================================
%  PANEL A — CLICK RESPONSE-TIME DISTRIBUTION WITH SENSITIVITY FIT
%  ========================================================================

fprintf('=== Panel A: Click Response Times ===\n');

%% A1. Load and prepare click data
data_click = readtable(clickFile);
clickTimes = data_click.ClickTime_s(data_click.Clicked == 1);
clickTimes = clickTimes(~isnan(clickTimes));
clickTimes = clickTimes(clickTimes >= 0 & clickTimes < 60);

T        = 60;           % Trial duration (s)
lambda   = exp(1);       % Rate parameter fixed at e
n_clicks = length(clickTimes);
n_dyads  = length(unique(data_click.DyadID));

fprintf('  %d valid clicks from %d dyads\n', n_clicks, n_dyads);

%% A2. Kernel density estimate
n_kde  = 500;
xi_kde = linspace(0, 60, n_kde);
[f_kde, xi_kde] = ksdensity(clickTimes, xi_kde, ...
    'Support', [-0.001 60.001], 'BoundaryCorrection', 'reflection');

%% A3. Histogram (for visual reference only)
bin_width   = 2;
edges       = 0:bin_width:T;
counts      = histcounts(clickTimes, edges);
bin_centres = edges(1:end-1) + bin_width / 2;

%% A4. Grid search: sweep onset lags and fit sensitivity function to KDE
fprintf('  Fitting sensitivity model S(x) = A*x*exp(-e*x) ...\n');

offsets = 0:0.1:15;
n_off   = length(offsets);
R2_vals = zeros(1, n_off);

for k = 1:n_off
    res        = fit_sensitivity_kde(offsets(k), xi_kde, f_kde, T, lambda);
    R2_vals(k) = res.R2;
end

[~, best_idx] = max(R2_vals);
best = fit_sensitivity_kde(offsets(best_idx), xi_kde, f_kde, T, lambda);

fprintf('  Best onset lag: %.1f s,  R^2 = %.4f,  peak = %.1f s\n', ...
    best.offset, best.R2, best.peak_time);

%% ========================================================================
%  PANEL B — HAPTIC CONTACT PROPORTION
%  ========================================================================

fprintf('\n=== Panel B: Haptic Feedback ===\n');

%% B1. Load haptic data
data_haptic = readtable(hapticFile, 'VariableNamingRule', 'preserve');
fprintf('  %d rows loaded.\n', height(data_haptic));

%% B2. Build time × segment matrix
segments = unique(data_haptic(:, {'DyadID','ParticipantID','TrialNum'}), 'rows');
nSeg = height(segments);

firstIdx = data_haptic.DyadID == segments.DyadID(1) & ...
           data_haptic.ParticipantID == segments.ParticipantID(1) & ...
           data_haptic.TrialNum == segments.TrialNum(1);
timeVec = data_haptic.Time_s(firstIdx);
nTime   = length(timeVec);
fs      = round(1 / median(diff(timeVec)));

fprintf('  %d segments, %d time points (%.1f s at %d Hz)\n', ...
    nSeg, nTime, timeVec(end), fs);

sigMatrix = NaN(nTime, nSeg);
for s = 1:nSeg
    idx = data_haptic.DyadID == segments.DyadID(s) & ...
          data_haptic.ParticipantID == segments.ParticipantID(s) & ...
          data_haptic.TrialNum == segments.TrialNum(s);
    sig = data_haptic.HapticFeedback(idx);
    nSamp = min(length(sig), nTime);
    sigMatrix(1:nSamp, s) = sig(1:nSamp);
end

rawBinMatrix = double(sigMatrix >= 0.5);

%% B3. Compute proportion and fit exponential rise
smoothWin = 5;   % seconds
nValid    = sum(~isnan(rawBinMatrix), 2);
propON    = sum(rawBinMatrix, 2, 'omitnan') ./ nValid;
semON     = sqrt(propON .* (1 - propON) ./ nValid);

smoothSamp = round(smoothWin * fs);
propSmooth = smoothdata(propON, 'gaussian', smoothSamp);

expRise = @(p, t) p(1) .* max(1 - exp(-(t - p(3)) ./ p(2)), 0);
p0   = [max(propSmooth), 10, 1];
lb   = [0.01, 0.5, 0];
ub   = [1.0,  60, 15];
opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
    'MaxFunctionEvaluations', 5000, 'MaxIterations', 1000);

[pFit, ~, residuals_h] = ...
    lsqcurvefit(expRise, p0, timeVec, propSmooth, lb, ub, opts);

haptic_A   = pFit(1);
haptic_tau = pFit(2);
haptic_t0  = pFit(3);
SS_res_h   = sum(residuals_h.^2);
SS_tot_h   = sum((propSmooth - mean(propSmooth)).^2);
haptic_R2  = 1 - SS_res_h / SS_tot_h;

target_1e = 1 / exp(1);   % 1/e ~ 0.3679

fprintf('  A = %.4f (1/e = %.4f), tau = %.2f s, R^2 = %.4f\n', ...
    haptic_A, target_1e, haptic_tau, haptic_R2);

%% ========================================================================
%  CREATE FIGURE — 2 panels side by side
%  ========================================================================

fprintf('\nCreating Figure 1 ...\n');

% --- Colour palette ---
col_data  = [0.20 0.40 0.73];    % Blue
col_kde   = [0.12 0.30 0.55];    % Dark blue
col_fit   = [0.80 0.15 0.15];    % Red
col_grey  = [0.50 0.50 0.50];    % Grey

font_sz       = 10;
font_sz_label = 11;
font_sz_title = 12;
font_sz_annot = 9;
font_sz_panel = 16;

% Figure dimensions (inches) — PNAS single-column = 3.42 in, double = 7.09 in
fig_w = 7.09;  fig_h = 3.2;
fig = figure('Units', 'inches', 'Position', [0.5 0.5 fig_w fig_h], ...
    'Color', 'w', 'PaperUnits', 'inches', ...
    'PaperSize', [fig_w fig_h], 'PaperPosition', [0 0 fig_w fig_h]);

%% ---- Panel A: Click-time distribution + sensitivity fit ----------------

ax_a = axes('Position', [0.08 0.17 0.40 0.72]);
hold on;

% Scale KDE to click counts for overlay
kde_scale = n_clicks * bin_width;

% Histogram
bar(bin_centres, counts, 1, ...
    'FaceColor', col_data, 'FaceAlpha', 0.30, ...
    'EdgeColor', 'w', 'LineWidth', 0.5);

% KDE
plot(xi_kde, f_kde * kde_scale, '-', ...
    'Color', col_kde, 'LineWidth', 1.5);

% Sensitivity model curve
t_mod = linspace(best.offset, T, 500);
x_mod = (t_mod - best.offset) / best.T_eff;
x_mod = max(x_mod, 1e-12);
y_mod = (best.A .* x_mod .* exp(-lambda .* x_mod)) * kde_scale;
plot(t_mod, y_mod, '-', 'Color', col_fit, 'LineWidth', 2);

% Peak marker
xline(best.peak_time, ':', 'Color', col_fit, 'LineWidth', 1, 'Alpha', 0.7);

xlim([0 T]);
ylim([0 max(counts) * 1.15]);
xlabel('Time (s)', 'FontSize', font_sz_label);
ylabel('Number of clicks', 'FontSize', font_sz_label);
title('Click response times', 'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'FontSize', font_sz, 'Box', 'on', 'TickDir', 'out');

legend({'Clicks', 'KDE', ...
        sprintf('S(x),  R^2 = %.3f', best.R2)}, ...
    'FontSize', font_sz_annot, 'Location', 'northeast', 'Box', 'off');

text(-0.14, 1.05, 'A', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');
hold off;

%% ---- Panel B: Haptic contact proportion --------------------------------

ax_b = axes('Position', [0.57 0.17 0.40 0.72]);
hold on;

% Raw proportion (faint)
plot(timeVec, propON, 'Color', [col_data 0.20], 'LineWidth', 0.5, ...
    'HandleVisibility', 'off');

% SEM band
fill([timeVec; flipud(timeVec)], ...
     [propSmooth + semON; flipud(propSmooth - semON)], ...
     col_data, 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
     'HandleVisibility', 'off');

% Smoothed line
plot(timeVec, propSmooth, 'Color', col_data, 'LineWidth', 1.5, ...
    'DisplayName', 'Smoothed proportion');

% Exponential fit
tPlot = linspace(0, 60, 1000)';
pPlot = expRise(pFit, tPlot);
plot(tPlot, pPlot, '-', 'Color', col_fit, 'LineWidth', 2, ...
    'DisplayName', sprintf('Exp. rise,  R^2 = %.3f', haptic_R2));

% Asymptote
yline(haptic_A, '--', 'Color', col_fit, 'LineWidth', 1, 'Alpha', 0.6, ...
    'HandleVisibility', 'off');

% 1/e reference line
yline(target_1e, ':', 'Color', col_grey, 'LineWidth', 1.2, 'Alpha', 0.6, ...
    'HandleVisibility', 'off');
text(57, target_1e + 0.012, '1/e', ...
    'FontSize', font_sz_annot, 'Color', col_grey, ...
    'HorizontalAlignment', 'right');

hold off;

xlim([0 60]);
ylim([0 max(propON + semON) * 1.1]);
xlabel('Time (s)', 'FontSize', font_sz_label);
ylabel('Proportion with active contact', 'FontSize', font_sz_label);
title('Haptic feedback activation', 'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'FontSize', font_sz, 'Box', 'on', 'TickDir', 'out');

legend('Location', 'southeast', 'Box', 'off', 'FontSize', font_sz_annot);

% Annotation: asymptote value
text(0.03, 0.95, sprintf('A = %.3f', haptic_A), ...
    'Units', 'normalized', 'FontSize', font_sz_annot, ...
    'VerticalAlignment', 'top', 'Color', col_fit);

text(-0.14, 1.05, 'B', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

%% ========================================================================
%  SAVE
%  ========================================================================

outFile = '../../results/Figure1_Behavioral.png';
print(fig, outFile, '-dpng', '-r600');
fprintf('  Saved: %s (600 dpi)\n', outFile);

%% ========================================================================
%  SUMMARY
%  ========================================================================

fprintf('\n==========================================================\n');
fprintf('  FIGURE 1 SUMMARY\n');
fprintf('==========================================================\n');
fprintf('  Panel A — Click Response Times\n');
fprintf('    Valid clicks:     %d (from %d dyads)\n', n_clicks, n_dyads);
fprintf('    Sensitivity fit:  S(x) = A*x*exp(-e*x)\n');
fprintf('    Best onset lag:   %.1f s\n', best.offset);
fprintf('    Model peak:       %.1f s\n', best.peak_time);
fprintf('    R^2:              %.4f\n', best.R2);
fprintf('    RMSE:             %.2e\n', best.RMSE);
fprintf('  --\n');
fprintf('  Panel B — Haptic Feedback\n');
fprintf('    Segments:         %d participant-trials\n', nSeg);
fprintf('    Exp. rise fit:    p(t) = A*(1 - exp(-(t-t0)/tau))\n');
fprintf('    Asymptote (A):    %.4f  (1/e = %.4f)\n', haptic_A, target_1e);
fprintf('    Time constant:    %.2f s\n', haptic_tau);
fprintf('    R^2:              %.4f\n', haptic_R2);
fprintf('==========================================================\n');
fprintf('Done.\n');

%% ========================================================================
%  LOCAL FUNCTION
%  ========================================================================

function res = fit_sensitivity_kde(offset, xi_kde, f_kde, T, lambda)
% FIT_SENSITIVITY_KDE  Fit sensitivity function to KDE density.
%
%   Model:  f(t) = A * x * exp(-lambda * x)
%           x    = (t - offset) / (T - offset)
%
%   This is |dP(0)/dlambda| = x * exp(-lambda*x), the sensitivity
%   of the exponential decay to rate perturbations, scaled by A.
%
%   A is estimated by scalar projection (closed-form OLS with one
%   regressor).  Lambda is fixed at e.  Onset lag (offset) is
%   selected externally by grid search.

    T_eff = T - offset;
    mask  = xi_kde > offset;
    t_fit = xi_kde(mask);
    f_fit = f_kde(mask);

    if sum(mask) < 10
        res.offset = offset; res.A = 0; res.R2 = 0;
        res.RMSE = Inf; res.T_eff = T_eff; res.peak_time = NaN;
        res.t_fit = []; res.f_fit = []; res.y_fitted = [];
        return;
    end

    % Normalised time
    x = (t_fit - offset) / T_eff;
    x = max(x, eps);

    % Sensitivity shape: x * exp(-lambda * x)
    y_shape = x .* exp(-lambda .* x);

    % Scalar projection: A = (y_shape . f) / (y_shape . y_shape)
    A        = dot(y_shape(:), f_fit(:)) / dot(y_shape(:), y_shape(:));
    y_fitted = A .* y_shape;

    % Goodness of fit
    SS_res = sum((f_fit(:) - y_fitted(:)).^2);
    SS_tot = sum((f_fit(:) - mean(f_fit)).^2);

    res.offset    = offset;
    res.A         = A;
    res.R2        = max(1 - SS_res / SS_tot, 0);
    res.RMSE      = sqrt(mean((f_fit(:) - y_fitted(:)).^2));
    res.T_eff     = T_eff;
    res.peak_time = offset + T_eff / lambda;
    res.t_fit     = t_fit;
    res.f_fit     = f_fit;
    res.y_fitted  = y_fitted;
end
