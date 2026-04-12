%% plotFigure1_Behavioral.m
% =========================================================================
% PNAS Figure 1: Bodily Evidence for Reliability Decay Dynamics
% =========================================================================
%
% Three-panel figure combining bodily/behavioral evidence:
%
%   Panel A — Click Response-Time Distribution
%     Click-response times follow the sensitivity function:
%       S(x) = A * x * exp(-e * x)
%     where x = (t - tau) / (T - tau) is normalised trial time,
%     lambda = e (Euler's number) is fixed, and tau is onset lag.
%     This is |dR/dlambda| evaluated at lambda = e: the sensitivity
%     of the reliability R(x) = exp(-lambda*x) to rate perturbations.
%
%   Panel B — Haptic Contact Proportion
%     The proportion of participant-trials with active haptic feedback
%     saturates at ~1/e, confirming that the interaction statistics
%     remain stable throughout the trial.
%
%   Panel C — Electrodermal Activity (EDA) — Reliability Decay
%     Sympathetic arousal (skin conductance) follows the reliability
%     function R(x) = A0*exp(-e*x) + B, with lambda fixed at e.
%     Provides converging autonomic evidence for the exponential decay.
%
% INPUT:
%   ../../data/preprocessed/ClickTimes/ClickResponseTimes.csv
%   ../../data/preprocessed/Haptics/HapticFeedback.csv.gz  (auto-decompressed)
%   ../../data/preprocessed/EDA/EDA_Task_Preprocessed.csv
%   ../../data/preprocessed/EDA/EDA_Rest_Preprocessed.csv
%   ../../data/preprocessed/EEG/globalScalpPotential_stats.mat  (for dissociation)
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
%  0. DATA PATHS AND AVAILABILITY CHECK
%  ========================================================================

clickFile   = '../../data/preprocessed/ClickTimes/ClickResponseTimes.csv';
hapticFile  = '../../data/preprocessed/Haptics/HapticFeedback.csv';
hapticGz    = '../../data/preprocessed/Haptics/HapticFeedback.csv.gz';
edaTaskFile = '../../data/preprocessed/EDA/EDA_Task_Preprocessed.csv';
edaRestFile = '../../data/preprocessed/EDA/EDA_Rest_Preprocessed.csv';
gspStatsFile = '../../data/preprocessed/EEG/globalScalpPotential_stats.mat';

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

if ~isfile(edaTaskFile)
    error('Missing: %s\nRun preprocessEDA.m first, or see README.', edaTaskFile);
end

%% ========================================================================
%  SHARED PARAMETERS
%  ========================================================================

T_trial = 60;          % Trial duration (s)
lambda  = exp(1);      % Rate parameter fixed at e

%% ========================================================================
%  PANEL A — CLICK RESPONSE-TIME DISTRIBUTION WITH SENSITIVITY FIT
%  ========================================================================

fprintf('=== Panel A: Click Response Times ===\n');

data_click = readtable(clickFile);
% Handle both column name formats
if ismember('ClickTime_s', data_click.Properties.VariableNames)
    clickTimes = data_click.ClickTime_s(data_click.Clicked == 1);
else
    clickTimes = data_click.ResponseTime_s(data_click.Clicked == 1);
end
clickTimes = clickTimes(~isnan(clickTimes));
clickTimes = clickTimes(clickTimes >= 0 & clickTimes < 60);

n_clicks = length(clickTimes);
fprintf('  %d valid clicks\n', n_clicks);

% KDE
n_kde  = 500;
xi_kde = linspace(0, 60, n_kde);
[f_kde, xi_kde] = ksdensity(clickTimes, xi_kde, ...
    'Support', [-0.001 60.001], 'BoundaryCorrection', 'reflection');

% Histogram
bin_width   = 2;
edges       = 0:bin_width:T_trial;
counts      = histcounts(clickTimes, edges);
bin_centres = edges(1:end-1) + bin_width / 2;

% Grid search: onset lag
fprintf('  Fitting sensitivity model S(x) = A*x*exp(-e*x) ...\n');
offsets = 0:0.1:15;
R2_vals = zeros(1, length(offsets));
for k = 1:length(offsets)
    res = fit_sensitivity_kde(offsets(k), xi_kde, f_kde, T_trial, lambda);
    R2_vals(k) = res.R2;
end
[~, best_idx] = max(R2_vals);
best = fit_sensitivity_kde(offsets(best_idx), xi_kde, f_kde, T_trial, lambda);

fprintf('  Onset lag: %.1f s,  R^2 = %.4f,  peak = %.1f s\n', ...
    best.offset, best.R2, best.peak_time);

%% ========================================================================
%  PANEL B — HAPTIC CONTACT PROPORTION
%  ========================================================================

fprintf('\n=== Panel B: Haptic Feedback ===\n');

data_haptic = readtable(hapticFile, 'VariableNamingRule', 'preserve');
segments = unique(data_haptic(:, {'DyadID','ParticipantID','TrialNum'}), 'rows');
nSeg = height(segments);

firstIdx = data_haptic.DyadID == segments.DyadID(1) & ...
           data_haptic.ParticipantID == segments.ParticipantID(1) & ...
           data_haptic.TrialNum == segments.TrialNum(1);
timeVec_h = data_haptic.Time_s(firstIdx);
nTime_h   = length(timeVec_h);
fs_h      = round(1 / median(diff(timeVec_h)));

fprintf('  %d segments, %d time points (%.1f s at %d Hz)\n', ...
    nSeg, nTime_h, timeVec_h(end), fs_h);

sigMatrix = NaN(nTime_h, nSeg);
for s = 1:nSeg
    idx = data_haptic.DyadID == segments.DyadID(s) & ...
          data_haptic.ParticipantID == segments.ParticipantID(s) & ...
          data_haptic.TrialNum == segments.TrialNum(s);
    sig = data_haptic.HapticFeedback(idx);
    nSamp = min(length(sig), nTime_h);
    sigMatrix(1:nSamp, s) = sig(1:nSamp);
end
rawBinMatrix = double(sigMatrix >= 0.5);

smoothWin  = 5;
nValid     = sum(~isnan(rawBinMatrix), 2);
propON     = sum(rawBinMatrix, 2, 'omitnan') ./ nValid;
semON      = sqrt(propON .* (1 - propON) ./ nValid);
smoothSamp = round(smoothWin * fs_h);
propSmooth = smoothdata(propON, 'gaussian', smoothSamp);

expRise = @(p, t) p(1) .* max(1 - exp(-(t - p(3)) ./ p(2)), 0);
p0   = [max(propSmooth), 10, 1];
lb   = [0.01, 0.5, 0];
ub   = [1.0, 60, 15];
opts_h = optimoptions('lsqcurvefit', 'Display', 'off', ...
    'MaxFunctionEvaluations', 5000, 'MaxIterations', 1000);
[pFit, ~, residuals_h] = lsqcurvefit(expRise, p0, timeVec_h, propSmooth, lb, ub, opts_h);

haptic_A  = pFit(1);
haptic_R2 = 1 - sum(residuals_h.^2) / sum((propSmooth - mean(propSmooth)).^2);
target_1e = 1 / exp(1);

fprintf('  A = %.4f (1/e = %.4f),  R^2 = %.4f\n', haptic_A, target_1e, haptic_R2);

%% ========================================================================
%  PANEL C — EDA RELIABILITY DECAY
%  ========================================================================

fprintf('\n=== Panel C: EDA Reliability Decay ===\n');

fs_eda       = 25;
task_samples = T_trial * fs_eda;

T_eda = readtable(edaTaskFile);
T_eda.PID = T_eda.DyadID * 10 + T_eda.ParticipantID;
pids = unique(T_eda.PID);

part_matrix = NaN(length(pids), task_samples);
for i = 1:length(pids)
    pid = pids(i);
    pdata = T_eda(T_eda.PID == pid, :);
    trials = unique(pdata.TrialNum);
    trial_sum = zeros(1, task_samples); trial_count = 0;
    for t = 1:length(trials)
        eda = pdata.EDA_uS(pdata.TrialNum == trials(t));
        if length(eda) == task_samples
            trial_sum = trial_sum + eda'; trial_count = trial_count + 1;
        end
    end
    if trial_count > 0
        part_matrix(i, :) = trial_sum / trial_count;
    end
end
part_matrix = part_matrix(~any(isnan(part_matrix), 2), :);
nPart_eda = size(part_matrix, 1);
time_eda  = (0:task_samples-1) / fs_eda;
grand_mean_eda = mean(part_matrix, 1);
grand_sem_eda  = std(part_matrix, 0, 1) / sqrt(nPart_eda);

% Fit R(x) = A0*exp(-e*x) + B with onset lag sweep
model_eda = @(params, x) params(1) * exp(-exp(1) * x) + params(2);
opts_eda = optimoptions('lsqcurvefit', 'Display', 'off', ...
    'MaxIterations', 10000, 'MaxFunctionEvaluations', 30000, ...
    'TolFun', 1e-12, 'TolX', 1e-12);

tau_sweep = 0:0.1:15;
R2_eda_sweep = NaN(1, length(tau_sweep));
A0_eda_sweep = NaN(1, length(tau_sweep));
B_eda_sweep  = NaN(1, length(tau_sweep));

for s = 1:length(tau_sweep)
    tau = tau_sweep(s);
    T_eff = T_trial - tau;
    if T_eff < 10; continue; end
    tau_samp = round(tau * fs_eda);
    idx = (tau_samp + 1):task_samples;
    y = grand_mean_eda(idx);
    x = (time_eda(idx) - tau) / T_eff;
    A0_init = y(1) - y(end); B_init = y(end);
    params = lsqcurvefit(model_eda, [A0_init, B_init], x, y, [0 0], [50 50], opts_eda);
    yfit = model_eda(params, x);
    SS_res = sum((y - yfit).^2);
    SS_tot = sum((y - mean(y)).^2);
    R2_eda_sweep(s) = 1 - SS_res / SS_tot;
    A0_eda_sweep(s) = params(1);
    B_eda_sweep(s)  = params(2);
end

[R2_eda_best, idx_eda_best] = max(R2_eda_sweep);
tau_eda_opt = tau_sweep(idx_eda_best);
T_eff_eda   = T_trial - tau_eda_opt;
A0_eda_best = A0_eda_sweep(idx_eda_best);
B_eda_best  = B_eda_sweep(idx_eda_best);

tau_samp_eda = round(tau_eda_opt * fs_eda);
idx_eda_opt  = (tau_samp_eda + 1):task_samples;
x_eda_opt    = (time_eda(idx_eda_opt) - tau_eda_opt) / T_eff_eda;
yfit_eda_opt = model_eda([A0_eda_best, B_eda_best], x_eda_opt);

fprintf('  N = %d,  tau = %.1f s,  R^2 = %.4f,  A0 = %.3f,  B = %.3f\n', ...
    nPart_eda, tau_eda_opt, R2_eda_best, A0_eda_best, B_eda_best);

%% ---- REST EDA (for dissociation analysis) ----

fprintf('\n  --- Rest EDA ---\n');

rest_duration = 180;
rest_samples  = rest_duration * fs_eda;

if isfile(edaRestFile)
    T_rest_eda  = readtable(edaRestFile);
    T_rest_eda.PID = T_rest_eda.DyadID * 10 + T_rest_eda.ParticipantID;
    rest_pids = unique(T_rest_eda.PID);

    rest_part_matrix = NaN(length(rest_pids), rest_samples);
    for i = 1:length(rest_pids)
        pid = rest_pids(i);
        pdata = T_rest_eda(T_rest_eda.PID == pid, :);
        periods = unique(pdata.RestNum);
        period_sum = zeros(1, rest_samples); period_count = 0;
        for r = 1:length(periods)
            eda = pdata.EDA_uS(pdata.RestNum == periods(r));
            if length(eda) == rest_samples
                period_sum = period_sum + eda'; period_count = period_count + 1;
            end
        end
        if period_count > 0
            rest_part_matrix(i, :) = period_sum / period_count;
        end
    end
    rest_part_matrix = rest_part_matrix(~any(isnan(rest_part_matrix), 2), :);
    nPart_rest_eda = size(rest_part_matrix, 1);
    time_rest_eda  = (0:rest_samples-1) / fs_eda;
    grand_mean_rest_eda = mean(rest_part_matrix, 1);

    % Fit R(x) to rest
    x_rest_eda = time_rest_eda / rest_duration;
    y_rest_eda = grand_mean_rest_eda;
    A0_init_r = y_rest_eda(1) - y_rest_eda(end);
    B_init_r  = y_rest_eda(end);
    params_rest = lsqcurvefit(model_eda, [A0_init_r, B_init_r], x_rest_eda, ...
        y_rest_eda, [0 0], [50 50], opts_eda);
    yfit_rest_eda = model_eda(params_rest, x_rest_eda);
    SS_res_r = sum((y_rest_eda - yfit_rest_eda).^2);
    SS_tot_r = sum((y_rest_eda - mean(y_rest_eda)).^2);
    R2_rest_eda = 1 - SS_res_r / SS_tot_r;

    fprintf('  Rest N = %d,  R^2 = %.4f\n', nPart_rest_eda, R2_rest_eda);
else
    fprintf('  Rest EDA file not found — skipping dissociation analysis.\n');
    R2_rest_eda = NaN;
end

% Load GSP rest R^2 for dissociation comparison
if isfile(gspStatsFile)
    gsp_stats = load(gspStatsFile, 'grandR2', 'grandRestR2');
    fprintf('\n  === DISSOCIATION SUMMARY ===\n');
    fprintf('  EDA  R(x) fit:   Task R^2 = %.3f,  Rest R^2 = %.3f\n', ...
        R2_eda_best, R2_rest_eda);
    fprintf('  GSP  S(x) fit:   Task R^2 = %.3f,  Rest R^2 = %.3f\n', ...
        gsp_stats.grandR2, gsp_stats.grandRestR2);
    fprintf('  -> EDA decay present in both conditions (autonomic clock)\n');
    fprintf('  -> GSP sensitivity only in task (neural engagement requires stakes)\n');
end

%% ========================================================================
%  CREATE FIGURE — 3 panels in a row
%  ========================================================================

fprintf('\nCreating Figure 1 ...\n');

% --- Colour palette ---
col_data  = [0.20 0.40 0.73];    % Blue
col_kde   = [0.12 0.30 0.55];    % Dark blue
col_fit   = [0.80 0.15 0.15];    % Red
col_grey  = [0.50 0.50 0.50];
col_eda   = [0.85 0.32 0.10];    % Orange for EDA
col_black = [0.10 0.10 0.10];    % Black for R(x) fit

font_sz       = 9;
font_sz_label = 10;
font_sz_title = 11;
font_sz_annot = 8;
font_sz_panel = 14;

% PNAS double-column width: 7.09 in
fig_w = 7.5;  fig_h = 3.0;
fig = figure('Units', 'inches', 'Position', [0.5 0.5 fig_w fig_h], ...
    'Color', 'w', 'PaperUnits', 'inches', ...
    'PaperSize', [fig_w fig_h], 'PaperPosition', [0 0 fig_w fig_h]);

% Panel positions
ml = 0.06;  mr = 0.02;  mb = 0.18;  mt = 0.08;
gap_ab = 0.07;  gap_bc = 0.07;
total_w = 1 - ml - mr - gap_ab - gap_bc;
pw_a = total_w * 0.36;
pw_b = total_w * 0.28;
pw_c = total_w * 0.36;
ph = 1 - mb - mt;

%% ---- Panel A: Click-time distribution + sensitivity fit ----------------

ax_a = axes('Position', [ml, mb, pw_a, ph]);
hold on;

kde_scale = n_clicks * bin_width;

bar(bin_centres, counts, 1, ...
    'FaceColor', col_data, 'FaceAlpha', 0.30, ...
    'EdgeColor', 'w', 'LineWidth', 0.5);

plot(xi_kde, f_kde * kde_scale, '-', 'Color', col_kde, 'LineWidth', 1.5);

t_mod = linspace(best.offset, T_trial, 500);
x_mod = (t_mod - best.offset) / best.T_eff;
x_mod = max(x_mod, 1e-12);
y_mod = (best.A .* x_mod .* exp(-lambda .* x_mod)) * kde_scale;
plot(t_mod, y_mod, '-', 'Color', col_fit, 'LineWidth', 2);

% Peak marker (label only — numeric value is in the annotation box)
xline(best.peak_time, ':', 'Color', col_fit, 'LineWidth', 1, 'Alpha', 0.7);
text(best.peak_time + 1, 3, 'T_{eff}/e', ...
    'FontSize', font_sz_annot, 'Color', col_fit);

% Onset lag marker
xline(best.offset, ':', 'Color', col_grey, 'LineWidth', 1.2);

xlim([0 T_trial]);
ylim([0 max(counts) * 1.35]);
xlabel('Time (s)', 'FontSize', font_sz_label);
ylabel('Number of clicks', 'FontSize', font_sz_label);
title('Click response times', 'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'FontSize', font_sz, 'Box', 'off', 'TickDir', 'out');

% Data legend (top-right — white box to prevent overlap with peak line)
lg_a = legend({'Clicks', 'KDE', 'S(x) = x \cdot exp(-ex)'}, ...
    'FontSize', font_sz_annot, 'Location', 'northeast');
lg_a.BoxFace.ColorType = 'truecoloralpha';
lg_a.BoxFace.ColorData = uint8([255; 255; 255; 255]);
lg_a.EdgeColor = [0.7 0.7 0.7];

% Annotation box (bottom-right)
text(0.97, 0.05, ...
    {sprintf('\\tau = %.1f s', best.offset), ...
     sprintf('R^2 = %.3f', best.R2), ...
     sprintf('peak = %.1f s', best.peak_time)}, ...
    'Units', 'normalized', 'FontSize', font_sz_annot, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
    'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'Margin', 2);

text(-0.16, 1.06, 'A', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');
hold off;

%% ---- Panel B: Haptic contact proportion --------------------------------

ax_b = axes('Position', [ml + pw_a + gap_ab, mb, pw_b, ph]);
hold on;

% Plot in percentage (proportion * 100)
plot(timeVec_h, propON * 100, 'Color', [col_data 0.20], 'LineWidth', 0.5, ...
    'HandleVisibility', 'off');
fill([timeVec_h; flipud(timeVec_h)], ...
     [(propSmooth + semON) * 100; flipud((propSmooth - semON) * 100)], ...
     col_data, 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
     'HandleVisibility', 'off');
plot(timeVec_h, propSmooth * 100, 'Color', col_data, 'LineWidth', 1.5, ...
    'DisplayName', 'Smoothed');

tPlot = linspace(0, 60, 1000)';
pPlot = expRise(pFit, tPlot);
plot(tPlot, pPlot * 100, '-', 'Color', col_fit, 'LineWidth', 2, ...
    'DisplayName', sprintf('Exp. rise,  R^2 = %.3f', haptic_R2));

% Asymptote line
yline(haptic_A * 100, '--', 'Color', col_fit, 'LineWidth', 1, 'Alpha', 0.6, ...
    'HandleVisibility', 'off');

% 1/e reference line (red dotted, evoking the 1/e peak of Panel A)
yline(target_1e * 100, ':', 'Color', col_fit, 'LineWidth', 1.2, 'Alpha', 0.6, ...
    'HandleVisibility', 'off');
text(57, target_1e * 100 + 1.2, '1/e', 'FontSize', font_sz_annot, ...
    'Color', col_fit, 'HorizontalAlignment', 'right');

% 95% of asymptote boot-up time: t = onset - tau_h * ln(1 - 0.95)
haptic_onset    = pFit(3);    % onset parameter
haptic_tau_h    = pFit(2);    % time constant
haptic_t95      = haptic_onset - haptic_tau_h * log(1 - 0.95);
fprintf('  Boot-up: onset = %.2f s, tau = %.2f s, t_95%% = %.1f s\n', ...
    haptic_onset, haptic_tau_h, haptic_t95);

% Boot-up marker (grey, consistent with tau markers in other panels)
xline(haptic_t95, ':', 'Color', col_grey, 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');

hold off;
xlim([0 60]);
ylim([0 max((propON + semON) * 100) * 1.1]);
xlabel('Time (s)', 'FontSize', font_sz_label);
ylabel('Trials with active contact (%)', 'FontSize', font_sz_label);
title('Haptic feedback activation', 'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'FontSize', font_sz, 'Box', 'off', 'TickDir', 'out');

% Compact legend (smaller ItemTokenSize to avoid line overlap with curve)
lg_b = legend('Location', 'southeast', 'Box', 'off', 'FontSize', font_sz_annot);
lg_b.ItemTokenSize = [15, 8];

% Boot-up annotation (grey, shifted right to clear the vertical line)
text(0.12, 0.95, sprintf('t_{95%%} = %.1f s', haptic_t95), ...
    'Units', 'normalized', 'FontSize', font_sz_annot, ...
    'VerticalAlignment', 'top', 'Color', col_grey);

text(-0.16, 1.06, 'B', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

%% ---- Panel C: EDA reliability decay ------------------------------------

ax_c = axes('Position', [ml + pw_a + gap_ab + pw_b + gap_bc, mb, pw_c, ph]);
hold on;

% SEM shading
fill([time_eda, fliplr(time_eda)], ...
    [grand_mean_eda + grand_sem_eda, fliplr(grand_mean_eda - grand_sem_eda)], ...
    col_eda, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Grand average
h_eda_data = plot(time_eda, grand_mean_eda, '-', 'Color', col_eda, 'LineWidth', 1.5);

% R(x) fit
h_eda_fit = plot(time_eda(idx_eda_opt), yfit_eda_opt, '--', ...
    'Color', col_black, 'LineWidth', 2);

% Onset lag
xline(tau_eda_opt, ':', 'Color', col_grey, 'LineWidth', 1.2);

hold off;

xlim([0 60]);
xlabel('Time (s)', 'FontSize', font_sz_label);
ylabel('Electrodermal activity, EDA (\muS)', 'FontSize', font_sz_label);
title(sprintf('Decay of arousal (N=%d)', nPart_eda), ...
    'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'FontSize', font_sz, 'Box', 'off', 'TickDir', 'out');

% Annotation box (tau, R^2, fitted amplitude — no R^2 duplication)
text(0.97, 0.95, ...
    {sprintf('\\tau = %.1f s', tau_eda_opt), ...
     sprintf('R^2 = %.3f', R2_eda_best), ...
     sprintf('A_0 = %.2f \\muS', A0_eda_best)}, ...
    'Units', 'normalized', 'FontSize', font_sz_annot, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'Margin', 2);

% Legend — equation form with lambda = e noted (no R^2 here, it's in the box)
legend([h_eda_data, h_eda_fit], ...
    {'\pmSEM', 'R(x) = A_0 exp(-ex) + B'}, ...
    'FontSize', font_sz_annot, 'Location', 'southwest', 'Box', 'off');

text(-0.16, 1.06, 'C', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

%% ========================================================================
%  SAVE
%  ========================================================================

outFile = '../../results/Figure1_Behavioral.png';
exportgraphics(fig, outFile, 'Resolution', 600);
fprintf('  Saved: %s (600 dpi)\n', outFile);

%% ========================================================================
%  SUMMARY
%  ========================================================================

fprintf('\n==========================================================\n');
fprintf('  FIGURE 1 SUMMARY (Bodily Evidence)\n');
fprintf('==========================================================\n');
fprintf('  Panel A — Click Response Times\n');
fprintf('    Clicks: %d,  S(x) R^2 = %.3f,  peak = %.1f s\n', ...
    n_clicks, best.R2, best.peak_time);
fprintf('  Panel B — Haptic Feedback\n');
fprintf('    Asymptote A = %.4f (1/e = %.4f),  R^2 = %.3f\n', ...
    haptic_A, target_1e, haptic_R2);
fprintf('  Panel C — EDA Reliability Decay\n');
fprintf('    Task:  N = %d,  tau = %.1f s,  R^2 = %.3f\n', ...
    nPart_eda, tau_eda_opt, R2_eda_best);
if ~isnan(R2_rest_eda)
    fprintf('    Rest:  R^2 = %.3f\n', R2_rest_eda);
end
fprintf('==========================================================\n');
fprintf('Done.\n');

%% ========================================================================
%  LOCAL FUNCTION
%  ========================================================================

function res = fit_sensitivity_kde(offset, xi_kde, f_kde, T, lambda)
% Fit sensitivity function S(x) = A * x * exp(-lambda * x) to KDE density.
% x = (t - offset) / (T - offset),  lambda fixed at e.
% A estimated by scalar projection.

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
