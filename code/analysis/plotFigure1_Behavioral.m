%% plotFigure1_Behavioral.m
% =========================================================================
% PNAS Figure 1: Bodily Evidence for Reliability Decay Dynamics
% =========================================================================
%
% Two-row, three-column figure combining bodily/behavioral evidence.
% TOP ROW (reduced height): original three panels, with tau locked to the
% joint value tau = 3.9 s for panels A and C.
% BOTTOM ROW (residual band): the residual of each panel's primary fit,
% with a slow-oscillation signature.
%
%   Panel A — Click Response-Time Distribution
%     Click-response times follow P(x) = A * x * exp(-e * x) where
%     x = (t - tau) / (T - tau), lambda = e fixed, tau locked to 3.9 s.
%     Residual: KDE minus the fit P(t), shown raw + lightly smoothed.
%     The residual is left blank during the boot-up (t < tau), since
%     the fit is undefined there.
%
%   Panel B — Haptic Contact Proportion
%     The proportion of participant-trials with active haptic feedback
%     saturates at ~1/e. Residual: smoothed proportion minus exp-rise
%     fit, left blank during the boot-up (t < t_{95%}) for visual
%     consistency with panels A and C.
%
%   Panel C — Electrodermal Activity (EDA) — Reliability Decay
%     R(t) = A0 * exp(-e*(t - tau)/(T - tau)) + B, with lambda = e and
%     tau locked to 3.9 s. Residual: grand-mean EDA minus R(t), shown
%     raw + lightly smoothed, blank during the boot-up (t < tau).
%     Residual traces in all three panels share the same colour so
%     they are visually comparable.
%
% INPUT:
%   ../../data/preprocessed/ClickTimes/ClickResponseTimes.csv
%   ../../data/preprocessed/Haptics/HapticFeedback.csv.gz  (auto-decompressed)
%   ../../data/preprocessed/EDA/EDA_Task_Preprocessed.csv
%   ../../data/preprocessed/EDA/EDA_Rest_Preprocessed.csv
%   ../../data/preprocessed/EEG/globalScalpPotential_stats.mat  (for dissociation)
%
% OUTPUT:
%   ../../results/Figure1_Behavioral.pdf  (vector PDF)
%
% DEPENDENCIES: Signal Processing Toolbox, Optimization Toolbox
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   April 2026
% =========================================================================

%% ========================================================================
%  0. DATA PATHS AND AVAILABILITY CHECK
%  ========================================================================

clickFile    = '../../data/preprocessed/ClickTimes/ClickResponseTimes.csv';
hapticFile   = '../../data/preprocessed/Haptics/HapticFeedback.csv';
hapticGz     = '../../data/preprocessed/Haptics/HapticFeedback.csv.gz';
edaTaskFile  = '../../data/preprocessed/EDA/EDA_Task_Preprocessed.csv';
edaRestFile  = '../../data/preprocessed/EDA/EDA_Rest_Preprocessed.csv';
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

T_trial  = 60;          % Trial duration (s)
lambda   = exp(1);      % Rate parameter fixed at e
TAU_LOCK = 3.9;         % Joint tau locked across panels A and C (and Fig 2)

%% ========================================================================
%  PANEL A — CLICK RESPONSE-TIME DISTRIBUTION WITH SENSITIVITY FIT
%  ========================================================================

fprintf('=== Panel A: Click Response Times (tau locked = %.2f s) ===\n', TAU_LOCK);

data_click = readtable(clickFile);
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

% Fit with tau LOCKED to 3.9 s
best = fit_sensitivity_kde(TAU_LOCK, xi_kde, f_kde, T_trial, lambda);

fprintf('  Onset lag (locked): %.1f s,  R^2 = %.4f,  peak = %.1f s,  A = %.3f\n', ...
    best.offset, best.R2, best.peak_time, best.A);

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
    sig   = data_haptic.HapticFeedback(idx);
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

% Boot-up (t_95%) of the exp-rise fit — precompute here so residual
% masking below can reference it.
haptic_onset = pFit(3);
haptic_tau_h = pFit(2);
haptic_t95   = haptic_onset - haptic_tau_h * log(1 - 0.95);

fprintf('  A = %.4f (1/e = %.4f),  R^2 = %.4f\n', haptic_A, target_1e, haptic_R2);
fprintf('  Boot-up: onset = %.2f s, tau = %.2f s, t_95%% = %.1f s\n', ...
    haptic_onset, haptic_tau_h, haptic_t95);

%% ========================================================================
%  PANEL C — EDA RELIABILITY DECAY (tau locked)
%  ========================================================================

fprintf('\n=== Panel C: EDA Reliability Decay (tau locked = %.2f s) ===\n', TAU_LOCK);

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
part_matrix    = part_matrix(~any(isnan(part_matrix), 2), :);
nPart_eda      = size(part_matrix, 1);
time_eda       = (0:task_samples-1) / fs_eda;
grand_mean_eda = mean(part_matrix, 1);
grand_sem_eda  = std(part_matrix, 0, 1) / sqrt(nPart_eda);

model_eda = @(params, x) params(1) * exp(-exp(1) * x) + params(2);
opts_eda  = optimoptions('lsqcurvefit', 'Display', 'off', ...
    'MaxIterations', 10000, 'MaxFunctionEvaluations', 30000, ...
    'TolFun', 1e-12, 'TolX', 1e-12);

% Locked-tau fit only
tau_eda_opt = TAU_LOCK;
T_eff_eda   = T_trial - tau_eda_opt;
tau_samp_eda = round(tau_eda_opt * fs_eda);
idx_eda_opt  = (tau_samp_eda + 1):task_samples;
x_eda_opt    = (time_eda(idx_eda_opt) - tau_eda_opt) / T_eff_eda;
y_eda_opt    = grand_mean_eda(idx_eda_opt);

A0_init = y_eda_opt(1) - y_eda_opt(end); B_init = y_eda_opt(end);
params  = lsqcurvefit(model_eda, [A0_init, B_init], x_eda_opt, y_eda_opt, ...
    [0 0], [50 50], opts_eda);
A0_eda_best = params(1);
B_eda_best  = params(2);
yfit_eda_opt = model_eda([A0_eda_best, B_eda_best], x_eda_opt);

SS_res = sum((y_eda_opt - yfit_eda_opt).^2);
SS_tot = sum((y_eda_opt - mean(y_eda_opt)).^2);
R2_eda_best = 1 - SS_res / SS_tot;

fprintf('  N = %d,  tau (locked) = %.1f s,  R^2 = %.4f,  A0 = %.3f,  B = %.3f\n', ...
    nPart_eda, tau_eda_opt, R2_eda_best, A0_eda_best, B_eda_best);

%% ---- REST EDA (for dissociation, unchanged) ----

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

if isfile(gspStatsFile)
    gsp_stats = load(gspStatsFile, 'grandR2', 'grandRestR2');
    fprintf('\n  === DISSOCIATION SUMMARY ===\n');
    fprintf('  EDA  R(x) fit:   Task R^2 = %.3f,  Rest R^2 = %.3f\n', ...
        R2_eda_best, R2_rest_eda);
    fprintf('  GSP  S(x) fit:   Task R^2 = %.3f,  Rest R^2 = %.3f\n', ...
        gsp_stats.grandR2, gsp_stats.grandRestR2);
end

%% ========================================================================
%  RESIDUALS — compute for all three panels
%  ========================================================================

fprintf('\n=== Residuals (bottom row) ===\n');

% ---- Panel A residual: KDE (scaled to counts) minus P(t) -----------------
% Evaluate P(t) on the KDE grid with locked tau
t_mod_full = xi_kde;
x_mod_full = (t_mod_full - best.offset) / best.T_eff;
x_mod_full = max(x_mod_full, 1e-12);
P_fit_full = best.A .* x_mod_full .* exp(-lambda .* x_mod_full);
P_fit_full(t_mod_full < best.offset) = 0;
% Scale everything to counts for comparability with top panel
kde_scale = n_clicks * bin_width;
resid_A_full = (f_kde - P_fit_full) * kde_scale;
% Raw, per-bin residual = actual counts minus the fit evaluated at bin centres
P_fit_bins   = interp1(xi_kde, P_fit_full, bin_centres, 'linear', 0);
resid_A_bins = counts - P_fit_bins * kde_scale;
% Blank out the boot-up window (t < tau) — the fit is not defined there,
% so showing a disjointed residual for that stretch is misleading.
resid_A_full(xi_kde       < TAU_LOCK) = NaN;
resid_A_bins(bin_centres  < TAU_LOCK) = NaN;
fprintf('  Panel A residual (t >= %.1f s):  SD = %.2f counts (KDE), %.2f counts (bins)\n', ...
    TAU_LOCK, std(resid_A_full, 'omitnan'), std(resid_A_bins, 'omitnan'));

% ---- Panel B residual: smoothed proportion minus exp-rise fit ------------
propFit_h      = expRise(pFit, timeVec_h);
resid_B        = (propSmooth - propFit_h) * 100;   % in percent
resid_B_raw    = (propON     - propFit_h) * 100;   % raw per-sample residual
resid_B_smooth = smoothdata(resid_B, 'gaussian', round(2 * fs_h));  % ~2 s window
% Blank out the boot-up (t < t_{95%}) for visual consistency with A and C
resid_B_raw(timeVec_h    < haptic_t95) = NaN;
resid_B_smooth(timeVec_h < haptic_t95) = NaN;
fprintf('  Panel B residual (t >= %.1f s):  SD = %.2f%% (raw), %.2f%% (smoothed)\n', ...
    haptic_t95, std(resid_B_raw, 'omitnan'), std(resid_B_smooth, 'omitnan'));

% ---- Panel C residual: grand-mean EDA minus R(t) fit ---------------------
% Evaluate fit on ALL samples, using x = (t - tau)/T_eff  clipped at 0 pre-tau
x_full_eda     = max((time_eda - tau_eda_opt) / T_eff_eda, 0);
yfit_eda_full  = model_eda([A0_eda_best, B_eda_best], x_full_eda);
resid_C_full   = grand_mean_eda - yfit_eda_full;
resid_C_smooth = smoothdata(resid_C_full, 'gaussian', round(2 * fs_eda));  % ~2 s window
% Blank out the boot-up window (t < tau)
resid_C_full(time_eda   < tau_eda_opt) = NaN;
resid_C_smooth(time_eda < tau_eda_opt) = NaN;
fprintf('  Panel C residual (t >= %.1f s):  SD = %.3f uS (raw), %.3f uS (smoothed)\n', ...
    tau_eda_opt, std(resid_C_full, 'omitnan'), std(resid_C_smooth, 'omitnan'));

%% ========================================================================
%  CREATE FIGURE — 2 rows x 3 panels
%  ========================================================================

fprintf('\nCreating Figure 1 ...\n');

% --- Colour palette ---
col_data  = [0.20 0.40 0.73];
col_kde   = [0.12 0.30 0.55];
col_fit   = [0.80 0.15 0.15];
col_grey  = [0.50 0.50 0.50];
col_eda   = [0.85 0.32 0.10];
col_black = [0.10 0.10 0.10];

font_sz       = 9;
font_sz_label = 10;
font_sz_title = 11;
font_sz_annot = 8;
font_sz_panel = 14;

% PNAS double-column width: 7.09 in
fig_w = 7.5;  fig_h = 4.2;
fig = figure('Units', 'inches', 'Position', [0.5 0.5 fig_w fig_h], ...
    'Color', 'w', 'PaperUnits', 'inches', ...
    'PaperSize', [fig_w fig_h], 'PaperPosition', [0 0 fig_w fig_h]);

% Margins (normalised) — top margin generous so that bold panel labels
% "A", "B", "C" (rendered at y = 1.08 of axes) and the panel titles are
% never clipped by the figure boundary.
ml = 0.07;  mr = 0.02;  mb = 0.08;  mt = 0.11;
gap_ab = 0.07;  gap_bc = 0.07;
gap_rows = 0.10;
total_w = 1 - ml - mr - gap_ab - gap_bc;
% Equal panel widths so the time axes of A, B, C line up one-to-one and
% the residuals below them are directly comparable.
pw_a = total_w / 3;
pw_b = total_w / 3;
pw_c = total_w / 3;

% Row heights
avail_h = 1 - mb - mt - gap_rows;
ph_top = avail_h * 0.70;       % top row is the larger band (2/3 of previous)
ph_bot = avail_h - ph_top;     % residual band gets ~30% of vertical space
row_bot_y = mb;
row_top_y = mb + ph_bot + gap_rows;

% Column x-positions
x_a = ml;
x_b = ml + pw_a + gap_ab;
x_c = ml + pw_a + gap_ab + pw_b + gap_bc;

%% ---- TOP-ROW Panel A: Click-time distribution + sensitivity fit --------

ax_a = axes('Position', [x_a, row_top_y, pw_a, ph_top]);
hold on;

bar(bin_centres, counts, 1, ...
    'FaceColor', col_data, 'FaceAlpha', 0.30, ...
    'EdgeColor', 'w', 'LineWidth', 0.5);

plot(xi_kde, f_kde * kde_scale, '-', 'Color', col_kde, 'LineWidth', 1.5);

t_mod = linspace(best.offset, T_trial, 500);
x_mod = (t_mod - best.offset) / best.T_eff;
x_mod = max(x_mod, 1e-12);
y_mod = (best.A .* x_mod .* exp(-lambda .* x_mod)) * kde_scale;
plot(t_mod, y_mod, '-', 'Color', col_fit, 'LineWidth', 2);

xline(best.peak_time, ':', 'Color', col_fit, 'LineWidth', 1, 'Alpha', 0.7);

xline(best.offset, ':', 'Color', col_grey, 'LineWidth', 1.2);

xlim([0 T_trial]);
ylim([0 max(counts) * 1.35]);
ylabel('Number of clicks', 'FontSize', font_sz_label);
title('Click response times', 'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'FontSize', font_sz, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', []);

% Place the T_eff/e label well below the peak of the red fit curve so
% the curve sits cleanly above the annotation. A white background keeps
% it readable over the histogram bars.
peak_y = best.A * (1/lambda) * exp(-1) * kde_scale;
text(best.peak_time + 0.4, peak_y * 0.55, 'T_{eff}/e', ...
    'FontSize', font_sz_annot, 'Color', col_fit, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'BackgroundColor', 'w', 'Margin', 1);

lg_a = legend({'Clicks', 'KDE', 'P(x) = A \cdot x \cdot exp(-ex)'}, ...
    'FontSize', font_sz_annot, 'Location', 'northwest');
lg_a.BoxFace.ColorType = 'truecoloralpha';
lg_a.BoxFace.ColorData = uint8([255; 255; 255; 255]);
lg_a.EdgeColor = [0.7 0.7 0.7];

text(0.97, 0.05, ...
    {sprintf('\\tau = %.1f s (locked)', best.offset), ...
     sprintf('R^2 = %.3f', best.R2), ...
     sprintf('peak = %.1f s', best.peak_time)}, ...
    'Units', 'normalized', 'FontSize', font_sz_annot, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
    'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'Margin', 2);

text(-0.16, 1.08, 'A', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');
hold off;

%% ---- TOP-ROW Panel B: Haptic contact proportion ------------------------

ax_b = axes('Position', [x_b, row_top_y, pw_b, ph_top]);
hold on;

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

yline(haptic_A * 100, '--', 'Color', col_fit, 'LineWidth', 1, 'Alpha', 0.6, ...
    'HandleVisibility', 'off');

yline(target_1e * 100, ':', 'Color', col_fit, 'LineWidth', 1.2, 'Alpha', 0.6, ...
    'HandleVisibility', 'off');
text(57, target_1e * 100 + 1.2, '1/e', 'FontSize', font_sz_annot, ...
    'Color', col_fit, 'HorizontalAlignment', 'right');

xline(haptic_t95, ':', 'Color', col_grey, 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');

hold off;
xlim([0 60]);
ylim([0 max((propON + semON) * 100) * 1.1]);
ylabel('Trials with active contact (%)', 'FontSize', font_sz_label);
title('Haptic feedback activation', 'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'FontSize', font_sz, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', []);

lg_b = legend('Location', 'southeast', 'Box', 'off', 'FontSize', font_sz_annot);
lg_b.ItemTokenSize = [15, 8];

text(0.12, 0.95, sprintf('t_{95%%} = %.1f s', haptic_t95), ...
    'Units', 'normalized', 'FontSize', font_sz_annot, ...
    'VerticalAlignment', 'top', 'Color', col_grey);

text(-0.16, 1.08, 'B', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

%% ---- TOP-ROW Panel C: EDA reliability decay ----------------------------

ax_c = axes('Position', [x_c, row_top_y, pw_c, ph_top]);
hold on;

fill([time_eda, fliplr(time_eda)], ...
    [grand_mean_eda + grand_sem_eda, fliplr(grand_mean_eda - grand_sem_eda)], ...
    col_eda, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

h_eda_data = plot(time_eda, grand_mean_eda, '-', 'Color', col_eda, 'LineWidth', 1.5);

h_eda_fit = plot(time_eda(idx_eda_opt), yfit_eda_opt, '--', ...
    'Color', col_black, 'LineWidth', 2);

xline(tau_eda_opt, ':', 'Color', col_grey, 'LineWidth', 1.2);

hold off;

xlim([0 60]);
ylabel('EDA (\muS)', 'FontSize', font_sz_label);
title(sprintf('Decay of arousal (N=%d)', nPart_eda), ...
    'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'FontSize', font_sz, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', []);

text(0.97, 0.95, ...
    {sprintf('\\tau = %.1f s (locked)', tau_eda_opt), ...
     sprintf('R^2 = %.3f', R2_eda_best), ...
     sprintf('A_0 = %.2f \\muS', A0_eda_best)}, ...
    'Units', 'normalized', 'FontSize', font_sz_annot, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'Margin', 2);

legend([h_eda_data, h_eda_fit], ...
    {'\pmSEM', 'R(x) = A_0 exp(-ex) + B'}, ...
    'FontSize', font_sz_annot, 'Location', 'southwest', 'Box', 'off');

text(-0.16, 1.08, 'C', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

%% ---- BOTTOM-ROW Panel A residual ---------------------------------------

ax_ar = axes('Position', [x_a, row_bot_y, pw_a, ph_bot]);
hold on;
yline(0, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
% Raw per-bin residual (counts minus fit at each 2-s bin) drawn as bars
% to mirror the histogram in the top panel; smoothed KDE-based residual
% overlaid in bold (same role as the KDE line up top).
bar(bin_centres, resid_A_bins, 1, ...
    'FaceColor', col_data, 'FaceAlpha', 0.30, ...
    'EdgeColor', 'w', 'LineWidth', 0.5);
plot(xi_kde, resid_A_full, '-', 'Color', col_data, 'LineWidth', 1.5);
% Boot-up marker (same dotted grey line as the top panel)
xline(TAU_LOCK, ':', 'Color', col_grey, 'LineWidth', 1.2);
hold off;
xlim([0 T_trial]);
yl = ylim; ylim([-max(abs(yl)) max(abs(yl))]);
xlabel('Time (s)', 'FontSize', font_sz_label);
ylabel('Residual (counts)', 'FontSize', font_sz_label);
set(gca, 'FontSize', font_sz, 'Box', 'off', 'TickDir', 'out');

%% ---- BOTTOM-ROW Panel B residual ---------------------------------------
% No boot-up marker: the haptic channel has no tau-offset boot-up.

ax_br = axes('Position', [x_b, row_bot_y, pw_b, ph_bot]);
hold on;
yline(0, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
% The haptic residual is sampled at 100 Hz, so the raw trace is ~200x
% denser than panel A's per-bin residual. Drop alpha and linewidth so
% the scatter reads as a faint background and the smoothed line stands
% out clearly.
plot(timeVec_h, resid_B_raw,    '-', 'Color', [col_data 0.10], 'LineWidth', 0.3);
plot(timeVec_h, resid_B_smooth, '-', 'Color', col_data,        'LineWidth', 1.5);
hold off;
xlim([0 60]);
yl = ylim; ylim([-max(abs(yl)) max(abs(yl))]);
xlabel('Time (s)', 'FontSize', font_sz_label);
ylabel('Residual (%)', 'FontSize', font_sz_label);
set(gca, 'FontSize', font_sz, 'Box', 'off', 'TickDir', 'out');

%% ---- BOTTOM-ROW Panel C residual ---------------------------------------

ax_cr = axes('Position', [x_c, row_bot_y, pw_c, ph_bot]);
hold on;
yline(0, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
% Use the same blue as panels A and B — residual traces are conceptually
% comparable across panels, so a shared colour makes that explicit.
plot(time_eda, resid_C_full,   '-', 'Color', [col_data 0.35], 'LineWidth', 0.6);
plot(time_eda, resid_C_smooth, '-', 'Color', col_data,       'LineWidth', 1.5);
% Boot-up marker
xline(tau_eda_opt, ':', 'Color', col_grey, 'LineWidth', 1.2);
hold off;
xlim([0 60]);
yl = ylim; ylim([-max(abs(yl)) max(abs(yl))]);
xlabel('Time (s)', 'FontSize', font_sz_label);
ylabel('Residual (\muS)', 'FontSize', font_sz_label);
set(gca, 'FontSize', font_sz, 'Box', 'off', 'TickDir', 'out');

%% ========================================================================
%  SAVE
%  ========================================================================

outFile = '../../results/Figure1_Behavioral.pdf';
exportgraphics(fig, outFile, 'ContentType', 'vector');
fprintf('  Saved: %s (vector PDF)\n', outFile);

%% ========================================================================
%  SUMMARY
%  ========================================================================

fprintf('\n==========================================================\n');
fprintf('  FIGURE 1 SUMMARY (Bodily Evidence)\n');
fprintf('==========================================================\n');
fprintf('  Panel A — Click Response Times  (tau locked = %.2f s)\n', TAU_LOCK);
fprintf('    Clicks: %d,  P(x) R^2 = %.3f,  peak = %.1f s,  A = %.3f\n', ...
    n_clicks, best.R2, best.peak_time, best.A);
fprintf('    Residual SD: %.2f counts (KDE), %.2f counts (bins)\n', ...
    std(resid_A_full, 'omitnan'), std(resid_A_bins, 'omitnan'));
fprintf('  Panel B — Haptic Feedback\n');
fprintf('    Asymptote A = %.4f (1/e = %.4f),  R^2 = %.3f\n', ...
    haptic_A, target_1e, haptic_R2);
fprintf('    Residual SD: %.2f%% (raw), %.2f%% (smoothed)\n', ...
    std(resid_B_raw, 'omitnan'), std(resid_B_smooth, 'omitnan'));
fprintf('  Panel C — EDA Reliability Decay  (tau locked = %.2f s)\n', TAU_LOCK);
fprintf('    Task:  N = %d,  R^2 = %.3f\n', nPart_eda, R2_eda_best);
fprintf('    Residual SD: %.3f uS (raw), %.3f uS (smoothed)\n', ...
    std(resid_C_full, 'omitnan'), std(resid_C_smooth, 'omitnan'));
if ~isnan(R2_rest_eda)
    fprintf('    Rest:  R^2 = %.3f\n', R2_rest_eda);
end
fprintf('==========================================================\n');
fprintf('Done.\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function res = fit_sensitivity_kde(offset, xi_kde, f_kde, T, lambda)
% Fit click-probability curve P(x) = A * x * exp(-lambda * x) to KDE density.
% x = (t - offset) / (T - offset), lambda fixed. A by scalar projection.

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

