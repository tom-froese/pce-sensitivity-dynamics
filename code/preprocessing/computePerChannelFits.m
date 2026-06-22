%% computePerChannelFits.m
% =========================================================================
% Per-electrode free-tau fit (sweep tau on grid, pick argmax R^2).
% =========================================================================
%
% Topographic maps of tau*, R^2(tau*), delta-R^2, |trough amplitude|,
% plus comparison to tau-locked (3.90 s) R^2.
%
% REQUIRES UNFILTERED, DC-COUPLED, uV-SCALE RAW EEG (allchannel_data.mat from
% extractAllChannels.m -> data/raw/EEG). Do NOT feed it the 1-40 Hz *-raw.fif
% (aperiodic-exponent input): its 1 Hz high-pass strips the slow S(x) trend.
% A guardrail (section 2b) errors out if filtered/volts-scale data is passed.
%
% Pipeline per channel:
%   - grand-median timecourse across participants
%   - 5-s moving average (valid)
%   - for each tau in TAU_GRID: OLS fit of y = A * s(x) + B  on mask
%     t >= tau, compute R^2
%   - select tau* = argmax R^2, store associated A, B, RMSE
%   - also record tau-locked (3.90 s) R^2 for direct comparison
%
% INPUTS:
%   ../../data/preprocessed/EEG/allchannel_data.mat
%     variables: tTask (1x600), chans (1x64 cell), allTC (62x64x600)
%
% OUTPUTS:
%   ../../results/Figure2_perchannel_fits.csv
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   May 2026
% =========================================================================

%% ========================================================================
%  0. PATHS AND PARAMETERS
%  ========================================================================

scriptDir = fileparts(mfilename('fullpath'));
ROOT      = fullfile(scriptDir, '..', '..');
dataFile  = fullfile(ROOT, 'data', 'preprocessed', 'EEG', 'allchannel_data.mat');
outCSV    = fullfile(ROOT, 'results', 'Figure2_perchannel_fits.csv');

if ~isfile(dataFile)
    error('Missing: %s\nRun extractAllChannels first.', dataFile);
end

E          = exp(1);
T_TRIAL    = 60.0;
TAU_LOCKED = 3.90;
TAU_GRID   = 0.0:0.1:15.0;      % 0 - 15 s, 0.1 s step
SMOOTH_WIN = 5.0;                % seconds
FS         = 10;                 % Hz
MIN_TEFF   = 10.0;              % require at least 10 s of fit window

%% ========================================================================
%  1. LOAD DATA
%  ========================================================================

d      = load(dataFile);
tTask  = d.tTask(:)';           % 1 x 600
chans  = d.chans;               % 1 x 64 cell
allTC  = d.allTC;               % 62 x 64 x 600
n_part = size(allTC, 1);
n_chan = numel(chans);

fprintf('Loaded: n_part = %d, n_chan = %d\n', n_part, n_chan);

%% ========================================================================
%  2. SMOOTHING KERNEL
%  ========================================================================

win   = round(SMOOTH_WIN * FS);
halfK = floor(win / 2);
kernel = ones(1, win) / win;

%% ========================================================================
%  2b. INPUT GUARDRAIL  --  fail loud on filtered / wrong-source data
%  ========================================================================
%  This slow-trend fit REQUIRES unfiltered, DC-coupled, uV-scale RAW EEG
%  (data/raw/EEG/*.mat via extractAllChannels.m). It must NOT be run on the
%  1-40 Hz *-raw.fif files (those are MNE/volts-scale and are exclusively for
%  the aperiodic-exponent analysis): a 1 Hz high-pass destroys the ~0.02-0.03 Hz
%  S(x) slow trend, collapsing every fit to R^2 ~ 0.1. (See dossier 2026-06-22:
%  a FIF-derived allchannel_data.mat silently replaced the raw one and broke the
%  topomap.) The two checks below catch that the moment it happens.

% (i) Scale: MNE/FIF data is in volts (~1e-4); raw .mat is in uV (~1e1-1e3).
maxAbs = max(abs(allTC(:)));
if maxAbs < 1e-2
    error('computePerChannelFits:voltsScale', ...
        ['Input looks VOLTS-scale (max|amp| = %.2e) = MNE/FIF-derived.\n' ...
         'This per-channel slow-trend fit needs uV-scale, DC-coupled RAW EEG via\n' ...
         'extractAllChannels.m, NOT the 1-40 Hz *-raw.fif. Regenerate ' ...
         'allchannel_data.mat from data/raw/EEG. (dossier 2026-06-22)'], maxAbs);
end

% (ii) Slow trend present: the across-channel mean (proto-GSP) must itself fit
%      S(x); a high-pass would flatten it. The true GSP fits at R^2 ~ 0.82, so a
%      0.3 floor is a wide safety margin that only a filtered input would fail.
gmean   = squeeze(mean(median(allTC, 1), 2))';        % 1 x 600 channel-mean of grand median
gmean_s = conv(gmean, kernel, 'valid');
t_gmean = tTask(halfK+1 : halfK+length(gmean_s));
r_gmean = fitAtTau(gmean_s, t_gmean, TAU_LOCKED, T_TRIAL, E, MIN_TEFF);
gmeanR2 = NaN; if ~isempty(r_gmean), gmeanR2 = r_gmean.R2; end
if ~(gmeanR2 >= 0.3)
    error('computePerChannelFits:noSlowTrend', ...
        ['Channel-mean S(x) fit R^2 = %.3f < 0.30 -- the slow trend is missing,\n' ...
         'consistent with high-passed input. This analysis needs UNFILTERED raw\n' ...
         'EEG (the 1-40 Hz FIF strips the slow trend). (dossier 2026-06-22)'], gmeanR2);
end
fprintf('Guardrail OK: uV-scale (max|amp|=%.1f uV), slow trend present (channel-mean R^2=%.2f).\n', ...
    maxAbs, gmeanR2);

%% ========================================================================
%  3. SWEEP ALL CHANNELS
%  ========================================================================

% Pre-allocate result vectors
tau_star_vec   = nan(n_chan, 1);
A_free_vec     = nan(n_chan, 1);
B_free_vec     = nan(n_chan, 1);
R2_free_vec    = nan(n_chan, 1);
RMSE_free_vec  = nan(n_chan, 1);
trough_amp_vec = nan(n_chan, 1);
t_trough_vec   = nan(n_chan, 1);
A_lock_vec     = nan(n_chan, 1);
B_lock_vec     = nan(n_chan, 1);
R2_lock_vec    = nan(n_chan, 1);
RMSE_lock_vec  = nan(n_chan, 1);

for ci = 1:n_chan
    tc_part = squeeze(allTC(:, ci, :));    % 62 x 600

    % Grand median across participants
    med = median(tc_part, 1);

    % 5-s boxcar smooth (valid / same as Python 'valid' mode)
    sm  = conv(med, kernel, 'valid');
    t_sm = tTask(halfK+1 : halfK+length(sm));

    % --- Sweep tau grid ---
    best_R2  = -Inf;
    best_tau = NaN;
    best_A   = NaN;
    best_B   = NaN;
    best_RMSE = NaN;

    for ti = 1:numel(TAU_GRID)
        tau = TAU_GRID(ti);
        r = fitAtTau(sm, t_sm, tau, T_TRIAL, E, MIN_TEFF);
        if isempty(r), continue; end
        if r.R2 > best_R2
            best_R2   = r.R2;
            best_tau  = tau;
            best_A    = r.A;
            best_B    = r.B;
            best_RMSE = r.RMSE;
        end
    end

    % --- Locked-tau fit ---
    r_lock = fitAtTau(sm, t_sm, TAU_LOCKED, T_TRIAL, E, MIN_TEFF);

    % --- Trough amplitude at tau* ---
    trough_x = 1.0 / E;
    t_trough = best_tau + trough_x * (T_TRIAL - best_tau);
    trough_amp = best_A * trough_x * exp(-E * trough_x);

    % --- Store ---
    tau_star_vec(ci)   = best_tau;
    A_free_vec(ci)     = best_A;
    B_free_vec(ci)     = best_B;
    R2_free_vec(ci)    = best_R2;
    RMSE_free_vec(ci)  = best_RMSE;
    trough_amp_vec(ci) = trough_amp;
    t_trough_vec(ci)   = t_trough;
    A_lock_vec(ci)     = r_lock.A;
    B_lock_vec(ci)     = r_lock.B;
    R2_lock_vec(ci)    = r_lock.R2;
    RMSE_lock_vec(ci)  = r_lock.RMSE;
end

dR2_vec = R2_free_vec - R2_lock_vec;

%% ========================================================================
%  4. SAVE CSV
%  ========================================================================

T = table(chans(:), tau_star_vec, A_free_vec, B_free_vec, R2_free_vec, ...
    RMSE_free_vec, trough_amp_vec, t_trough_vec, ...
    A_lock_vec, B_lock_vec, R2_lock_vec, RMSE_lock_vec, dR2_vec, ...
    'VariableNames', {'channel','tau_star','A_free','B_free','R2_free', ...
    'RMSE_free','trough_amp_free','t_trough_free', ...
    'A_lock','B_lock','R2_lock','RMSE_lock','dR2'});

writetable(T, outCSV);
fprintf('\nPer-channel free/locked fits saved: %s\n', outCSV);

%% ========================================================================
%  5. CONSOLE DIAGNOSTICS
%  ========================================================================

fprintf('\nDistribution of tau* across %d channels:\n', n_chan);
fprintf('  median = %.2f s\n', median(tau_star_vec));
fprintf('  mean   = %.2f s\n', mean(tau_star_vec));
fprintf('  IQR    = [%.2f, %.2f] s\n', ...
    quantile(tau_star_vec, 0.25), quantile(tau_star_vec, 0.75));

bins = [0 2; 2 5; 5 8; 8 11; 11 15];
for b = 1:size(bins,1)
    m = tau_star_vec >= bins(b,1) & tau_star_vec < bins(b,2);
    fprintf('  tau* in [%4.1f, %4.1f) s: %2d channels  (median R^2 = %.3f)\n', ...
        bins(b,1), bins(b,2), sum(m), median(R2_free_vec(m)));
end

fprintf('\nTop 10 channels by free-tau R^2:\n');
[~, sortIdx] = sort(R2_free_vec, 'descend');
top10 = sortIdx(1:10);
fprintf('  %-6s  %6s  %6s  %6s  %7s  %8s  %9s  %6s\n', ...
    'Chan', 'tau*', 'R2_f', 'R2_l', 'dR2', 'A_free', 'trough', 't_tr');
for k = 1:10
    i = top10(k);
    fprintf('  %-6s  %6.2f  %6.3f  %6.3f  %+7.3f  %+8.2f  %+9.2f  %6.2f\n', ...
        chans{i}, tau_star_vec(i), R2_free_vec(i), R2_lock_vec(i), ...
        dR2_vec(i), A_free_vec(i), trough_amp_vec(i), t_trough_vec(i));
end

fprintf('\nDone.\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function r = fitAtTau(sm, t_sm, tau, T_trial, E, minTeff)
% fitAtTau  OLS fit of S(x) = A*x*exp(-e*x) + B for t >= tau.
%
%   Returns struct with A, B, R2, RMSE, n  or  [] if Teff too short.

    Teff = T_trial - tau;
    if Teff < minTeff
        r = [];
        return;
    end

    mask  = t_sm >= tau;
    t_fit = t_sm(mask);
    y_fit = sm(mask);

    x = (t_fit - tau) / Teff;
    s = x .* exp(-E * x);
    X = [s(:), ones(numel(s), 1)];

    coef   = X \ y_fit(:);
    A      = coef(1);
    B      = coef(2);
    yhat   = A * s + B;
    ss_res = sum((y_fit(:) - yhat(:)).^2);
    ss_tot = sum((y_fit(:) - mean(y_fit)).^2);

    r.A    = A;
    r.B    = B;
    r.R2   = 1 - ss_res / ss_tot;
    r.RMSE = sqrt(ss_res / numel(y_fit));
    r.n    = numel(y_fit);

    if ss_tot == 0
        r.R2 = NaN;
    end
end

