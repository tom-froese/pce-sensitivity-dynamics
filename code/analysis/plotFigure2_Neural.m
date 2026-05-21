%% plotFigure2_Neural.m
% =========================================================================
% PNAS Figure 2 — Neural Evidence
% =========================================================================
%
% Two-row figure:
%
%   TOP ROW (A) — Scalp-map triptych with progressive filtering:
%              (i)   R^2_free (all 64 channels) — where does the model fit?
%              (ii)  tau* (filtered: tau* in [2, 6] s) — timing gradient
%              (iii) Signed trough magnitude (filtered: tau* in [2, 6] s
%                    AND R^2_lock >= 0.70) — effect size
%
%   BOTTOM ROW (B) — ROI-mean scalp potential of the left and right
%              parietal hemispheres (L: P1/P3/P5/P7; R: P2/P4/P6/P8),
%              overlaid with locked-tau S(x) = A*x*exp(-e*x) + B fits
%              (tau = 3.90 s, fixed from the full-parietal free-tau
%              optimum).
%
% INPUTS:
%   ../../data/preprocessed/EEG/parietal_hemisphere_data.mat
%       variables: tTask, left_TC, mid_TC, right_TC  (62x600 each)
%   ../../data/preprocessed/EEG/globalScalpPotential_data.mat
%       variables: roi.names, roiParticipantTC  (cell of 62x600 per ROI)
%   ../../results/FigureS2_GSP_TopoMap_FreeTau_perchannel.csv
%
% OUTPUTS:
%   ../../results/Figure2_Neural.pdf  (vector PDF)
%   ../../results/Figure2_Neural.png  (300 dpi raster for inspection)
%
% DEPENDENCIES: EEGLAB (topoplot function)
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   May 2026
% =========================================================================


 
%% ========================================================================
%  0. PATHS AND PARAMETERS
%  ========================================================================

scriptDir = fileparts(mfilename('fullpath'));
ROOT      = fullfile(scriptDir, '..', '..');
hemFile   = fullfile(ROOT, 'data', 'preprocessed', 'EEG', ...
                     'parietal_hemisphere_data.mat');
gspFile   = fullfile(ROOT, 'data', 'preprocessed', 'EEG', ...
                     'globalScalpPotential_data.mat');
csvFile   = fullfile(ROOT, 'results', ...
                     'FigureS2_GSP_TopoMap_FreeTau_perchannel.csv');
outPDF    = fullfile(ROOT, 'results', 'Figure2_Neural.pdf');
outPNG    = fullfile(ROOT, 'results', 'Figure2_Neural.png');

for f = {hemFile, gspFile, csvFile}
    if ~isfile(f{1})
        error('Missing: %s', f{1});
    end
end

% EEGLAB for topographic maps
if ~exist('topoplot', 'file')
    eeglabDir = fullfile(getenv('HOME'), 'Documents', 'MATLAB', 'eeglab2026.0.0');
    if isfolder(eeglabDir), addpath(eeglabDir); eeglab('nogui');
    else, error('EEGLAB not found. Install to ~/Documents/MATLAB/eeglab2026.0.0/');
    end
end
locFile      = fullfile(fileparts(which('eeglab')), 'plugins', 'dipfit', ...
    'standard_BEM', 'elec', 'standard_1005.elc');
allLocs      = readlocs(locFile);
% standard_1005.elc uses ASA convention: +X = right ear, +Y = nose, +Z = vertex.
% EEGLAB internally expects:              +X = nose,     +Y = left ear, +Z = vertex.
% Mapping:  EEGLAB_X = file_Y,  EEGLAB_Y = -file_X,  EEGLAB_Z = file_Z.
for iLoc = 1:length(allLocs)
    origX = allLocs(iLoc).X;
    allLocs(iLoc).X =  allLocs(iLoc).Y;    % nose
    allLocs(iLoc).Y = -origX;               % left ear (negate right ear)
end
allLocs      = convertlocs(allLocs, 'cart2all');   % populate theta/radius for topoplot
allLabelsStd = {allLocs.labels};

E            = exp(1);
T_TRIAL      = 60.0;
TAU_LOCKED   = 3.90;
SMOOTH_WIN_S = 5.0;
FS           = 10;

% Strict-filter criteria for the top-row topomaps
TAU_LO = 2.0;
TAU_HI = 6.0;
R2_MIN = 0.70;

%% ========================================================================
%  1. LOAD HEMISPHERE DATA
%  ========================================================================

d_hem  = load(hemFile);
tTask  = d_hem.tTask(:)';       % 1 x 600
left   = d_hem.left_TC;         % 62 x 600
right  = d_hem.right_TC;        % 62 x 600
n_part = size(left, 1);

fprintf('Loaded hemisphere data: n_part = %d, samples = %d\n', ...
    n_part, size(left, 2));

%% ========================================================================
%  2. LOAD GSP DATA (full parietal ROI)
%  ========================================================================

d_gsp      = load(gspFile);
roi_names  = d_gsp.roi.names;
par_idx    = find(strcmpi(roi_names, 'Parietal'), 1);
par_full   = d_gsp.roiParticipantTC{par_idx};   % 62 x 600

fprintf('Loaded GSP data: Parietal ROI = index %d\n', par_idx);

%% ========================================================================
%  3. SMOOTHED GRAND MEDIANS + LOCKED-TAU FITS
%  ========================================================================

[t_sm, sm_left]  = smoothMedian(left,     tTask, SMOOTH_WIN_S, FS);
[~,    sm_right] = smoothMedian(right,    tTask, SMOOTH_WIN_S, FS);
[~,    sm_full]  = smoothMedian(par_full, tTask, SMOOTH_WIN_S, FS);

fit_full  = fitLocked(t_sm, sm_full,  TAU_LOCKED, T_TRIAL, E);
fit_left  = fitLocked(t_sm, sm_left,  TAU_LOCKED, T_TRIAL, E);
fit_right = fitLocked(t_sm, sm_right, TAU_LOCKED, T_TRIAL, E);

for pair = {'Full parietal', fit_full;
            'Left  parietal', fit_left;
            'Right parietal', fit_right}'
    name = pair{1}; fit = pair{2};
    fprintf('%s (locked tau = %.2f s):  A = %+.2f, B = %+.2f, R^2 = %.3f, trough at t = %.2f s\n', ...
        name, TAU_LOCKED, fit.A, fit.B, fit.R2, fit.t_trough);
end

%% ========================================================================
%  4. SCALP-MAP DATA (progressive filtering across three panels)
%  ========================================================================
%
%  Panel 1: R^2_free for ALL channels (no filter) — where does the model fit?
%  Panel 2: tau* for channels with tau* in [2, 6] s — timing gradient
%  Panel 3: Trough magnitude for channels passing both tau*-range AND
%           R^2_lock >= 0.70 — effect size in the best-fit region
%
%  Each panel adds one filter; annotations explain the progressive narrowing.

df      = readtable(csvFile);
N_TOTAL = height(df);

% --- Matching helper: build chanlocs + aligned data for a channel subset ---
buildChanlocs = @(names) matchToMontage(names, allLocs, allLabelsStd);

% Panel 1: ALL channels
[locs_all, mask_all] = buildChanlocs(df.channel);
R2_free_all          = df.R2_free(mask_all);
N_ALL                = sum(mask_all);

fprintf('\nPanel 1 (R^2 free): %d/%d channels matched to montage\n', N_ALL, N_TOTAL);

% Panel 2: filter by tau* range only
mask_tau          = (df.tau_star >= TAU_LO) & (df.tau_star <= TAU_HI);
df_tau            = df(mask_tau, :);
N_DROP_TAU        = sum(~mask_tau);
N_AFTER_TAU       = sum(mask_tau);
[locs_tau, mask_tau_m] = buildChanlocs(df_tau.channel);
tau_vals          = df_tau.tau_star(mask_tau_m);
N_TAU             = sum(mask_tau_m);

fprintf('Panel 2 (tau*):     %d channels after tau*-range filter (dropped %d)\n', N_TAU, N_DROP_TAU);

% Panel 3: filter by tau* range AND R^2_lock
mask_strict       = mask_tau & (df.R2_lock >= R2_MIN);
df_strict         = df(mask_strict, :);
N_DROP_R2         = sum(mask_tau & ~(df.R2_lock >= R2_MIN));
N_STRICT          = sum(mask_strict);
[locs_strict, mask_str_m] = buildChanlocs(df_strict.channel);
A_lock_vals       = df_strict.A_lock(mask_str_m);
N_KEPT            = sum(mask_str_m);

trough_shape      = (1.0 / E) * exp(-1.0);
trough_signed     = A_lock_vals * trough_shape;

fprintf('Panel 3 (trough):   %d channels after R^2_lock >= %.2f filter (dropped %d more)\n', ...
    N_KEPT, R2_MIN, N_DROP_R2);

%% ========================================================================
%  5. FIGURE
%  ========================================================================

% --- Styling (match existing figures) ---
font_sz       = 9;
font_sz_label = 10;
font_sz_title = 11;
font_sz_annot = 8;
font_sz_panel = 14;

fig_w = 9.0;  fig_h = 9.5;
fig = figure('Units', 'inches', 'Position', [0.5 0.5 fig_w fig_h], ...
    'Color', 'w', 'PaperUnits', 'inches', ...
    'PaperSize', [fig_w fig_h], 'PaperPosition', [0 0 fig_w fig_h]);

set(fig, 'DefaultAxesFontSize', font_sz, ...
         'DefaultAxesTickDir', 'out', ...
         'DefaultAxesBox', 'off', ...
         'DefaultAxesColor', 'w');

% --- Layout ---
% Top row: 3 topomap panels
% Bottom row: single wide panel
ml = 0.07;  mr = 0.04;  mb = 0.08;  mt = 0.08;
gap_h = 0.02;          % horizontal gap between top panels
gap_v = 0.10;           % vertical gap between rows

avail_w = 1 - ml - mr - 2*gap_h;
pw_top  = avail_w / 3;

avail_h = 1 - mb - mt - gap_v;
ph_top  = avail_h * 0.58;
ph_bot  = avail_h - ph_top;

row_bot_y = mb;
row_top_y = mb + ph_bot + gap_v;

x1 = ml;
x2 = ml + pw_top + gap_h;
x3 = ml + 2*(pw_top + gap_h);

%% ---- TOP ROW: three scalp maps (EEGLAB topoplot) ----

% Strategy: EEGLAB topoplot uses the *figure-level* colormap, so per-axes
% colormaps have no effect.  freezeTopoColors (converting indexed CData to
% true-colour RGB) fails on some of topoplot's patch types.
%
% Robust workaround: render each topoplot in its own temporary figure with
% the correct colormap, capture as an RGB image via print(), then place
% the rasterised image into the main figure.  Each temp figure is an
% independent colormap scope, so there is no cross-panel contamination.

TOPO_DPI = 300;

topoSpecs = { ...
    R2_free_all, locs_all,    viridis_cmap(256),  [0 1],               ...
        'R^2_{free}  (best-fit quality)', 'R^2', ...
        sprintf('all %d channels; free \\tau^{*} per channel', N_ALL);
    tau_vals,    locs_tau,     coolwarm_cmap(256),  [TAU_LO TAU_HI],    ...
        '\tau^{*}  (s)', '\tau^{*} (s)', ...
        sprintf('filter: \\tau^{*} \\in [%.0f, %.0f] s  (%d/%d ch.)', TAU_LO, TAU_HI, N_TAU, N_TOTAL);
    trough_signed, locs_strict, magma_cmap(256),   [-max(abs(trough_signed)) 0], ...
        'Signed trough magnitude  A \cdot s(1/e)', '\muV (trough, < 0)', ...
        sprintf('+ R^2_{lock} \\geq %.2f  (%d/%d ch.)', R2_MIN, N_KEPT, N_TOTAL) ...
};

ax_topos = gobjects(3,1);
xPositions = [x1, x2, x3];

for ip = 1:3
    vals   = topoSpecs{ip, 1};
    locs   = topoSpecs{ip, 2};
    cmap   = topoSpecs{ip, 3};
    mlims  = topoSpecs{ip, 4};
    ttl    = topoSpecs{ip, 5};
    cblbl  = topoSpecs{ip, 6};
    annot  = topoSpecs{ip, 7};

    % --- Render in a temporary figure ---
    tmpFig = figure('Units', 'pixels', 'Position', [100 100 500 550], ...
        'Color', 'w', 'Visible', 'off');
    colormap(tmpFig, cmap);
    tmpAx = axes('Parent', tmpFig, 'Position', [0.05 0.15 0.90 0.78]);
    topoplot(vals, locs, 'maplimits', mlims, ...
        'numcontour', 6, 'electrodes', 'ptslabels', 'efontsize', 8, ...
        'style', 'fill', 'conv', 'off');
    title(ttl, 'FontSize', font_sz_title);
    cb = colorbar('southoutside');
    cb.Label.String = cblbl; cb.FontSize = 8;

    % Capture as RGB image
    drawnow;
    frame = print(tmpFig, '-RGBImage', sprintf('-r%d', TOPO_DPI));
    close(tmpFig);

    % --- Place captured image into main figure ---
    ax_topos(ip) = axes('Position', [xPositions(ip), row_top_y, pw_top, ph_top], ...
        'Parent', fig);
    image(frame, 'Parent', ax_topos(ip));
    axis(ax_topos(ip), 'image', 'off');

    % Annotation below the image
    text(0.5, -0.02, annot, ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', font_sz_annot, 'Parent', ax_topos(ip));
end

%% ---- BOTTOM ROW: parietal L/R time courses + S(x) fits ----

ax_bot = axes('Position', [ml, row_bot_y, 1-ml-mr, ph_bot]);
hold(ax_bot, 'on');

col_L = [0.12 0.47 0.71];    % #1f77b4
col_R = [0.84 0.15 0.16];    % #d62728

% Baseline-subtract: remove fit offset B so the trough dips negative,
% consistent with the signed trough magnitudes in the top-right panel.
sm_left_bc  = sm_left  - fit_left.B;
sm_right_bc = sm_right - fit_right.B;
yhat_left_bc  = fit_left.yhat  - fit_left.B;
yhat_right_bc = fit_right.yhat - fit_right.B;

% Raw grand-median time courses (de-emphasised)
plot(ax_bot, t_sm, sm_left_bc,  '-', 'Color', [col_L 0.35], 'LineWidth', 1.0, ...
    'DisplayName', 'Left parietal  (P1/P3/P5/P7)');
plot(ax_bot, t_sm, sm_right_bc, '-', 'Color', [col_R 0.35], 'LineWidth', 1.0, ...
    'DisplayName', 'Right parietal (P2/P4/P6/P8)');

% Locked-tau S(x) fits (only for t >= TAU_LOCKED)
plot(ax_bot, fit_left.t_fit, yhat_left_bc, '-', ...
    'Color', col_L, 'LineWidth', 2.2, ...
    'DisplayName', sprintf('L fit: A=%+.1f \\muV,  R^2=%.2f', ...
    fit_left.A, fit_left.R2));
plot(ax_bot, fit_right.t_fit, yhat_right_bc, '-', ...
    'Color', col_R, 'LineWidth', 2.2, ...
    'DisplayName', sprintf('R fit: A=%+.1f \\muV,  R^2=%.2f', ...
    fit_right.A, fit_right.R2));

% Per-hemisphere trough markers (B already subtracted -> trough = A*s(1/e))
for pair = {'L', fit_left, col_L; 'R', fit_right, col_R}'
    fit = pair{2};
    col = pair{3};
    y_trough = fit.A * (1.0/E) * exp(-1.0);
    plot(ax_bot, fit.t_trough, y_trough, 'o', ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', ...
        'MarkerSize', 7, 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
end

% Locked-tau landmark + zero line
xline(ax_bot, TAU_LOCKED, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, ...
    'HandleVisibility', 'off');
yline(ax_bot, 0, '-', 'Color', [0 0 0 0.3], 'LineWidth', 0.4, ...
    'HandleVisibility', 'off');

xlim(ax_bot, [t_sm(1) t_sm(end)]);
xlabel(ax_bot, 'Time (s)', 'FontSize', font_sz_label);
ylabel(ax_bot, 'ROI-mean scalp potential, baseline-corrected (\muV)', 'FontSize', font_sz_label);
title(ax_bot, sprintf('Parietal hemispheres -- grand median with S(x) = A x e^{-ex} fits  (\\tau = %.2f s)', TAU_LOCKED), ...
    'FontSize', font_sz_title);

lg = legend(ax_bot, 'Location', 'north', 'NumColumns', 2, ...
    'FontSize', font_sz_annot, 'Box', 'off');
lg.Position(2) = lg.Position(2) + 0.01;

set(ax_bot, 'FontSize', font_sz, 'Box', 'off', 'TickDir', 'out');
hold(ax_bot, 'off');

%% ---- Panel letters and suptitle ----

text(ax_topos(1), -0.08, 1.08, 'A', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');
text(ax_bot, -0.04, 1.08, 'B', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

sgtitle(sprintf('Figure 2 -- Neural evidence  (N = %d participants)', n_part), ...
    'FontSize', 12, 'FontWeight', 'bold');

%% ========================================================================
%  6. SAVE
%  ========================================================================

exportgraphics(fig, outPDF, 'ContentType', 'vector', 'BackgroundColor', 'w');
exportgraphics(fig, outPNG, 'Resolution', 300, 'BackgroundColor', 'w');
fprintf('Saved: %s (vector PDF)\n', outPDF);
fprintf('Saved: %s (300 dpi PNG)\n', outPNG);
fprintf('\nDone.\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function [t_sm, sm] = smoothMedian(mat, tTask, winSec, fs)
% smoothMedian  Boxcar-smoothed grand median of (nPart x nSamples).
    med = median(mat, 1);
    win = round(winSec * fs);
    halfK = floor(win / 2);
    kernel = ones(1, win) / win;
    sm  = conv(med, kernel, 'valid');
    t_sm = tTask(halfK+1 : halfK+length(sm));
end

function fit = fitLocked(t_sm, sm, tau, T_trial, E)
% fitLocked  Least-squares fit of S(x) = A*x*exp(-e*x) + B past tau.
    Teff = T_trial - tau;
    mask = t_sm >= tau;
    t_fit = t_sm(mask);
    y_fit = sm(mask);

    x = (t_fit - tau) / Teff;
    s = x .* exp(-E * x);
    X = [s(:), ones(numel(s), 1)];

    coef = X \ y_fit(:);
    A = coef(1);
    B = coef(2);
    yhat = (A * s + B)';

    ss_res = sum((y_fit(:) - yhat(:)).^2);
    ss_tot = sum((y_fit(:) - mean(y_fit)).^2);
    R2 = 1 - ss_res / ss_tot;

    trough_x = 1.0 / E;
    t_trough = tau + trough_x * Teff;

    fit.A        = A;
    fit.B        = B;
    fit.R2       = R2;
    fit.t_trough = t_trough;
    fit.yhat     = yhat;
    fit.t_fit    = t_fit;
end

function cmap = coolwarm_cmap(n)
% Diverging blue-white-red colormap (approximation of matplotlib coolwarm).
    if nargin < 1, n = 256; end
    half = ceil(n/2);
    r1 = linspace(0.23, 0.97, half)';
    g1 = linspace(0.30, 0.96, half)';
    b1 = linspace(0.75, 0.95, half)';
    r2 = linspace(0.97, 0.71, n - half)';
    g2 = linspace(0.96, 0.15, n - half)';
    b2 = linspace(0.95, 0.16, n - half)';
    cmap = [r1 g1 b1; r2 g2 b2];
end

function cmap = viridis_cmap(n)
% Perceptually uniform colormap — dark end truncated for readability.
    if nargin < 1, n = 256; end
    anchors = [ ...
        0.267 0.005 0.329;
        0.283 0.141 0.458;
        0.254 0.265 0.530;
        0.207 0.372 0.553;
        0.164 0.471 0.558;
        0.128 0.567 0.551;
        0.135 0.659 0.518;
        0.267 0.749 0.441;
        0.478 0.821 0.318;
        0.741 0.873 0.150;
        0.993 0.906 0.144];
    xi = linspace(0, 1, size(anchors,1));
    % Sample from 25% onward to skip the near-black region
    xq = linspace(0.25, 1, n);
    cmap = interp1(xi, anchors, xq);
end

function cmap = magma_cmap(n)
% Sequential dark-to-light colormap — dark end truncated for readability.
    if nargin < 1, n = 256; end
    anchors = [ ...
        0.001 0.000 0.014;
        0.101 0.048 0.210;
        0.280 0.084 0.398;
        0.478 0.108 0.430;
        0.676 0.168 0.380;
        0.838 0.282 0.336;
        0.946 0.453 0.356;
        0.988 0.653 0.448;
        0.996 0.838 0.580;
        0.987 0.991 0.750];
    xi = linspace(0, 1, size(anchors,1));
    % Sample from 30% onward to skip the near-black region
    xq = linspace(0.30, 1, n);
    cmap = interp1(xi, anchors, xq);
end

function [locs, matched] = matchToMontage(ch_names, allLocs, allLabelsStd)
% matchToMontage  Build EEGLAB chanlocs for a channel subset.
%   Returns the matched chanlocs struct and a logical mask into ch_names
%   indicating which channels were found in the standard montage.
    matched = false(numel(ch_names), 1);
    locs    = allLocs([]);
    for i = 1:numel(ch_names)
        idx = find(strcmpi(allLabelsStd, ch_names{i}), 1);
        if ~isempty(idx)
            matched(i)  = true;
            locs(end+1) = allLocs(idx); %#ok<AGROW>
        end
    end
end

