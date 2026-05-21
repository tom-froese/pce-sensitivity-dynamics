%% plotFigureS2_TopoMap.m
% =========================================================================
% Four-panel topographic map of per-channel S(x) fits (EEGLAB topoplot)
% =========================================================================
%
% Generates the supplementary topomap figure from per-channel fit results.
% Requires the CSV output of computePerChannelFits.m and the allchannel
% data for channel labels.
%
% ARCHIVED — PNAS Brief Report format does not include supplementary
% figures.  Retained here for potential use in journal resubmission or
% extended manuscript.
%
% INPUTS:
%   ../../results/FigureS2_GSP_TopoMap_FreeTau_perchannel.csv
%   ../../data/preprocessed/EEG/allchannel_data.mat  (for channel labels)
%
% OUTPUTS:
%   ../../results/FigureS2_GSP_TopoMap_FreeTau.pdf
%
% DEPENDENCIES: EEGLAB (topoplot function)
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   May 2026
% =========================================================================

close all;

scriptDir = fileparts(mfilename('fullpath'));
ROOT      = fullfile(scriptDir, '..', '..');
csvFile   = fullfile(ROOT, 'results', 'FigureS2_GSP_TopoMap_FreeTau_perchannel.csv');
dataFile  = fullfile(ROOT, 'data', 'preprocessed', 'EEG', 'allchannel_data.mat');
outPDF    = fullfile(ROOT, 'results', 'FigureS2_GSP_TopoMap_FreeTau.pdf');

TAU_LOCKED = 3.90;

if ~isfile(csvFile)
    error('Missing: %s\nRun computePerChannelFits first.', csvFile);
end

%% Load per-channel fit results
T = readtable(csvFile);

tau_star_vec   = T.tau_star;
R2_free_vec    = T.R2_free;
R2_lock_vec    = T.R2_lock;
dR2_vec        = T.dR2;
trough_amp_vec = T.trough_amp_free;

% Get channel labels from data file
d     = load(dataFile, 'chans', 'allTC');
chans = d.chans;
n_part = size(d.allTC, 1);
n_chan = numel(chans);

%% EEGLAB channel locations
if ~exist('topoplot', 'file')
    eeglabDir = fullfile(getenv('HOME'), 'Documents', 'MATLAB', 'eeglab2026.0.0');
    if isfolder(eeglabDir), addpath(eeglabDir); eeglab('nogui');
    else, error('EEGLAB not found. Install to ~/Documents/MATLAB/eeglab2026.0.0/');
    end
end
locFile   = fullfile(fileparts(which('eeglab')), 'plugins', 'dipfit', ...
    'standard_BEM', 'elec', 'standard_1005.elc');
allLocs   = readlocs(locFile);
allLabels = {allLocs.labels};
chanlocs  = allLocs([]);
for i = 1:n_chan
    idx = find(strcmpi(allLabels, chans{i}), 1);
    if ~isempty(idx), chanlocs(end+1) = allLocs(idx); end %#ok<SAGROW>
end
fprintf('Matched %d/%d channels to standard 10-05 montage\n', length(chanlocs), n_chan);

%% Four-panel topographic figure

fprintf('Creating topomap figure ...\n');

fig = figure('Units', 'inches', 'Position', [0.5 0.5 18.0 4.8], ...
    'Color', 'w', 'PaperUnits', 'inches', ...
    'PaperSize', [18.0 4.8], 'PaperPosition', [0 0 18.0 4.8]);

% --- Panel 1: tau* (diverging around 3.9 s) ---
ax1 = subplot(1, 4, 1);
topoplot(tau_star_vec, chanlocs, 'maplimits', [0 12], ...
    'numcontour', 6, 'electrodes', 'on', 'style', 'fill', 'conv', 'on');
colormap(ax1, coolwarm(256));
title('\tau^{*}  (s, free fit)', 'FontSize', 10);
cb1 = colorbar('southoutside'); cb1.Label.String = '\tau^{*} (s)';

% --- Panel 2: R^2 (free) ---
ax2 = subplot(1, 4, 2);
topoplot(R2_free_vec, chanlocs, 'maplimits', [0 max(R2_free_vec)], ...
    'numcontour', 6, 'electrodes', 'on', 'style', 'fill', 'conv', 'on');
colormap(ax2, viridis(256));
title('R^2 at \tau^{*}', 'FontSize', 10);
cb2 = colorbar('southoutside'); cb2.Label.String = 'R^2';

% --- Panel 3: delta-R^2 (free - locked) ---
ax3 = subplot(1, 4, 3);
dR2_max = max(prctile(abs(dR2_vec), 95), 0.01);
topoplot(dR2_vec, chanlocs, 'maplimits', [-dR2_max dR2_max], ...
    'numcontour', 6, 'electrodes', 'on', 'style', 'fill', 'conv', 'on');
colormap(ax3, PiYG(256));
title(sprintf('\\DeltaR^2 (free - lock, \\tau=%.2f)', TAU_LOCKED), 'FontSize', 10);
cb3 = colorbar('southoutside'); cb3.Label.String = '\DeltaR^2';

% --- Panel 4: |trough amp| at tau* ---
ax4 = subplot(1, 4, 4);
abs_trough = abs(trough_amp_vec);
abs_max = prctile(abs_trough, 90);
topoplot(abs_trough, chanlocs, 'maplimits', [0 abs_max], ...
    'numcontour', 6, 'electrodes', 'on', 'style', 'fill', 'conv', 'on');
colormap(ax4, flipud(magma(256)));
title(sprintf('|A s(1/e)| at \\tau^{*} (\\muV, 90th=%.0f)', abs_max), 'FontSize', 10);
cb4 = colorbar('southoutside'); cb4.Label.String = '\muV';

sgtitle(sprintf('Per-electrode S(x) fit with free \\tau  (grid 0-15 s @ 0.1 s; n = %d participants)', ...
    n_part), 'FontSize', 12);

%% Save figure

exportgraphics(fig, outPDF, 'ContentType', 'vector');
fprintf('Saved: %s\n', outPDF);
fprintf('\nDone.\n');

%% ========================================================================
%  LOCAL COLORMAP FUNCTIONS
%  ========================================================================

function cmap = coolwarm(n)
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

function cmap = viridis(n)
    if nargin < 1, n = 256; end
    anchors = [ ...
        0.267 0.005 0.329; 0.283 0.141 0.458; 0.254 0.265 0.530;
        0.207 0.372 0.553; 0.164 0.471 0.558; 0.128 0.567 0.551;
        0.135 0.659 0.518; 0.267 0.749 0.441; 0.478 0.821 0.318;
        0.741 0.873 0.150; 0.993 0.906 0.144];
    xi = linspace(0, 1, size(anchors,1));
    xq = linspace(0, 1, n);
    cmap = interp1(xi, anchors, xq);
end

function cmap = PiYG(n)
    if nargin < 1, n = 256; end
    anchors = [ ...
        0.56 0.00 0.42; 0.78 0.35 0.62; 0.91 0.60 0.80;
        0.97 0.80 0.91; 0.97 0.96 0.97; 0.85 0.94 0.68;
        0.65 0.86 0.42; 0.40 0.74 0.24; 0.15 0.58 0.10];
    xi = linspace(0, 1, size(anchors,1));
    xq = linspace(0, 1, n);
    cmap = interp1(xi, anchors, xq);
end

function cmap = magma(n)
    if nargin < 1, n = 256; end
    anchors = [ ...
        0.001 0.000 0.014; 0.101 0.048 0.210; 0.280 0.084 0.398;
        0.478 0.108 0.430; 0.676 0.168 0.380; 0.838 0.282 0.336;
        0.946 0.453 0.356; 0.988 0.653 0.448; 0.996 0.838 0.580;
        0.987 0.991 0.750];
    xi = linspace(0, 1, size(anchors,1));
    xq = linspace(0, 1, n);
    cmap = interp1(xi, anchors, xq);
end
