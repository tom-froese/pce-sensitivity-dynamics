%% globalScalpPotential.m
% =========================================================================
% Global Scalp Potential: Hierarchical Extraction from Raw EEG
% =========================================================================
%
% PURPOSE:
%   Extracts the global scalp potential (mean across all 64 EEG channels)
%   from raw .mat files, using a proper hierarchical aggregation:
%     1. Downsample each trial from 1000 Hz to 10 Hz (block-averaging)
%     2. Average across all 64 channels at each time point (global potential)
%     3. Average (median) across trials within each participant
%     4. Baseline-correct: subtract mean of first 2 s
%
%   This script produces participant-level time courses that can then be
%   analyzed by globalScalpPotential_stats.m and plotted by
%   globalScalpPotential_figures.m.
%
%   No common average re-referencing is applied, preserving the global
%   scalp potential (which would be zeroed by CAR). No filtering beyond
%   the implicit low-pass from downsampling. Median aggregation across
%   trials is robust to single-trial outliers.
%
% INPUT:
%   Raw EEG .mat files in: <dataDir>/pceXXYYMMDD/pceXX_P{1,2}_Trial{1-18}.mat
%   Each file: 64 channels x 60000 samples (float32, 1000 Hz)
%
%   Optional: Rest files pceXX_P{1,2}_Rest{1-4}.mat (180s, 64 x 180000)
%
% OUTPUT (saved to <outputDir>):
%   globalScalpPotential_data.mat containing:
%     - taskParticipantTC:  [nParticipants x nSampTask] participant-level
%                           baseline-corrected median time courses (task)
%     - restParticipantTC:  [nParticipants x nSampRest] (rest)
%     - participantInfo:    table with DyadID, ParticipantID, nTrials, etc.
%     - roiParticipantTC:   {nROI x 1} cell, each [nPart x nSampTask]
%     - tTask, tRest:       time vectors (s)
%     - cfg:                configuration struct
%     - roi:                ROI definitions
%
% USAGE:
%   Run from the Matlab scripts directory, or set paths in cfg below.
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   April 2026
% =========================================================================

clearvars; close all; clc;

%% ========================================================================
%  1. CONFIGURATION
%  ========================================================================

scriptDir      = fileparts(mfilename('fullpath'));
cfg.dataDir    = fullfile(scriptDir, '..', '..', 'data', 'raw', 'EEG');
cfg.outputDir  = fullfile(scriptDir, '..', '..', 'data', 'preprocessed', 'EEG');

cfg.nChan      = 64;        % Total channels in raw files
cfg.fsOrig     = 1000;      % Original sampling rate (Hz)
cfg.fsOut      = 10;        % Output sampling rate (Hz)
cfg.dsRatio    = cfg.fsOrig / cfg.fsOut;  % = 100

cfg.trialDur_s = 60;        % Trial duration (s)
cfg.restDur_s  = 180;       % Rest duration (s)
cfg.baselineEnd = 2.0;      % Baseline window [0, baselineEnd) seconds

cfg.excludeDyad = 31;       % Dyad 31 excluded (incomplete recording)

cfg.participants = [1, 2];
cfg.trials       = 1:18;
cfg.restBlocks   = 1:4;

fprintf('==========================================================\n');
fprintf('  Global Scalp Potential: Hierarchical Extraction\n');
fprintf('==========================================================\n');
fprintf('  Data directory:   %s\n', cfg.dataDir);
fprintf('  Output directory:  %s\n', cfg.outputDir);
fprintf('  Channels:          %d (all)\n', cfg.nChan);
fprintf('  Downsample:        %d Hz -> %d Hz\n', cfg.fsOrig, cfg.fsOut);
fprintf('  Trial duration:    %d s\n', cfg.trialDur_s);
fprintf('  Baseline:          [0, %.1f) s\n', cfg.baselineEnd);
fprintf('==========================================================\n\n');

%% ========================================================================
%  2. ROI DEFINITIONS (64-channel extended 10-10 montage)
%  ========================================================================
%  Channel order in the raw .mat files follows the actiCAP snap layout.
%  We define ROIs by label and map to hardware indices below.

% Full 64-channel label order (matching the raw .mat file rows)
cfg.chanLabels = { ...
    'Fp1','Fp2','F7','F3','Fz','F4','F8','FC5', ...
    'FC1','FC2','FC6','T7','C3','Cz','C4','T8', ...
    'TP9','CP5','CP1','CP2','CP6','TP10','P7','P3', ...
    'Pz','P4','P8','PO9','O1','Oz','O2','PO10', ...
    'AF7','AF3','AF4','AF8','F5','F1','F2','F6', ...
    'FT9','FT7','FC3','FC4','FT8','FT10','C5','C1', ...
    'C2','C6','TP7','CP3','CPz','CP4','TP8','P5', ...
    'P1','P2','P6','PO7','PO3','POz','PO4','PO8'  ...
};

roi.names  = {'Prefrontal','Frontal','Fronto-Central','Central', ...
              'Centro-Parietal','Parietal','Occipital'};
roi.chans  = { ...
    {'Fp1','Fp2','AF3','AF4','AF7','AF8'}, ...
    {'F7','F3','Fz','F4','F8','F1','F2','F5','F6'}, ...
    {'FC5','FC1','FC2','FC6','FC3','FC4','FT7','FT8','FT9','FT10'}, ...
    {'C3','Cz','C4','C1','C2','C5','C6','T7','T8'}, ...
    {'CP5','CP1','CP2','CP6','CP3','CP4','CPz','TP7','TP8','TP9','TP10'}, ...
    {'P7','P3','Pz','P4','P8','P1','P2','P5','P6'}, ...
    {'O1','Oz','O2','PO3','PO4','PO7','PO8','PO9','PO10','POz'} ...
};
nROI = length(roi.names);

% Map ROI channel labels to row indices in the 64-channel raw data
roi.indices = cell(nROI, 1);
for ri = 1:nROI
    idx = [];
    for ci = 1:length(roi.chans{ri})
        match = find(strcmp(cfg.chanLabels, roi.chans{ri}{ci}));
        if ~isempty(match)
            idx = [idx, match]; %#ok<AGROW>
        else
            warning('Channel %s not found in chanLabels', roi.chans{ri}{ci});
        end
    end
    roi.indices{ri} = idx;
    fprintf('  ROI %-18s: %2d channels\n', roi.names{ri}, length(idx));
end
fprintf('\n');

%% ========================================================================
%  3. SETUP
%  ========================================================================

if ~exist(cfg.outputDir, 'dir'), mkdir(cfg.outputDir); end

nSampTask = cfg.trialDur_s * cfg.fsOut;   % 600 samples
nSampRest = cfg.restDur_s  * cfg.fsOut;   % 1800 samples
blIdx     = 1:round(cfg.baselineEnd * cfg.fsOut);  % baseline indices

tTask = (0:nSampTask-1) / cfg.fsOut;
tRest = (0:nSampRest-1) / cfg.fsOut;

%% ========================================================================
%  4. DISCOVER DYAD FOLDERS
%  ========================================================================

allEntries = dir(cfg.dataDir);
allEntries = allEntries([allEntries.isdir]);
dyadFolders = struct('dyadNum',{},'folderName',{},'folderPath',{});

for i = 1:length(allEntries)
    tok = regexp(allEntries(i).name, '^pce(\d{2})\d{6}$', 'tokens');
    if ~isempty(tok)
        dN = str2double(tok{1}{1});
        if dN == cfg.excludeDyad, continue; end
        e.dyadNum    = dN;
        e.folderName = allEntries(i).name;
        e.folderPath = fullfile(cfg.dataDir, allEntries(i).name);
        dyadFolders(end+1) = e; %#ok<SAGROW>
    end
end
[~, si] = sort([dyadFolders.dyadNum]);
dyadFolders = dyadFolders(si);
nDyads = length(dyadFolders);
fprintf('Found %d dyad folders (excluding dyad %d)\n\n', nDyads, cfg.excludeDyad);

%% ========================================================================
%  5. EXTRACT AND AGGREGATE
%  ========================================================================
%  For each participant:
%    - Load all trial files, downsample to 10 Hz, compute global mean
%    - Take MEDIAN across trials -> robust participant-level time course
%    - Baseline-correct: subtract mean of [0, 2) s
%    - Also compute per-ROI averages

% Pre-allocate storage
taskGlobal   = {};   % will become [nPart x nSampTask]
restGlobal   = {};   % will become [nPart x nSampRest]
taskROI      = cell(nROI, 1);  % {roi}{partIdx} -> vec
for ri = 1:nROI, taskROI{ri} = {}; end

partInfo = table('Size', [0 5], ...
    'VariableTypes', {'double','double','double','double','double'}, ...
    'VariableNames', {'DyadID','ParticipantID','nTrialsLoaded','nRestLoaded','nTrialsAvailable'});

nPartOK = 0;  nSkipped = 0;
tStart = tic;

for di = 1:nDyads
    df = dyadFolders(di);
    fprintf('[%2d/%2d] Dyad %02d (%s)\n', di, nDyads, df.dyadNum, df.folderName);

    for p = cfg.participants
        fprintf('  P%d: ', p);
        pTic = tic;

        % --- Load all TASK trials ---
        trialVecs = nan(length(cfg.trials), nSampTask);  % global
        roiVecs   = cell(nROI, 1);
        for ri = 1:nROI
            roiVecs{ri} = nan(length(cfg.trials), nSampTask);
        end
        nLoaded = 0;

        for ti = 1:length(cfg.trials)
            tNum = cfg.trials(ti);
            fn = sprintf('pce%02d_P%d_Trial%d.mat', df.dyadNum, p, tNum);
            fp = fullfile(df.folderPath, fn);

            if ~exist(fp, 'file'), continue; end

            try
                raw = loadRawEEG(fp, cfg.nChan);
                if isempty(raw), continue; end

                % Downsample: block-average 100 samples -> 1 sample
                nFull = min(size(raw, 2), cfg.trialDur_s * cfg.fsOrig);
                raw   = raw(:, 1:nFull);
                nBins = floor(nFull / cfg.dsRatio);
                raw   = raw(:, 1:nBins * cfg.dsRatio);
                ds    = squeeze(mean(reshape(raw, cfg.nChan, cfg.dsRatio, nBins), 2));
                % ds is [64 x nBins]

                nOut = min(size(ds, 2), nSampTask);

                % Global: mean across all 64 channels
                trialVecs(ti, 1:nOut) = mean(ds(:, 1:nOut), 1);

                % Per ROI: mean across ROI channels
                for ri = 1:nROI
                    roiVecs{ri}(ti, 1:nOut) = mean(ds(roi.indices{ri}, 1:nOut), 1);
                end

                nLoaded = nLoaded + 1;
            catch ME
                fprintf('[ERR:%s] ', ME.message);
            end
        end

        if nLoaded < 3
            fprintf('skipped (only %d trials)\n', nLoaded);
            nSkipped = nSkipped + 1;
            continue;
        end

        % Median across trials -> participant-level time course
        pGlobal = median(trialVecs, 1, 'omitnan');

        % Baseline-correct: subtract mean of [0, 2) s
        pGlobal = pGlobal - mean(pGlobal(blIdx), 'omitnan');

        taskGlobal{end+1} = pGlobal; %#ok<SAGROW>

        % Per ROI
        for ri = 1:nROI
            rMed = median(roiVecs{ri}, 1, 'omitnan');
            rMed = rMed - mean(rMed(blIdx), 'omitnan');
            taskROI{ri}{end+1} = rMed;
        end

        % --- Load REST blocks ---
        restVecs = [];
        nRestLoaded = 0;
        for r = cfg.restBlocks
            fn = sprintf('pce%02d_P%d_Rest%d.mat', df.dyadNum, p, r);
            fp = fullfile(df.folderPath, fn);
            if ~exist(fp, 'file'), continue; end

            try
                raw = loadRawEEG(fp, cfg.nChan);
                if isempty(raw), continue; end

                nFull = min(size(raw, 2), cfg.restDur_s * cfg.fsOrig);
                raw   = raw(:, 1:nFull);
                nBins = floor(nFull / cfg.dsRatio);
                raw   = raw(:, 1:nBins * cfg.dsRatio);
                ds    = squeeze(mean(reshape(raw, cfg.nChan, cfg.dsRatio, nBins), 2));

                nOut = min(size(ds, 2), nSampRest);
                vec  = nan(1, nSampRest);
                vec(1:nOut) = mean(ds(:, 1:nOut), 1);
                restVecs = [restVecs; vec]; %#ok<AGROW>
                nRestLoaded = nRestLoaded + 1;
            catch
                % skip silently
            end
        end

        if nRestLoaded > 0
            pRest = median(restVecs, 1, 'omitnan');
            pRest = pRest - mean(pRest(blIdx), 'omitnan');
            restGlobal{end+1} = pRest; %#ok<SAGROW>
        end

        nPartOK = nPartOK + 1;
        partInfo(end+1, :) = {df.dyadNum, p, nLoaded, nRestLoaded, length(cfg.trials)};

        fprintf('%d trials, %d rest (%.1f s)\n', nLoaded, nRestLoaded, toc(pTic));
    end
end
elapsedMin = toc(tStart) / 60;

%% ========================================================================
%  6. ASSEMBLE MATRICES
%  ========================================================================

taskParticipantTC = vertcat(taskGlobal{:});   % [nPart x nSampTask]
nPart = size(taskParticipantTC, 1);

if ~isempty(restGlobal)
    restParticipantTC = vertcat(restGlobal{:});
else
    restParticipantTC = [];
end

roiParticipantTC = cell(nROI, 1);
for ri = 1:nROI
    roiParticipantTC{ri} = vertcat(taskROI{ri}{:});
end

fprintf('\n=== EXTRACTION COMPLETE ===\n');
fprintf('  Participants (task):  %d\n', nPart);
fprintf('  Participants (rest):  %d\n', size(restParticipantTC, 1));
fprintf('  Skipped:              %d\n', nSkipped);
fprintf('  Elapsed:              %.1f min\n', elapsedMin);

%% ========================================================================
%  7. SAVE
%  ========================================================================

outputFile = fullfile(cfg.outputDir, 'globalScalpPotential_data.mat');
save(outputFile, ...
    'taskParticipantTC', 'restParticipantTC', ...
    'roiParticipantTC', 'partInfo', ...
    'tTask', 'tRest', 'nSampTask', 'nSampRest', ...
    'cfg', 'roi', 'nROI', 'nPart', ...
    '-v7.3');
fprintf('  Saved: %s\n', outputFile);
fprintf('  Matrix size (task): %d x %d\n', size(taskParticipantTC));
fprintf('Done.\n');

%% ========================================================================
%  LOCAL FUNCTION: loadRawEEG
%  ========================================================================

function raw = loadRawEEG(fpath, nCh)
%LOADRAWEEG  Load a raw EEG .mat file and return [nCh x nSamples].
%   Searches for a numeric matrix with the expected number of channels.
    S = load(fpath);
    fns = fieldnames(S);
    raw = [];
    for i = 1:length(fns)
        v = S.(fns{i});
        if isnumeric(v) && ismatrix(v)
            if size(v, 1) == nCh
                raw = double(v);
                return;
            elseif size(v, 2) == nCh
                raw = double(v');
                return;
            end
        end
    end
    % If no exact match, try the first numeric matrix
    for i = 1:length(fns)
        v = S.(fns{i});
        if isnumeric(v) && ismatrix(v) && min(size(v)) > 1
            if size(v, 1) < size(v, 2)
                raw = double(v);  % assume channels x samples
            else
                raw = double(v');
            end
            return;
        end
    end
end
