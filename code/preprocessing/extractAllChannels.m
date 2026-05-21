%% extractAllChannels.m
% =========================================================================
% Extract Per-Participant, Per-Channel Timecourses from Raw EEG
% =========================================================================
%
% PURPOSE:
%   Extracts per-participant timecourses for all 64 EEG channels from raw
%   .mat files. Pipeline per participant:
%     1. Load each trial file (64 channels x 60000 samples at 1000 Hz)
%     2. Downsample 1000 Hz -> 10 Hz by block-averaging (blocks of 100)
%     3. Take median across trials (robust to single-trial outliers)
%     4. Baseline-correct: subtract mean of first 2 s per channel
%   Participants with fewer than 3 valid trials are excluded.
%
% INPUT:
%   Raw EEG .mat files in: <dataDir>/pceXXYYMMDD/pceXX_P{1,2}_Trial{1-18}.mat
%   Each file contains a 64 x 60000 numeric matrix (float, 1000 Hz).
%
% OUTPUT (saved to data/preprocessed/EEG/allchannel_data.mat):
%   tTask    : [600 x 1]            time axis (s)
%   chans    : {64 x 1} cell        channel label strings
%   allTC    : [nPart x 64 x 600]   baseline-corrected median timecourses
%   dyad_id  : [nPart x 1]          dyad number per participant
%   pid      : [nPart x 1]          participant ID within dyad (1 or 2)
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   May 2026
% =========================================================================

clearvars; close all; clc;

%% ========================================================================
%  1. CONFIGURATION
%  ========================================================================

scriptDir  = fileparts(mfilename('fullpath'));
dataDir    = fullfile(scriptDir, '..', '..', 'data', 'raw', 'EEG');
outputDir  = fullfile(scriptDir, '..', '..', 'data', 'preprocessed', 'EEG');

nChan       = 64;
fsOrig      = 1000;      % Original sampling rate (Hz)
fsOut       = 10;        % Output sampling rate (Hz)
dsRatio     = fsOrig / fsOut;   % = 100
trialDur_s  = 60;        % Trial duration (s)
nSampOut    = trialDur_s * fsOut;   % 600 samples
baselineEnd = 2.0;       % Baseline window [0, baselineEnd) seconds
blIdx       = 1:round(baselineEnd * fsOut);   % baseline sample indices
excludeDyad = 31;        % Dyad 31 excluded (incomplete recording)

% Full 64-channel label order (matching the raw .mat file rows)
chanLabels = { ...
    'Fp1','Fp2','F7','F3','Fz','F4','F8','FC5', ...
    'FC1','FC2','FC6','T7','C3','Cz','C4','T8', ...
    'TP9','CP5','CP1','CP2','CP6','TP10','P7','P3', ...
    'Pz','P4','P8','PO9','O1','Oz','O2','PO10', ...
    'AF7','AF3','AF4','AF8','F5','F1','F2','F6', ...
    'FT9','FT7','FC3','FC4','FT8','FT10','C5','C1', ...
    'C2','C6','TP7','CP3','CPz','CP4','TP8','P5', ...
    'P1','P2','P6','PO7','PO3','POz','PO4','PO8'  ...
};

tTask = ((0:nSampOut-1) / fsOut)';   % [600 x 1] column vector

fprintf('==========================================================\n');
fprintf('  Extract All Channels: Per-Participant Timecourses\n');
fprintf('==========================================================\n');
fprintf('  Data directory:    %s\n', dataDir);
fprintf('  Output directory:  %s\n', outputDir);
fprintf('  Channels:          %d\n', nChan);
fprintf('  Downsample:        %d Hz -> %d Hz\n', fsOrig, fsOut);
fprintf('  Baseline:          [0, %.1f) s\n', baselineEnd);
fprintf('==========================================================\n\n');

%% ========================================================================
%  2. DISCOVER DYAD FOLDERS
%  ========================================================================

allEntries  = dir(dataDir);
allEntries  = allEntries([allEntries.isdir]);
dyadFolders = struct('dyadNum',{},'folderName',{},'folderPath',{});

for i = 1:length(allEntries)
    tok = regexp(allEntries(i).name, '^pce(\d{2})\d{6}$', 'tokens');
    if ~isempty(tok)
        dN = str2double(tok{1}{1});
        if dN == excludeDyad, continue; end
        e.dyadNum    = dN;
        e.folderName = allEntries(i).name;
        e.folderPath = fullfile(dataDir, allEntries(i).name);
        dyadFolders(end+1) = e; %#ok<SAGROW>
    end
end
[~, si]     = sort([dyadFolders.dyadNum]);
dyadFolders = dyadFolders(si);
nDyads      = length(dyadFolders);
fprintf('Found %d dyad folders (excluding dyad %d)\n\n', nDyads, excludeDyad);

%% ========================================================================
%  3. EXTRACT AND AGGREGATE
%  ========================================================================

allTC_cell  = {};    % each entry: [64 x 600]
dyad_id_vec = [];
pid_vec     = [];
nPartOK     = 0;
nSkipped    = 0;
tStart      = tic;

for di = 1:nDyads
    df = dyadFolders(di);
    fprintf('[%2d/%2d] Dyad %02d (%s)\n', di, nDyads, df.dyadNum, df.folderName);

    for p = [1 2]
        fprintf('  P%d: ', p);

        % Pre-allocate trial stack: [nTrials x 64 x 600], filled with NaN
        trialStack = nan(18, nChan, nSampOut);
        nLoaded = 0;

        for ti = 1:18
            fn = sprintf('pce%02d_P%d_Trial%d.mat', df.dyadNum, p, ti);
            fp = fullfile(df.folderPath, fn);

            if ~exist(fp, 'file'), continue; end

            try
                raw = loadRawEEG_allChan(fp, nChan);
                if isempty(raw), continue; end

                % Clip to at most 60 s and trim to whole blocks
                nFull = min(size(raw, 2), trialDur_s * fsOrig);
                nBins = floor(nFull / dsRatio);
                raw   = raw(:, 1:nBins * dsRatio);

                % Downsample: block-average [64 x dsRatio x nBins] -> [64 x nBins]
                ds = squeeze(mean(reshape(raw, nChan, dsRatio, nBins), 2));

                nOut = min(nBins, nSampOut);
                trialStack(ti, :, 1:nOut) = ds(:, 1:nOut);
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

        % Median across trials -> [64 x 600]
        med = squeeze(median(trialStack(1:18, :, :), 1, 'omitnan'));

        % Baseline-correct per channel: subtract mean of [0, 2) s
        baseMean = mean(med(:, blIdx), 2, 'omitnan');   % [64 x 1]
        med = med - baseMean;

        % Store
        nPartOK = nPartOK + 1;
        allTC_cell{nPartOK}    = med;          %#ok<SAGROW>
        dyad_id_vec(nPartOK)   = df.dyadNum;   %#ok<SAGROW>
        pid_vec(nPartOK)       = p;             %#ok<SAGROW>

        fprintf('%d trials  (elapsed %.1f s)\n', nLoaded, toc(tStart));
    end
end

%% ========================================================================
%  4. ASSEMBLE AND SAVE
%  ========================================================================

if ~exist(outputDir, 'dir'), mkdir(outputDir); end

% Stack into 3-D array: [nPart x 64 x 600]
allTC   = zeros(nPartOK, nChan, nSampOut);
for k = 1:nPartOK
    allTC(k, :, :) = allTC_cell{k};
end

dyad_id = dyad_id_vec(:);   % column vector
pid     = pid_vec(:);        % column vector
chans   = chanLabels(:);     % {64 x 1} cell array of strings

outFile = fullfile(outputDir, 'allchannel_data.mat');
save(outFile, 'tTask', 'chans', 'allTC', 'dyad_id', 'pid', '-v7.3');

fprintf('\n==========================================================\n');
fprintf('  Saved: %s\n', outFile);
fprintf('  allTC shape: [%d x %d x %d]\n', size(allTC));
fprintf('  Participants: %d  (skipped: %d)\n', nPartOK, nSkipped);
fprintf('  Total time: %.1f s\n', toc(tStart));
fprintf('==========================================================\n');


%% ========================================================================
%  LOCAL FUNCTION
%  ========================================================================

function raw = loadRawEEG_allChan(filepath, expectedChan)
%LOADRAWEEG_ALLCHAN  Load the first 2-D numeric variable with expectedChan rows.
%   Returns double matrix [expectedChan x nSamples], or [] if not found.
    S = load(filepath);
    fnames = fieldnames(S);
    raw = [];
    for fi = 1:length(fnames)
        v = S.(fnames{fi});
        if isnumeric(v) && ismatrix(v) && size(v,1) == expectedChan
            raw = double(v);
            return;
        end
    end
end
