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

% Full 64-channel label order MATCHING THE RAW .mat FILE ROWS.
% ------------------------------------------------------------------------
% CORRECTED 2026-07-13 (montage bug). The OSF/dataset per-trial `.mat` rows are
% in the dataset authors' ANATOMICAL order (documented by the shipped
% `actiCAP-10-20-Cap64.ced`), NOT the BrainProducts/actiCAP acquisition order.
% The previous array here was a manufacturer-default template — wrong ORDER
% (61/64 rows mislabelled) AND wrong MEMBERSHIP (it listed PO9/PO10; this cap
% recorded AFz/Iz — verified from the raw BrainVision `.vhdr` montage
% CMA-64_REF_HS2.bvef and from signal physiology: alpha posterior, blinks at
% AFz, midline Iz). Using the old array scrambled every per-channel/topographic
% result (e.g. Figure2_perchannel_fits.csv -> Fig 4B). The channel-MEAN GSP
% (Panel A) and 1/f exponent (Panel C) are permutation-invariant and were
% unaffected. Full provenance: pce-master-loop docs/2026-07-13-montage-cap-provenance.md.
chanLabels = { ...
    'Fp1','Fp2','AF7','AF3','AFz','AF4','AF8','F7', ...
    'F5','F3','F1','Fz','F2','F4','F6','F8', ...
    'FT9','FT7','FC5','FC3','FC1','FC2','FC4','FC6', ...
    'FT8','FT10','T7','C5','C3','C1','Cz','C2', ...
    'C4','C6','T8','TP9','TP7','CP5','CP3','CP1', ...
    'CPz','CP2','CP4','CP6','TP8','TP10','P7','P5', ...
    'P3','P1','Pz','P2','P4','P6','P8','PO7', ...
    'PO3','POz','PO4','PO8','O1','Oz','O2','Iz'  ...
};
% GUARDRAIL: this cap recorded AFz/Iz, never PO9/PO10 — fail loud if that ever
% regresses to a manufacturer-default template.
assert(any(strcmp(chanLabels,'AFz')) && any(strcmp(chanLabels,'Iz')) && ...
       ~any(strcmp(chanLabels,'PO9')) && ~any(strcmp(chanLabels,'PO10')) && ...
       numel(chanLabels)==64, ...
       'extractAllChannels:badMontage', ...
       ['chanLabels must be the .ced anatomical order with AFz/Iz (not PO9/PO10). ' ...
        'See docs/2026-07-13-montage-cap-provenance.md.']);

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

% GUARDRAIL: if no raw dyad folders are found, FAIL rather than silently writing
% an empty allchannel_data.mat (which would then be "fixed" from the wrong source
% -- e.g. the 1-40 Hz FIF, see dossier 2026-06-22). The raw .mat live in the
% master loop: pce-master-loop/data/raw/EEG (or fetch from OSF). On Tom's machine
% the spoke path data/raw/EEG is a symlink to it.
if nDyads == 0
    error('extractAllChannels:noRawData', ...
        ['No raw dyad folders (pceXXYYMMDD/) under:\n  %s\n' ...
         'Raw EEG lives in pce-master-loop/data/raw/EEG (symlink it here) or fetch\n' ...
         'from OSF. This unfiltered raw is REQUIRED for the slow-trend GSP/per-channel\n' ...
         'fits; the 1-40 Hz *-raw.fif is NOT a substitute. (dossier 2026-06-22)'], dataDir);
end

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
