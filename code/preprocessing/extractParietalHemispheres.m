%% extractParietalHemispheres.m
% =========================================================================
% Extract Per-Participant Parietal Hemisphere Timecourses from Raw EEG
% =========================================================================
%
% PURPOSE:
%   Extracts per-participant timecourses for three parietal sub-regions
%   (Left, Midline, Right) from raw EEG .mat files. Pipeline per
%   participant:
%     1. Load each trial file (64 channels x 60000 samples at 1000 Hz)
%     2. Select parietal channels, average within each hemisphere cluster
%     3. Downsample 1000 Hz -> 10 Hz by block-averaging (blocks of 100)
%     4. Take median across trials (robust to single-trial outliers)
%     5. Baseline-correct: subtract mean of first 2 s
%   Participants with fewer than 3 valid trials are excluded.
%
% CHANNEL DEFINITIONS (1-based indices into the 64-channel montage):
%   Parietal-Left  (4 chans): P7 (23), P3 (24), P5 (56), P1 (57)
%   Parietal-Mid   (1 chan):  Pz (25)
%   Parietal-Right (4 chans): P4 (26), P8 (27), P2 (58), P6 (59)
%
% INPUT:
%   Raw EEG .mat files in: <dataDir>/pceXXYYMMDD/pceXX_P{1,2}_Trial{1-18}.mat
%   Each file contains a 64 x 60000 numeric matrix (float, 1000 Hz).
%
% OUTPUT (saved to data/preprocessed/EEG/parietal_hemisphere_data.mat):
%   tTask    : [600 x 1]        time axis (s)
%   left_TC  : [nPart x 600]    left-parietal mean timecourses
%   mid_TC   : [nPart x 600]    midline-parietal (Pz) timecourses
%   right_TC : [nPart x 600]    right-parietal mean timecourses
%   dyad_id  : [nPart x 1]      dyad number per participant
%   pid      : [nPart x 1]      participant ID within dyad (1 or 2)
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

% Parietal channel indices (1-based) into the 64-channel montage
leftIdx  = [23, 24, 56, 57];   % P7, P3, P5, P1
midIdx   = [25];               % Pz
rightIdx = [26, 27, 58, 59];   % P4, P8, P2, P6

tTask = ((0:nSampOut-1) / fsOut)';   % [600 x 1] column vector

fprintf('==========================================================\n');
fprintf('  Extract Parietal Hemispheres: Per-Participant Timecourses\n');
fprintf('==========================================================\n');
fprintf('  Data directory:    %s\n', dataDir);
fprintf('  Output directory:  %s\n', outputDir);
fprintf('  Downsample:        %d Hz -> %d Hz\n', fsOrig, fsOut);
fprintf('  Left channels:     P7, P3, P5, P1 (indices %s)\n', mat2str(leftIdx));
fprintf('  Mid channel:       Pz (index %d)\n', midIdx);
fprintf('  Right channels:    P4, P8, P2, P6 (indices %s)\n', mat2str(rightIdx));
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

left_cell   = {};
mid_cell    = {};
right_cell  = {};
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

        % Pre-allocate trial stacks: [nTrials x 600], filled with NaN
        leftTrials  = nan(18, nSampOut);
        midTrials   = nan(18, nSampOut);
        rightTrials = nan(18, nSampOut);
        nLoaded = 0;

        for ti = 1:18
            fn = sprintf('pce%02d_P%d_Trial%d.mat', df.dyadNum, p, ti);
            fp = fullfile(df.folderPath, fn);

            if ~exist(fp, 'file'), continue; end

            try
                raw = loadRawEEG_parietal(fp, nChan);
                if isempty(raw), continue; end

                % Clip to at most 60 s and trim to whole blocks
                nFull = min(size(raw, 2), trialDur_s * fsOrig);
                nBins = floor(nFull / dsRatio);
                raw   = raw(:, 1:nBins * dsRatio);

                % Downsample: block-average [64 x dsRatio x nBins] -> [64 x nBins]
                ds = squeeze(mean(reshape(raw, nChan, dsRatio, nBins), 2));

                nOut = min(nBins, nSampOut);

                % Average within hemisphere clusters, then store
                leftTrials(ti, 1:nOut)  = mean(ds(leftIdx,  1:nOut), 1);
                midTrials(ti, 1:nOut)   = mean(ds(midIdx,   1:nOut), 1);
                rightTrials(ti, 1:nOut) = mean(ds(rightIdx, 1:nOut), 1);
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

        % Median across trials -> [1 x 600]
        leftMed  = median(leftTrials,  1, 'omitnan');
        midMed   = median(midTrials,   1, 'omitnan');
        rightMed = median(rightTrials, 1, 'omitnan');

        % Baseline-correct: subtract mean of [0, 2) s
        leftMed  = leftMed  - mean(leftMed(blIdx),  'omitnan');
        midMed   = midMed   - mean(midMed(blIdx),   'omitnan');
        rightMed = rightMed - mean(rightMed(blIdx),  'omitnan');

        % Store
        nPartOK = nPartOK + 1;
        left_cell{nPartOK}     = leftMed;       %#ok<SAGROW>
        mid_cell{nPartOK}      = midMed;        %#ok<SAGROW>
        right_cell{nPartOK}    = rightMed;       %#ok<SAGROW>
        dyad_id_vec(nPartOK)   = df.dyadNum;    %#ok<SAGROW>
        pid_vec(nPartOK)       = p;              %#ok<SAGROW>

        fprintf('%d trials  (elapsed %.1f s)\n', nLoaded, toc(tStart));
    end
end

%% ========================================================================
%  4. ASSEMBLE AND SAVE
%  ========================================================================

if ~exist(outputDir, 'dir'), mkdir(outputDir); end

% Stack into 2-D arrays: [nPart x 600]
left_TC  = vertcat(left_cell{:});
mid_TC   = vertcat(mid_cell{:});
right_TC = vertcat(right_cell{:});

dyad_id = dyad_id_vec(:);   % column vector
pid     = pid_vec(:);        % column vector

outFile = fullfile(outputDir, 'parietal_hemisphere_data.mat');
save(outFile, 'tTask', 'left_TC', 'mid_TC', 'right_TC', 'dyad_id', 'pid', '-v7.3');

fprintf('\n==========================================================\n');
fprintf('  Saved: %s\n', outFile);
fprintf('  left_TC:  [%d x %d]\n', size(left_TC));
fprintf('  mid_TC:   [%d x %d]\n', size(mid_TC));
fprintf('  right_TC: [%d x %d]\n', size(right_TC));
fprintf('  Participants: %d  (skipped: %d)\n', nPartOK, nSkipped);
fprintf('  Total time: %.1f s\n', toc(tStart));
fprintf('==========================================================\n');


%% ========================================================================
%  LOCAL FUNCTION
%  ========================================================================

function raw = loadRawEEG_parietal(filepath, expectedChan)
%LOADRAWEEG_PARIETAL  Load the first 2-D numeric variable with expectedChan rows.
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
