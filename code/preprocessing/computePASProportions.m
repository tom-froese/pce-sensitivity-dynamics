%% computePASProportions.m
% =========================================================================
% Per-PAS-Level Click Proportions in Disjoint 1-Second Bins
% =========================================================================
%
% Computes the proportion of clicks at each PAS level (1-4) within
% disjoint 1-s time bins across a 60-s trial, then averages across
% participants (grand mean and SEM).
%
% Pipeline mirrors computePASCrossover.m:
%   - Load preprocessed click-response times
%   - Load raw per-participant PAS ratings (click_presence question)
%   - Exclude dyad 31
%   - Merge by DyadNum + ParticipantID + TrialNum
%   - For each participant, compute the proportion of clicks in each
%     disjoint 1-s bin [t, t+1) that fall at each PAS level (1-4)
%   - Grand mean and SEM across participants per bin, per level
%
% INPUT:
%   data/preprocessed/ClickTimes/ClickResponseTimes.csv
%   data/raw/Behavior/pce*/questionnaires/pair_*_PAS_confidence_absence.csv
%
% OUTPUT:
%   results/pas_unsmoothed_proportions.csv
%     columns: bin_center, pas_level, grand_prop, grand_sem, n_participants
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   May 2026
% =========================================================================

scriptDir = fileparts(mfilename('fullpath'));
ROOT      = fullfile(scriptDir, '..', '..');
dataDir   = fullfile(ROOT, 'data');
resDir    = fullfile(ROOT, 'results');

T_TRIAL   = 60;
BIN_W     = 1;                          % disjoint 1-s bins
BIN_EDGES = 0:BIN_W:(T_TRIAL + BIN_W);
BIN_CENTERS = BIN_EDGES(1:end-1) + BIN_W / 2;
nB = length(BIN_CENTERS);

%% ========================================================================
%  1. LOAD CLICK TIMES
%  ========================================================================

CT = readtable(fullfile(dataDir, 'preprocessed', 'ClickTimes', ...
    'ClickResponseTimes.csv'));
CT.DyadNum = floor(CT.DyadID / 1e6);
CT = CT(CT.DyadNum ~= 31, :);
CT = CT(CT.Clicked == 1 & ~isnan(CT.ClickTime_s), :);
CT = CT(:, {'DyadNum', 'ParticipantID', 'TrialNum', 'ClickTime_s'});

fprintf('Clicks retained: %d\n', height(CT));

%% ========================================================================
%  2. LOAD PAS FROM RAW QUESTIONNAIRES
%  ========================================================================

rawBehDir = fullfile(dataDir, 'raw', 'Behavior');
allDyads  = dir(rawBehDir);
allDyads  = allDyads([allDyads.isdir]);

pasRows = [];
for i = 1:length(allDyads)
    tok = regexp(allDyads(i).name, '^pce(\d{2})\d{6}$', 'tokens');
    if isempty(tok), continue; end
    dN = str2double(tok{1}{1});
    if dN == 31, continue; end
    for p = [1 2]
        fn = sprintf('pair_%02d_P%d_PAS_confidence_absence.csv', dN, p);
        fp = fullfile(rawBehDir, allDyads(i).name, 'questionnaires', fn);
        if ~exist(fp, 'file'), continue; end
        tbl = readtable(fp, 'TextType', 'string');
        pasMask = tbl.question_id == "click_presence";
        if sum(pasMask) == 0, continue; end
        trials_0 = tbl.trial_id(pasMask);
        pasVals  = tbl.answer(pasMask);
        for ti = 1:length(trials_0)
            pasRows = [pasRows; dN, p, trials_0(ti)+1, pasVals(ti)]; %#ok<AGROW>
        end
    end
end
PAS = array2table(pasRows, 'VariableNames', ...
    {'DyadNum', 'ParticipantID', 'TrialNum', 'PAS'});

fprintf('PAS rows: %d\n', height(PAS));

%% ========================================================================
%  3. MERGE
%  ========================================================================

CT.Key  = CT.DyadNum * 1000 + CT.ParticipantID * 100 + CT.TrialNum;
PAS.Key = PAS.DyadNum * 1000 + PAS.ParticipantID * 100 + PAS.TrialNum;
[~, iCT, iPAS] = intersect(CT.Key, PAS.Key);

merged = table();
merged.DyadNum       = CT.DyadNum(iCT);
merged.ParticipantID = CT.ParticipantID(iCT);
merged.TrialNum      = CT.TrialNum(iCT);
merged.ClickTime     = CT.ClickTime_s(iCT);
merged.PAS           = PAS.PAS(iPAS);
merged.PartID        = merged.DyadNum * 10 + merged.ParticipantID;

fprintf('Merged click-PAS trials: %d\n', height(merged));
fprintf('Unique participants: %d\n', length(unique(merged.PartID)));

%% ========================================================================
%  4. PER-PARTICIPANT PROPORTIONS IN DISJOINT 1-S BINS
%  ========================================================================

uParts = unique(merged.PartID);
nP     = length(uParts);

% prop(lv, pi, bi) = proportion of clicks at PAS level lv for participant
% pi in bin bi; NaN if participant had 0 clicks in that bin.
prop = nan(4, nP, nB);

for pi = 1:nP
    mask = merged.PartID == uParts(pi);
    t    = merged.ClickTime(mask);
    pv   = merged.PAS(mask);
    binIdx = floor(t / BIN_W) + 1;          % 1-based bin index

    for b = 1:nB
        inBin = binIdx == b;
        nHere = sum(inBin);
        if nHere < 1, continue; end         % leave NaN
        for lv = 1:4
            prop(lv, pi, b) = mean(pv(inBin) == lv);
        end
    end
end

% Grand mean and SEM across participants
gm     = squeeze(nanmean(prop, 2));                     % (4 x nB)
gsd    = squeeze(nanstd(prop, 0, 2));                    % ddof=1 default
nInBin = squeeze(sum(~isnan(prop), 2));                  % (4 x nB)
gsem   = nan(4, nB);
gsem(nInBin > 1) = gsd(nInBin > 1) ./ sqrt(nInBin(nInBin > 1));

%% ========================================================================
%  5. WRITE RESULTS
%  ========================================================================

outRows = [];
for lv = 1:4
    for bi = 1:nB
        outRows = [outRows; ...
            BIN_CENTERS(bi), lv, gm(lv, bi), gsem(lv, bi), nInBin(lv, bi)]; %#ok<AGROW>
    end
end

outTbl = array2table(outRows, 'VariableNames', ...
    {'bin_center', 'pas_level', 'grand_prop', 'grand_sem', 'n_participants'});

outFile = fullfile(resDir, 'pas_unsmoothed_proportions.csv');
writetable(outTbl, outFile);
fprintf('\nWrote: %s\n', outFile);
disp(outTbl(1:12, :));
