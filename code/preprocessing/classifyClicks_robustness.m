%% classifyClicks_robustness.m
% =========================================================================
% Robustness Check for Click-Target Classification
% =========================================================================
%
% PURPOSE:
%   Reruns the Froese et al. (2014) click-target classification under
%   multiple variants to confirm that the primary conclusion (wrong clicks
%   are NOT left-shifted relative to correct clicks) is stable with respect
%   to specific algorithmic choices.
%
%   Variants tested:
%     V0 (primary):   lookback = 1.0 s, range = +/- 70  -> avatar > shadow > static > none
%     V1:             lookback = 0.5 s, range = +/- 70
%     V2:             lookback = 2.0 s, range = +/- 70
%     V3:             last-contact-in-window rule
%                       For each click, examine motor activity (haptic contact)
%                       during the 2-second window preceding the click.
%                       ClickTarget = identity of the most recent object the
%                       clicker made haptic contact with (using +/- 20 unit
%                       collision range, priority avatar > shadow > static).
%                       If no contact in window, fall back to V0.
%
%   For each variant, reports:
%     - Target distribution (avatar / shadow / static / none)
%     - Overall accuracy (fraction classified as avatar)
%     - Shadow-vs-avatar KS statistic and p-value
%     - Median and peak-of-KDE difference (shadow - avatar)
%
% OUTPUT:
%   ../../results/verification/classification_robustness.csv
%   Plus console summary.
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   April 2026
% =========================================================================

clear; clc;

%% Parameters
L             = 600;                       % Line length (wrap)
LOOKBACKS     = [0.5, 1.0, 2.0];           % seconds
CLASS_RANGE   = 70;                        % primary range
HAPTIC_RANGE  = 20;                        % haptic contact radius
HAPTIC_WINDOW = 2.0;                       % last-contact window (seconds)
T_trial       = 60;

scriptDir   = fileparts(mfilename('fullpath'));
rawDir      = fullfile(scriptDir, '..', '..', 'data', 'raw', 'Behavior');
clickDir    = fullfile(scriptDir, '..', '..', 'data', 'preprocessed', 'ClickTimes');
clickCsv    = fullfile(clickDir, 'ClickResponseTimes.csv');
resultsDir  = fullfile(scriptDir, '..', '..', 'results', 'verification');
if ~isfolder(resultsDir); mkdir(resultsDir); end

fprintf('==========================================================\n');
fprintf('  Click-Target Classification - Robustness Check\n');
fprintf('==========================================================\n\n');

%% Load click response times
fprintf('Loading %s ...\n', clickCsv);
clicksTbl = readtable(clickCsv);
nRows = height(clicksTbl);

%% Build raw-trial path lookup
folders = dir(fullfile(rawDir, 'pce*'));
folders = folders([folders.isdir]);
dyadPathMap = containers.Map('KeyType','double','ValueType','char');
for f = 1:length(folders)
    tok = regexp(folders(f).name, 'pce(\d+)', 'tokens');
    if ~isempty(tok)
        dyad_id = str2double(tok{1}{1});
        dyadPathMap(dyad_id) = fullfile(rawDir, folders(f).name, 'trials');
    end
end

%% Read each trial CSV once into a cache keyed by (dyad,trial)
fprintf('Caching raw trial CSVs (this may take a minute) ...\n');
cache = containers.Map('KeyType','char','ValueType','any');
uniqDyadTrial = unique(clicksTbl(:,{'DyadID','TrialNum'}), 'rows');
for i = 1:height(uniqDyadTrial)
    dyad_id = uniqDyadTrial.DyadID(i);
    trial_num = uniqDyadTrial.TrialNum(i);
    if ~isKey(dyadPathMap, dyad_id); continue; end
    trialsDir = dyadPathMap(dyad_id);
    candidate = {sprintf('pair_%02d_trial_%d.csv', dyad_id, trial_num), ...
                 sprintf('pair_%d_trial_%d.csv',  dyad_id, trial_num)};
    csv_path = '';
    for c = 1:length(candidate)
        p = fullfile(trialsDir, candidate{c});
        if isfile(p); csv_path = p; break; end
    end
    if isempty(csv_path)
        hits = dir(fullfile(trialsDir, sprintf('*trial_%d.csv', trial_num)));
        if ~isempty(hits); csv_path = fullfile(trialsDir, hits(1).name); end
    end
    if ~isempty(csv_path) && isfile(csv_path)
        try
            cache(sprintf('%d_%d', dyad_id, trial_num)) = readtable(csv_path);
        catch
        end
    end
end
fprintf('  Cached %d trial CSVs.\n\n', cache.Count);

%% Apply variants
variantNames = {'V0: lookback=1.0s, range=70 (primary)', ...
                'V1: lookback=0.5s, range=70', ...
                'V2: lookback=2.0s, range=70', ...
                'V3: last-contact in 2.0s window (range=20)'};
variantTgt   = cell(length(variantNames), 1);
for v = 1:length(variantNames); variantTgt{v} = repmat("none", nRows, 1); end

fprintf('Running 4 classification variants ...\n');

for i = 1:nRows
    dyad_id   = clicksTbl.DyadID(i);
    part_id   = clicksTbl.ParticipantID(i);
    trial_num = clicksTbl.TrialNum(i);
    clicked   = clicksTbl.Clicked(i);
    if clicked ~= 1; continue; end

    key = sprintf('%d_%d', dyad_id, trial_num);
    if ~isKey(cache, key); continue; end
    T = cache(key);

    % V0, V1, V2: priority-rule classifier with different lookbacks
    for v = 1:3
        lbl = classify_priority(T, part_id, LOOKBACKS(v), CLASS_RANGE, L);
        variantTgt{v}(i) = lbl;
    end

    % V3: last-contact-in-window rule (2 s window, 20 unit range)
    lbl3 = classify_last_contact(T, part_id, HAPTIC_WINDOW, HAPTIC_RANGE, L);
    if lbl3 == "none"
        % fall back to V0
        lbl3 = classify_priority(T, part_id, 1.0, CLASS_RANGE, L);
    end
    variantTgt{4}(i) = lbl3;
end

%% Summarise and compare each variant
clickMask = clicksTbl.Clicked == 1 & ~isnan(clicksTbl.ClickTime_s) & ...
            clicksTbl.ClickTime_s >= 0 & clicksTbl.ClickTime_s < T_trial;
clickTimes = clicksTbl.ClickTime_s(clickMask);

rows_out = {};
fprintf('\n==========================================================\n');
fprintf('Variant summaries:\n');
fprintf('==========================================================\n');
for v = 1:length(variantNames)
    tgt = variantTgt{v}(clickMask);
    nTotal = length(tgt);
    nAv = sum(tgt == "avatar");
    nSh = sum(tgt == "shadow");
    nSt = sum(tgt == "static");
    nNo = sum(tgt == "none");
    accuracy = nAv / nTotal;

    rt_av = clickTimes(tgt == "avatar");
    rt_sh = clickTimes(tgt == "shadow");

    if length(rt_av) > 3 && length(rt_sh) > 3
        [~, p_ks, ks_stat] = kstest2(rt_av, rt_sh);
        med_diff = median(rt_sh) - median(rt_av);

        % KDE peaks
        xi = linspace(0, T_trial, 241);
        fa = ksdensity(rt_av, xi, 'Support', [-0.001 T_trial+0.001], ...
            'BoundaryCorrection', 'reflection');
        fs = ksdensity(rt_sh, xi, 'Support', [-0.001 T_trial+0.001], ...
            'BoundaryCorrection', 'reflection');
        [~, ia] = max(fa); [~, is] = max(fs);
        peak_diff = xi(is) - xi(ia);
    else
        ks_stat = NaN; p_ks = NaN; med_diff = NaN; peak_diff = NaN;
    end

    fprintf('\n%s\n', variantNames{v});
    fprintf('  Avatar: %4d (%5.1f%%)  Shadow: %4d (%5.1f%%)  Static: %4d (%5.1f%%)  None: %4d (%5.1f%%)\n', ...
        nAv, 100*nAv/nTotal, nSh, 100*nSh/nTotal, nSt, 100*nSt/nTotal, nNo, 100*nNo/nTotal);
    fprintf('  Accuracy (avatar frac): %.3f\n', accuracy);
    if ~isnan(ks_stat)
        fprintf('  Shadow vs Avatar RT:  KS D=%.3f, p=%.3f\n', ks_stat, p_ks);
        fprintf('  Median diff (shadow - avatar): %+.2f s\n', med_diff);
        fprintf('  Peak diff    (shadow - avatar): %+.2f s\n', peak_diff);
    end

    rows_out(end+1, :) = {variantNames{v}, nAv, nSh, nSt, nNo, accuracy, ...
        ks_stat, p_ks, med_diff, peak_diff}; %#ok<SAGROW>
end

%% Write CSV summary
outCsv = fullfile(resultsDir, 'classification_robustness.csv');
fid = fopen(outCsv, 'w');
fprintf(fid, 'Variant,NAvatar,NShadow,NStatic,NNone,Accuracy,KS_D,KS_p,MedianDiff_s,PeakDiff_s\n');
for r = 1:size(rows_out, 1)
    fprintf(fid, '"%s",%d,%d,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f\n', rows_out{r,:});
end
fclose(fid);
fprintf('\nSaved: %s\n', outCsv);

%% Cross-variant agreement
fprintf('\n==========================================================\n');
fprintf('Cross-variant agreement with primary (V0):\n');
fprintf('==========================================================\n');
v0_tgt = variantTgt{1};
for v = 2:length(variantNames)
    both_valid = (v0_tgt ~= "") & (variantTgt{v} ~= "") & clickMask;
    agree = sum(v0_tgt(both_valid) == variantTgt{v}(both_valid));
    total = sum(both_valid);
    fprintf('  V%d vs V0: %d / %d = %.1f%% agreement\n', ...
        v-1, agree, total, 100*agree/total);

    % Avatar-specific agreement
    v0_av = v0_tgt == "avatar" & both_valid;
    vv_av = variantTgt{v} == "avatar" & both_valid;
    bothAv = sum(v0_av & vv_av);
    eitherAv = sum(v0_av | vv_av);
    if eitherAv > 0
        fprintf('    Avatar classifications: Jaccard = %.3f\n', bothAv/eitherAv);
    end
end

fprintf('\nDone.\n');


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function lbl = classify_priority(T, part_id, lookback_s, range_u, L)
    lbl = "none";
    if part_id == 1; btn = T.button0; else; btn = T.button1; end
    ci = find(btn == 1, 1, 'first');
    if isempty(ci); return; end
    t_target = T.timestamp(ci) - lookback_s;
    ts = T.timestamp;
    if t_target < ts(1); idx = 1; else; [~, idx] = min(abs(ts - t_target)); end
    if part_id == 1
        self_p   = T.pos0(idx);
        part_p   = T.pos1(idx);
        shad_p   = mod(T.pos1(idx) + T.shadow_delta1(idx), L);
        stat_p   = T.static_object_0(idx);
    else
        self_p   = T.pos1(idx);
        part_p   = T.pos0(idx);
        shad_p   = mod(T.pos0(idx) + T.shadow_delta0(idx), L);
        stat_p   = T.static_object_1(idx);
    end
    d_av = wrap_d(self_p, part_p, L);
    d_sh = wrap_d(self_p, shad_p, L);
    d_st = wrap_d(self_p, stat_p, L);
    if d_av <= range_u;     lbl = "avatar";
    elseif d_sh <= range_u; lbl = "shadow";
    elseif d_st <= range_u; lbl = "static";
    else;                   lbl = "none";
    end
end

function lbl = classify_last_contact(T, part_id, window_s, range_u, L)
% Find the most-recent haptic contact in the window [t_click - window_s, t_click].
% Returns 'none' if no contact.
    lbl = "none";
    if part_id == 1; btn = T.button0; else; btn = T.button1; end
    ci = find(btn == 1, 1, 'first');
    if isempty(ci); return; end
    t_click = T.timestamp(ci);
    t_start = t_click - window_s;
    ts = T.timestamp;
    inWin = ts >= t_start & ts <= t_click;
    if ~any(inWin); return; end
    idxW = find(inWin);

    if part_id == 1
        self_p  = T.pos0(idxW);
        part_p  = T.pos1(idxW);
        shad_p  = mod(T.pos1(idxW) + T.shadow_delta1(idxW), L);
        stat_p  = T.static_object_0(idxW);
    else
        self_p  = T.pos1(idxW);
        part_p  = T.pos0(idxW);
        shad_p  = mod(T.pos0(idxW) + T.shadow_delta0(idxW), L);
        stat_p  = T.static_object_1(idxW);
    end

    d_av = wrap_vec(self_p, part_p, L);
    d_sh = wrap_vec(self_p, shad_p, L);
    d_st = wrap_vec(self_p, stat_p, L);

    in_av = find(d_av <= range_u, 1, 'last');
    in_sh = find(d_sh <= range_u, 1, 'last');
    in_st = find(d_st <= range_u, 1, 'last');

    ts_av = -Inf; ts_sh = -Inf; ts_st = -Inf;
    if ~isempty(in_av); ts_av = ts(idxW(in_av)); end
    if ~isempty(in_sh); ts_sh = ts(idxW(in_sh)); end
    if ~isempty(in_st); ts_st = ts(idxW(in_st)); end

    [t_max, k] = max([ts_av, ts_sh, ts_st]);
    if t_max == -Inf; return; end
    names = ["avatar","shadow","static"];
    lbl = names(k);
end

function d = wrap_d(a, b, L)
    d = mod(abs(a - b), L);
    d = min(d, L - d);
end

function d = wrap_vec(a, b, L)
    d = mod(abs(a - b), L);
    d = min(d, L - d);
end
