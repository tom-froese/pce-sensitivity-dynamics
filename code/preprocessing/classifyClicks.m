%% classifyClicks.m
% =========================================================================
% Click-Target Classification for Perceptual Crossing Experiment
% =========================================================================
%
% PURPOSE:
%   Augments the existing ClickResponseTimes.csv produced by
%   preprocessClicks.m with two additional columns:
%
%       ClickTarget  - {'avatar','shadow','static','none'}
%       Correct      - 1 if ClickTarget == 'avatar', else 0
%
%   Classification follows the algorithm of Froese, Iizuka & Ikegami (2014),
%   "Embodied social interaction constitutes social cognition in pairs of
%   humans: A minimalist virtual reality experiment", Sci. Rep. 4:3672.
%
% ALGORITHM (Froese et al. 2014):
%   - For each detected click at trial time t_click, examine the position
%     state at t_click - LOOKBACK (1.0 s before the click).
%   - Compute the (wrap-around) distance between the clicker's position
%     and three candidate targets:
%         (i)  partner avatar
%         (ii) partner shadow
%         (iii) own static object
%   - Classify the click by priority rule:
%         avatar  (|d_avatar| <= RANGE)
%         shadow  (|d_shadow| <= RANGE)
%         static  (|d_static| <= RANGE)
%         none    (no target within RANGE)
%   - RANGE = 70 units (so avatar and shadow cannot both be within
%     range simultaneously, since shadow separation is 150 units and
%     2 * 70 = 140 < 150).
%
% COORDINATE CONVENTIONS (verified against raw data, 2026-04-19):
%   - Line length L = 600 units (wraps).
%   - Wrap-around distance: min(|a-b| mod L, L - |a-b| mod L).
%   - Partner shadows are at (pos_partner + shadow_delta_partner) mod L.
%       * P1 sees partner shadow at (pos1 + shadow_delta1) mod L
%       * P2 sees partner shadow at (pos0 + shadow_delta0) mod L
%   - Each participant has ONE static object:
%       * P1 sees only static_object_0
%       * P2 sees only static_object_1
%   - Effective collision radius for haptic feedback is ~20 units;
%     classification range (70) is larger than the haptic range so that
%     the categorisation captures the "target being attended to" rather
%     than the moment of contact.
%
% INPUT:
%   data/preprocessed/ClickTimes/ClickResponseTimes.csv (from preprocessClicks.m)
%   data/raw/Behavior/pce*/trials/pair_XX_trial_Y.csv   (raw trial data)
%
% OUTPUT:
%   data/preprocessed/ClickTimes/ClickResponseTimes.csv (augmented with
%       ClickTarget, Correct columns)
%   data/preprocessed/ClickTimes/ClickResponseTimes.json (updated sidecar)
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   April 2026
% =========================================================================

%% Parameters
LINE_LENGTH   = 600;    % Virtual line length (wraps)
CLASS_RANGE   = 70;     % Classification range (Froese et al. 2014)
LOOKBACK_S    = 1.0;    % Pre-click look-back window (s)

fprintf('==========================================================\n');
fprintf('  Click-Target Classification Script (PCE)\n');
fprintf('  Algorithm: Froese, Iizuka & Ikegami (2014) Sci. Rep.\n');
fprintf('==========================================================\n');
fprintf('  Line length L:            %d units\n', LINE_LENGTH);
fprintf('  Classification range:     +/- %d units\n', CLASS_RANGE);
fprintf('  Look-back window:         %.2f s\n', LOOKBACK_S);
fprintf('==========================================================\n\n');

%% Resolve paths (script-relative)
scriptDir  = fileparts(mfilename('fullpath'));
rawDir     = fullfile(scriptDir, '..', '..', 'data', 'raw', 'Behavior');
clickDir   = fullfile(scriptDir, '..', '..', 'data', 'preprocessed', 'ClickTimes');
clickCsv   = fullfile(clickDir, 'ClickResponseTimes.csv');
clickJson  = fullfile(clickDir, 'ClickResponseTimes.json');

if ~isfile(clickCsv)
    error(['ClickResponseTimes.csv not found at:\n  %s\n' ...
           'Run preprocessClicks.m first.'], clickCsv);
end

if ~isfolder(rawDir)
    error(['Raw behavioural data not found at:\n  %s\n' ...
           'Download from OSF and place in data/raw/Behavior/.'], rawDir);
end

%% Load existing click response times
fprintf('Loading %s ...\n', clickCsv);
clicksTbl = readtable(clickCsv);
nRows = height(clicksTbl);
fprintf('  Loaded %d rows.\n\n', nRows);

% Preallocate new columns
clickTarget = repmat("none", nRows, 1);
isCorrect   = nan(nRows, 1);

%% Build a lookup of raw trial paths by dyad ID
folders = dir(fullfile(rawDir, 'pce*'));
folders = folders([folders.isdir]);

dyadPathMap = containers.Map('KeyType', 'double', 'ValueType', 'char');
for f = 1:length(folders)
    tok = regexp(folders(f).name, 'pce(\d+)', 'tokens');
    if ~isempty(tok)
        dyad_id = str2double(tok{1}{1});
        dyadPathMap(dyad_id) = fullfile(rawDir, folders(f).name, 'trials');
    end
end
fprintf('Found %d dyad folders.\n\n', dyadPathMap.Count);

%% Cache: re-read each raw CSV at most once per (dyad, trial)
% Keyed by 'dyad_trial' string.
cache = containers.Map('KeyType', 'char', 'ValueType', 'any');

% Counters
nClassified = 0;
nAvatar     = 0;
nShadow     = 0;
nStatic     = 0;
nNone       = 0;
nSkipped    = 0;

fprintf('Classifying clicks ...\n');

for i = 1:nRows
    dyad_id   = clicksTbl.DyadID(i);
    part_id   = clicksTbl.ParticipantID(i);
    trial_num = clicksTbl.TrialNum(i);
    clicked   = clicksTbl.Clicked(i);

    if clicked ~= 1
        % No click - ClickTarget stays "none", Correct stays NaN
        nSkipped = nSkipped + 1;
        continue;
    end

    % Locate the raw CSV for (dyad, trial). File naming pattern:
    %   pair_<dyad_zero_padded>_trial_<trial_num>.csv
    if ~isKey(dyadPathMap, dyad_id)
        warning('  Row %d: no folder for Dyad %d. Leaving as "none".', i, dyad_id);
        nSkipped = nSkipped + 1;
        continue;
    end

    cacheKey = sprintf('%d_%d', dyad_id, trial_num);

    if isKey(cache, cacheKey)
        T = cache(cacheKey);
    else
        trialsDir = dyadPathMap(dyad_id);
        candidate_names = { ...
            sprintf('pair_%02d_trial_%d.csv', dyad_id, trial_num), ...
            sprintf('pair_%d_trial_%d.csv',  dyad_id, trial_num)};
        csv_path = '';
        for c = 1:length(candidate_names)
            p = fullfile(trialsDir, candidate_names{c});
            if isfile(p)
                csv_path = p;
                break;
            end
        end
        if isempty(csv_path)
            % Fallback: glob
            hits = dir(fullfile(trialsDir, sprintf('*trial_%d.csv', trial_num)));
            if ~isempty(hits)
                csv_path = fullfile(trialsDir, hits(1).name);
            end
        end

        if isempty(csv_path) || ~isfile(csv_path)
            warning('  Row %d: trial CSV not found for Dyad %d trial %d. Leaving as "none".', ...
                i, dyad_id, trial_num);
            nSkipped = nSkipped + 1;
            continue;
        end

        try
            T = readtable(csv_path);
        catch ME
            warning('  Row %d: error reading %s (%s). Leaving as "none".', ...
                i, csv_path, ME.message);
            nSkipped = nSkipped + 1;
            continue;
        end
        cache(cacheKey) = T;
    end

    % Classify this click using the 2014 algorithm
    [label, ok] = classify_single_click(T, part_id, LOOKBACK_S, CLASS_RANGE, LINE_LENGTH);

    if ~ok
        nSkipped = nSkipped + 1;
        continue;
    end

    clickTarget(i) = label;
    isCorrect(i)   = double(label == "avatar");
    nClassified    = nClassified + 1;

    switch char(label)
        case 'avatar';  nAvatar = nAvatar + 1;
        case 'shadow';  nShadow = nShadow + 1;
        case 'static';  nStatic = nStatic + 1;
        case 'none';    nNone   = nNone   + 1;
    end

    if mod(i, 100) == 0
        fprintf('  processed %d / %d rows\n', i, nRows);
    end
end

fprintf('\nClassification complete.\n');
fprintf('  Classified clicks: %d\n', nClassified);
fprintf('  Skipped rows:      %d\n', nSkipped);

%% Summary
if nClassified > 0
    fprintf('\n  Target breakdown (of classified clicks):\n');
    fprintf('    avatar  (correct): %4d (%.1f%%)\n', nAvatar, 100 * nAvatar / nClassified);
    fprintf('    shadow:            %4d (%.1f%%)\n', nShadow, 100 * nShadow / nClassified);
    fprintf('    static:            %4d (%.1f%%)\n', nStatic, 100 * nStatic / nClassified);
    fprintf('    none:              %4d (%.1f%%)\n', nNone,   100 * nNone   / nClassified);
end

%% Attach new columns and write augmented CSV
clicksTbl.ClickTarget = clickTarget;
clicksTbl.Correct     = isCorrect;

% Ensure column order: DyadID, ParticipantID, TrialNum, ClickTime_s, Clicked, ClickTarget, Correct
newOrder = {'DyadID','ParticipantID','TrialNum','ClickTime_s','Clicked','ClickTarget','Correct'};
clicksTbl = clicksTbl(:, newOrder);

%% Export CSV (overwrites input with augmented version)
fprintf('\nWriting augmented %s ...\n', clickCsv);

fid = fopen(clickCsv, 'w');
fprintf(fid, 'DyadID,ParticipantID,TrialNum,ClickTime_s,Clicked,ClickTarget,Correct\n');

for i = 1:nRows
    dyad_id   = clicksTbl.DyadID(i);
    part_id   = clicksTbl.ParticipantID(i);
    trial_num = clicksTbl.TrialNum(i);
    click_t   = clicksTbl.ClickTime_s(i);
    clicked   = clicksTbl.Clicked(i);
    tgt       = char(clicksTbl.ClickTarget(i));
    corr      = clicksTbl.Correct(i);

    if isnan(click_t)
        ct_str = 'NaN';
    else
        ct_str = sprintf('%.6f', click_t);
    end

    if isnan(corr)
        corr_str = 'NaN';
    else
        corr_str = sprintf('%d', corr);
    end

    fprintf(fid, '%d,%d,%d,%s,%d,%s,%s\n', ...
        dyad_id, part_id, trial_num, ct_str, clicked, tgt, corr_str);
end
fclose(fid);

csv_info = dir(clickCsv);
fprintf('  -> %s (%.1f KB)\n', csv_info.name, csv_info.bytes / 1e3);

%% Update JSON sidecar
fprintf('\nUpdating %s ...\n', clickJson);

timestamp_str = datestr(now, 'yyyy-mm-ddTHH:MM:SS');

meta = struct();

meta.Name = 'Perceptual Crossing Experiment - Click Response Times (with Target Classification)';
meta.Description = ['First button-press response times extracted from ' ...
    'perceptual crossing experiment trials, augmented with click-target ' ...
    'classification (avatar / shadow / static / none) following the ' ...
    'algorithm of Froese, Iizuka & Ikegami (2014) Sci. Rep. 4:3672. ' ...
    'Each row represents one participant''s response in one trial. ' ...
    'Trials last 60 seconds.'];

meta.Columns = {'DyadID','ParticipantID','TrialNum','ClickTime_s','Clicked','ClickTarget','Correct'};

meta.DyadID = struct( ...
    'LongName', 'Dyad Identifier', ...
    'Description', ['Numeric identifier for the experimental dyad. ' ...
        'Extracted from the pce folder name (e.g., pce01230807 -> Dyad 1).']);

meta.ParticipantID = struct( ...
    'LongName', 'Participant Identifier', ...
    'Description', ['Participant number within the dyad. ' ...
        '1 = button0 column in raw data; 2 = button1 column.']);

meta.TrialNum = struct( ...
    'LongName', 'Trial Number', ...
    'Description', 'Trial number within the session (1-18). Parsed from the CSV filename.');

meta.ClickTime_s = struct( ...
    'LongName', 'Click Response Time', ...
    'Description', ['Timestamp (in seconds from trial onset) of the ' ...
        'participant''s first button press. NaN if no click was ' ...
        'recorded during the trial.'], ...
    'Units', 's');

meta.Clicked = struct( ...
    'LongName', 'Click Detected', ...
    'Description', 'Binary flag indicating whether a button press was detected within the trial.', ...
    'Levels', struct('x0', 'No click detected (missed trial)', ...
                     'x1', 'Click detected'));

meta.ClickTarget = struct( ...
    'LongName', 'Classified Click Target', ...
    'Description', ['Categorical label for the target toward which the ' ...
        'participant was clicking, determined from the participant''s ' ...
        'position state 1.0 s before the button press. Classification ' ...
        'follows the priority rule and +/-70 unit range described in ' ...
        'Froese, Iizuka & Ikegami (2014) Sci. Rep. 4:3672. Set to "none" ' ...
        'when no click was recorded.'], ...
    'Levels', struct( ...
        'avatar', 'Partner''s avatar was within +/-70 units (correct click)', ...
        'shadow', 'Partner''s shadow was within +/-70 units (and avatar was not)', ...
        'static', 'Own static object was within +/-70 units (and neither avatar nor shadow was)', ...
        'none',   'No target was within +/-70 units, or no click was recorded'));

meta.Correct = struct( ...
    'LongName', 'Correct Click', ...
    'Description', ['Binary flag: 1 if ClickTarget == "avatar", 0 otherwise. ' ...
        'NaN for rows with no detected click.'], ...
    'Levels', struct('x0','Click classified as shadow, static, or none', ...
                     'x1','Click classified as partner avatar'));

% Source data information
meta.SourceData = struct( ...
    'Description', 'Raw CSV trial recordings from perceptual crossing experiment.', ...
    'FileFormat', '.csv', ...
    'TrialDuration', 60, ...
    'TrialDurationUnit', 's', ...
    'TrialsPerDyad', 18, ...
    'ParticipantsPerDyad', 2, ...
    'FolderStructure', 'pceXXYYMMDD/trials/pair_XX_trial_Y.csv', ...
    'RawColumns', {{'index','timestamp','static_object_0','static_object_1', ...
        'motor_0_vibrate_software','motor_1_vibrate_software', ...
        'pos0','button0','shadow_delta0','pos1','button1','shadow_delta1'}}, ...
    'ButtonMapping', struct( ...
        'button0','Participant 1 response', ...
        'button1','Participant 2 response'));

% Classification method
meta.ClassificationMethod = struct( ...
    'Reference', 'Froese T, Iizuka H, Ikegami T (2014). Embodied social interaction constitutes social cognition in pairs of humans: A minimalist virtual reality experiment. Scientific Reports 4:3672. doi:10.1038/srep03672', ...
    'LookbackWindowSeconds', LOOKBACK_S, ...
    'ClassificationRangeUnits', CLASS_RANGE, ...
    'LineLengthUnits', LINE_LENGTH, ...
    'WrapAround', true, ...
    'PriorityOrder', {{'avatar','shadow','static','none'}}, ...
    'ShadowConvention', 'pos_partner + shadow_delta_partner (mod L)', ...
    'StaticOwnership', struct( ...
        'Participant1','static_object_0', ...
        'Participant2','static_object_1'), ...
    'RangeJustification', ['The +/-70 unit range was chosen (following ' ...
        'Froese et al. 2014) so that the partner''s avatar and shadow ' ...
        'cannot be within range at the same time: shadow separation ' ...
        'is 150 units and 2 * 70 = 140 < 150.'], ...
    'Steps', {{ ...
        '1. Identify click sample: first row where buttonX == 1.', ...
        '2. Compute target time: t_click - 1.0 s.', ...
        '3. Find row with timestamp closest to target time (or first row if earlier than trial start).', ...
        '4. At that row, compute wrap-around distance from self position to: partner avatar, partner shadow, own static object.', ...
        '5. Apply priority rule: avatar if |d| <= 70; else shadow if |d| <= 70; else static if |d| <= 70; else none.', ...
        '6. Mark Correct = 1 iff ClickTarget == "avatar".'}});

% Extraction method
meta.ExtractionMethod = struct( ...
    'Description', ['Click times were extracted by preprocessClicks.m. ' ...
        'This script (classifyClicks.m) re-opens each raw trial CSV to ' ...
        'extract the position state 1.0 s before the recorded click and ' ...
        'classifies the target under the Froese et al. (2014) priority rule.'], ...
    'Pipeline', {{'preprocessClicks.m (extracts first click per participant per trial)', ...
                  'classifyClicks.m  (augments with ClickTarget and Correct)'}});

% Data summary
meta.DataSummary = struct( ...
    'NumDyads', dyadPathMap.Count, ...
    'NumRows', nRows, ...
    'NumClassifiedClicks', nClassified, ...
    'NumAvatar', nAvatar, ...
    'NumShadow', nShadow, ...
    'NumStatic', nStatic, ...
    'NumNone',   nNone, ...
    'OverallAccuracy', nAvatar / max(nClassified,1));

% File provenance
meta.GeneratedBy = struct( ...
    'Name', 'classifyClicks.m', ...
    'Description', ['MATLAB script that augments ClickResponseTimes.csv ' ...
        'with click-target classification using the Froese et al. (2014) algorithm.'], ...
    'GenerationDateTime', timestamp_str);

json_str = jsonencode(meta);
json_str = prettify_json(json_str);
fid = fopen(clickJson, 'w');
fprintf(fid, '%s', json_str);
fclose(fid);
fprintf('  -> %s\n', clickJson);

%% Final summary
fprintf('\n==========================================================\n');
fprintf('  SUMMARY\n');
fprintf('==========================================================\n');
fprintf('  Rows processed:       %d\n', nRows);
fprintf('  Clicks classified:    %d\n', nClassified);
if nClassified > 0
    fprintf('  Avatar (correct):     %d (%.1f%%)\n', nAvatar, 100 * nAvatar / nClassified);
    fprintf('  Shadow:               %d (%.1f%%)\n', nShadow, 100 * nShadow / nClassified);
    fprintf('  Static:               %d (%.1f%%)\n', nStatic, 100 * nStatic / nClassified);
    fprintf('  None:                 %d (%.1f%%)\n', nNone,   100 * nNone   / nClassified);
end
fprintf('  Output files:\n');
fprintf('    %s\n', clickCsv);
fprintf('    %s\n', clickJson);
fprintf('==========================================================\n');
fprintf('Done.\n');


%% ========================================================================
%  LOCAL FUNCTION: classify_single_click
%  ========================================================================
function [label, ok] = classify_single_click(T, part_id, lookback_s, range_u, L)
% Classify the click target for a single click row.
%
%   T          : table containing the full trial's raw recording
%   part_id    : 1 or 2 (participant ID; P1 = button0, P2 = button1)
%   lookback_s : seconds to look back from the click (e.g., 1.0)
%   range_u    : classification range in position units (e.g., 70)
%   L          : line length for wrap-around (e.g., 600)
%
% Returns:
%   label : string in {"avatar","shadow","static","none"}
%   ok    : logical, true if classification succeeded

    label = "none";
    ok = false;

    % Select the correct button column
    if part_id == 1
        btn = T.button0;
    else
        btn = T.button1;
    end

    click_idx = find(btn == 1, 1, 'first');
    if isempty(click_idx)
        return;  % should not happen when Clicked == 1, but guard anyway
    end

    t_click = T.timestamp(click_idx);
    t_target = t_click - lookback_s;

    ts = T.timestamp;

    if t_target < ts(1)
        idx_pre = 1;
    else
        [~, idx_pre] = min(abs(ts - t_target));
    end

    % Extract positions at look-back row
    if part_id == 1
        self_pos      = T.pos0(idx_pre);
        partner_pos   = T.pos1(idx_pre);
        partner_shad  = mod(T.pos1(idx_pre) + T.shadow_delta1(idx_pre), L);
        own_static    = T.static_object_0(idx_pre);
    else
        self_pos      = T.pos1(idx_pre);
        partner_pos   = T.pos0(idx_pre);
        partner_shad  = mod(T.pos0(idx_pre) + T.shadow_delta0(idx_pre), L);
        own_static    = T.static_object_1(idx_pre);
    end

    d_avatar = wrap_dist(self_pos, partner_pos,  L);
    d_shadow = wrap_dist(self_pos, partner_shad, L);
    d_static = wrap_dist(self_pos, own_static,   L);

    if d_avatar <= range_u
        label = "avatar";
    elseif d_shadow <= range_u
        label = "shadow";
    elseif d_static <= range_u
        label = "static";
    else
        label = "none";
    end
    ok = true;
end


%% ========================================================================
%  LOCAL FUNCTION: wrap_dist
%  ========================================================================
function d = wrap_dist(a, b, L)
% Shortest distance between two points on a ring of length L.
    d = mod(abs(a - b), L);
    d = min(d, L - d);
end


%% ========================================================================
%  LOCAL FUNCTION: prettify_json
%  ========================================================================
%  MATLAB's jsonencode produces a single-line string. This function adds
%  indentation and line breaks for human readability.
%  ========================================================================
function pretty = prettify_json(json_str)
    indent = 0;
    indent_str = '    ';
    pretty = '';
    in_string = false;
    i = 1;
    n = length(json_str);

    while i <= n
        ch = json_str(i);

        if ch == '"' && (i == 1 || json_str(i-1) ~= '\')
            in_string = ~in_string;
            pretty = [pretty, ch]; %#ok<AGROW>
            i = i + 1;
            continue;
        end

        if in_string
            pretty = [pretty, ch]; %#ok<AGROW>
            i = i + 1;
            continue;
        end

        switch ch
            case {'{', '['}
                pretty = [pretty, ch, newline]; %#ok<AGROW>
                indent = indent + 1;
                pretty = [pretty, repmat(indent_str, 1, indent)]; %#ok<AGROW>

            case {'}', ']'}
                pretty = [pretty, newline]; %#ok<AGROW>
                indent = indent - 1;
                pretty = [pretty, repmat(indent_str, 1, indent), ch]; %#ok<AGROW>

            case ','
                pretty = [pretty, ',', newline]; %#ok<AGROW>
                pretty = [pretty, repmat(indent_str, 1, indent)]; %#ok<AGROW>

            case ':'
                pretty = [pretty, ': ']; %#ok<AGROW>

            otherwise
                if ~isspace(ch)
                    pretty = [pretty, ch]; %#ok<AGROW>
                end
        end

        i = i + 1;
    end
end
