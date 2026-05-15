%% sanity_check_classification.m
% =========================================================================
% Sanity-Check of Click-Target Classification Output
% =========================================================================
%
% Verifies that classifyClicks.m has produced a sensible augmented
% ClickResponseTimes.csv. Checks:
%   1. Row counts agree with source
%   2. Correct column is coherent with ClickTarget
%   3. Per-dyad accuracy distribution
%   4. Response-time summary by target
%   5. Histogram of accuracy across dyads printed to console
%
% =========================================================================

scriptDir = fileparts(mfilename('fullpath'));
clickCsv  = fullfile(scriptDir, '..', '..', 'data', 'preprocessed', ...
    'ClickTimes', 'ClickResponseTimes.csv');

T = readtable(clickCsv);

fprintf('\n===== Sanity check of classification output =====\n');
fprintf('File: %s\n', clickCsv);
fprintf('Rows: %d\n', height(T));

% Coherence check: Correct == 1 iff ClickTarget == 'avatar'
tgt = string(T.ClickTarget);
corr = T.Correct;

valid = ~isnan(corr);
sub_tgt = tgt(valid);
sub_corr = corr(valid);

nCorrOne  = sum(sub_corr == 1);
nAvatar   = sum(sub_tgt == "avatar");
nMismatch = sum((sub_corr == 1) ~= (sub_tgt == "avatar"));
fprintf('\nCoherence of Correct vs ClickTarget:\n');
fprintf('  Correct == 1 rows:      %d\n', nCorrOne);
fprintf('  ClickTarget == avatar:  %d\n', nAvatar);
fprintf('  Mismatches:             %d  (should be 0)\n', nMismatch);

% NaN handling
nClicked    = sum(T.Clicked == 1);
nNonClicked = sum(T.Clicked == 0);
nCorrNaN    = sum(isnan(T.Correct));
fprintf('\nClicked == 0 rows:        %d\n', nNonClicked);
fprintf('Correct == NaN rows:      %d  (should match unclicked)\n', nCorrNaN);

% Target breakdown
fprintf('\nTarget distribution (classified rows):\n');
for lbl = ["avatar","shadow","static","none"]
    n = sum(sub_tgt == lbl);
    fprintf('  %-8s %4d  (%.1f%% of classified)\n', char(lbl), n, 100*n/length(sub_tgt));
end

% Per-dyad accuracy
dyads = unique(T.DyadID);
accs  = nan(length(dyads), 1);
nClicks = nan(length(dyads), 1);
for k = 1:length(dyads)
    mask = T.DyadID == dyads(k) & T.Clicked == 1;
    accs(k) = mean(T.Correct(mask), 'omitnan');
    nClicks(k) = sum(mask);
end

fprintf('\nPer-dyad accuracy:\n');
fprintf('  Dyads: %d\n', length(dyads));
fprintf('  Median: %.3f | Mean: %.3f | Min: %.3f | Max: %.3f\n', ...
    median(accs,'omitnan'), mean(accs,'omitnan'), min(accs), max(accs));
fprintf('  Dyads with acc > 0.9: %d\n', sum(accs > 0.9));
fprintf('  Dyads with acc < 0.5: %d\n', sum(accs < 0.5));

% Print dyad-level table (sorted by accuracy)
[accs_sorted, order] = sort(accs, 'descend');
fprintf('\nDyad accuracy table (descending):\n');
fprintf('  Dyad  NClicks  Acc\n');
for k = 1:length(dyads)
    fprintf('   %2d     %3d     %.3f\n', dyads(order(k)), nClicks(order(k)), accs_sorted(k));
end

% Response-time summary by target
fprintf('\nClick-time (s) by target:\n');
fprintf('  Target    N     Median     Mean     Std\n');
for lbl = ["avatar","shadow","static","none"]
    mask = tgt == lbl & ~isnan(T.ClickTime_s);
    rt = T.ClickTime_s(mask);
    if isempty(rt)
        fprintf('  %-8s  %3d    n/a       n/a      n/a\n', char(lbl), 0);
    else
        fprintf('  %-8s  %3d    %5.2f     %5.2f    %5.2f\n', ...
            char(lbl), length(rt), median(rt), mean(rt), std(rt));
    end
end

fprintf('\nDone.\n');
