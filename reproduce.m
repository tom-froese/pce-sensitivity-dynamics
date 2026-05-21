function reproduce(fromRaw)
%REPRODUCE  Master script to reproduce all results from Froese (2026).
%
%   reproduce            — from preprocessed data (default)
%   reproduce('raw')     — from raw data (requires data/raw/)
%
% Prerequisites:
%   MATLAB R2022a or later with Statistics and Optimization toolboxes
%   EEGLAB (https://sccn.ucsd.edu/eeglab/)
%
% Data:
%   Preprocessed: https://doi.org/10.5281/zenodo.19425014
%   Raw:          https://osf.io/47n3p

close all; clc;

if nargin < 1, fromRaw = ''; end
fromRaw = strcmpi(fromRaw, 'raw');

ROOT = fileparts(mfilename('fullpath'));
prepDir = fullfile(ROOT, 'code', 'preprocessing');
anaDir  = fullfile(ROOT, 'code', 'analysis');

if fromRaw
    mode = 'from raw data';
else
    mode = 'from preprocessed data';
end

fprintf('=============================================\n');
fprintf('  Reproducing: Froese (2026) PNAS Brief Report\n');
fprintf('  Mode: %s\n', mode);
fprintf('=============================================\n\n');

%% ====================================================================
%  STEP 1: Preprocessing
%  ====================================================================

if fromRaw
    fprintf('--- Step 1: Preprocessing raw data ---\n');

    if ~isfolder(fullfile(ROOT, 'data', 'raw'))
        error('data/raw/ not found. Download from https://osf.io/47n3p');
    end

    runStep(prepDir, 'preprocessClicks',            '1a', 'Preprocessing click response times');
    runStep(prepDir, 'preprocessHaptics',            '1b', 'Preprocessing haptic feedback');
    runStep(prepDir, 'preprocessEDA',                '1c', 'Preprocessing EDA');
    runStep(prepDir, 'preprocessGSP',                '1d', 'Preprocessing global scalp potential');
    runStep(prepDir, 'computeGSPStats',              '1e', 'Computing GSP statistics');
    runStep(prepDir, 'extractAllChannels',            '1f', 'Extracting all 64 EEG channels');
    runStep(prepDir, 'extractParietalHemispheres',    '1g', 'Extracting parietal hemisphere data');
    runStep(prepDir, 'computePerChannelFits',         '1h', 'Computing per-channel sensitivity fits');
    runStep(prepDir, 'computePASProportions',         '1i', 'Computing PAS proportions');

    fprintf('\n  Step 1 complete.\n\n');
else
    fprintf('--- Step 1: Skipped (using preprocessed data) ---\n');
    fprintf('  Checking preprocessed data exists...\n');

    required = {
        'data/preprocessed/ClickTimes/ClickResponseTimes.csv'
        'data/preprocessed/Haptics/HapticFeedback.csv.gz'
        'data/preprocessed/EDA/EDA_Task_Preprocessed.csv'
        'data/preprocessed/EEG/globalScalpPotential_data.mat'
        'data/preprocessed/EEG/globalScalpPotential_stats.mat'
        'data/preprocessed/EEG/allchannel_data.mat'
        'data/preprocessed/EEG/parietal_hemisphere_data.mat'
        'results/pas_unsmoothed_proportions.csv'
    };

    missing = false;
    for i = 1:numel(required)
        fp = fullfile(ROOT, required{i});
        if ~isfile(fp)
            fprintf('  MISSING: %s\n', required{i});
            missing = true;
        end
    end

    if missing
        fprintf('\n  Some preprocessed files are missing. Either:\n');
        fprintf('    1. Download from https://doi.org/10.5281/zenodo.19425014\n');
        fprintf('    2. Run reproduce(''raw'') (requires data/raw/)\n');
        error('Missing preprocessed data.');
    end
    fprintf('  All preprocessed data present.\n\n');
end

%% ====================================================================
%  STEP 2: Generate figures
%  ====================================================================

fprintf('--- Step 2: Derived statistics ---\n');

runStep(prepDir, 'computePerChannelFits', '2a', 'Per-channel sensitivity fits');

fprintf('\n--- Step 3: Generating manuscript figures ---\n');

runStep(anaDir, 'plotFigure1_Behavioral', '3a', 'Figure 1: Behavioral and bodily evidence');
runStep(anaDir, 'plotFigure2_Neural',     '3b', 'Figure 2: Neural evidence');
runStep(anaDir, 'computePASCrossover',    '3c', 'PAS crossover statistics');
runStep(anaDir, 'plotFigure3_Perceptual', '3d', 'Figure 3: Perceptual evidence');

fprintf('\n=============================================\n');
fprintf('  Done. Figures saved to results/\n');
fprintf('=============================================\n');

for pat = {'Figure1_*.pdf', 'Figure1_*.png', ...
           'Figure2_*.pdf', 'Figure2_*.png', ...
           'Figure3_*.pdf', 'Figure3_*.png'}
    d = dir(fullfile(ROOT, 'results', pat{1}));
    for i = 1:numel(d)
        fprintf('  %s  (%d KB)\n', d(i).name, round(d(i).bytes/1024));
    end
end

end


function runStep(scriptDir, scriptName, stepId, description)
    fprintf('  [%s] %s...\n', stepId, description);
    t0 = tic;
    oldDir = cd(scriptDir);
    try
        runIsolated(scriptName);
    catch ME
        cd(oldDir);
        fprintf('  ERROR in step %s: %s\n', stepId, ME.message);
        rethrow(ME);
    end
    cd(oldDir);
    fprintf('  [%s] Done (%.1f s)\n\n', stepId, toc(t0));
end


function runIsolated(scriptName)
%RUNISOLATED  Execute a script in an isolated workspace.
%   Prevents clear/clear all inside the script from wiping out
%   the caller's variables (oldDir, t0, etc.).
    run(scriptName);
end
