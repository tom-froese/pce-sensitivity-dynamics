function reproduce(fromRaw)
%REPRODUCE  Master script to reproduce all results from Froese (2026).
%
%   reproduce            — from preprocessed data (default)
%   reproduce('raw')     — from raw data (requires data/raw/)
%
% Prerequisites:
%   MATLAB R2022a or later with Statistics and Optimization toolboxes
%   Python 3.9+ with MNE-Python, matplotlib, numpy, pandas, scipy
%
% Data:
%   Preprocessed: https://doi.org/10.5281/zenodo.19425013
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
fprintf('  Reproducing: Froese (2026) PNAS Nexus Research Report\n');
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
    runStep(prepDir, 'computePerChannelFits',         '1h', 'Computing per-channel sensitivity fits');
    runStep(prepDir, 'computePASProportions',         '1i', 'Computing PAS proportions');
    runStep(anaDir,  'computePASCrossover',           '1j', 'Computing PAS crossover + click-time/PAS Spearman statistics');
    runPythonStep(fullfile(anaDir, 'computeClickPAS.py'), '1k', ...
                  'Building the per-click ClickPAS table');
    runPythonStep(fullfile(prepDir, 'preprocessEEGForExponent.py'), '1l --all', ...
                  'Preprocessing 250 Hz cleaned EEG for the aperiodic-exponent analysis');

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
        'data/preprocessed/EEG/pce01/pce01_P1_task-raw.fif'
        'data/preprocessed/ClickPAS/ClickPAS.csv'
        'results/Figure3_pas_proportions.csv'
        'results/Figure3_crossover_stats.csv'
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
        fprintf(['    1. Download preprocessed.zip (-> data/) and eeg_task_raw_fif.zip\n' ...
                 '       (-> data/preprocessed/EEG/) from https://doi.org/10.5281/zenodo.19425013\n']);
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
runPythonStep(fullfile(anaDir, 'computeAperiodicExponent.py'), '2b', ...
              'Aperiodic exponent per (participant x within-trial bin) -- FOOOF');
runPythonStep(fullfile(anaDir, 'fitExponentSensitivity.py'),    '2c', ...
              'Aperiodic exponent S(x) fit + bootstrap peak CI');

fprintf('\n--- Step 3: Generating manuscript figures ---\n');

runStep(anaDir, 'plotFigure1_Behavioral', '3a', 'Figure 1: Behavioral and bodily evidence');
runPythonStep(fullfile(anaDir, 'plotFigure2_Neural.py'), '3b', 'Figure 2: Neural evidence');

runStep(anaDir, 'plotFigure3_Perceptual', '3c', 'Figure 3: Perceptual evidence');
runPythonStep(fullfile(anaDir, 'plotFigConcept_Sensitivity.py'), '3e', ...
              'Sensitivity-framework concept figure (R(x), S(x) curves)');

fprintf('\n--- Step 4: Derived statistics beyond the figures ---\n');

runPythonStep(fullfile(anaDir, 'computeEDAFreeLambda.py'), '4a', ...
              'EDA exponential-vs-linear test (free-lambda ~ 2.62, dR2 = +0.07)');
runPythonStep(fullfile(anaDir, 'computeWithinTrialComplexity.py'), '4b', ...
              'Within-trial neural complexity (multichannel LZc slope + spectral entropy)');
runPythonStep(fullfile(anaDir, 'computeDecouplingBF.py'), '4c', ...
              'Neural<->phenomenal decoupling Bayes factors (0.17, 0.25, 0.28)');

fprintf('\n=============================================\n');
fprintf('  Done. Figures saved to results/\n');
fprintf('=============================================\n');

for pat = {'Figure1_*.pdf', 'Figure1_*.png', ...
           'Figure2_*.pdf', 'Figure2_*.png', ...
           'Figure3_*.pdf', 'Figure3_*.png', ...
           'FigConcept_*.pdf', 'FigConcept_*.png'}
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


function runPythonStep(pyScript, stepId, description)
%RUNPYTHONSTEP  Invoke one of the spoke's Python pipeline scripts.
%   stepId may include CLI args after a space, e.g. '1j --all'.
    parts = strsplit(stepId, ' ');
    label = parts{1};
    if numel(parts) > 1
        cliArgs = strjoin(parts(2:end), ' ');
    else
        cliArgs = '';
    end

    fprintf('  [%s] %s (Python)...\n', label, description);
    t0 = tic;
    cmd = sprintf('python3 "%s" %s', pyScript, cliArgs);
    [status, cmdout] = system(cmd);
    if status ~= 0
        error('Python step %s failed:\n%s', label, cmdout);
    end
    fprintf('%s', cmdout);
    fprintf('  [%s] Done (%.1f s)\n\n', label, toc(t0));
end
