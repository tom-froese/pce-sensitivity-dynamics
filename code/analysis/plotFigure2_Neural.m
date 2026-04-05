%% plotFigure2_Neural.m
% =========================================================================
% PNAS Figure 2: Neural Evidence for Sensitivity Dynamics
% =========================================================================
%
% Four-panel figure combining neural evidence from the Global Scalp
% Potential (GSP):
%
%   Panel A — Grand-Average GSP with Sensitivity Fit
%     The inverted, smoothed GSP (grand median across N=62 participants)
%     conforms to the sensitivity function S(x) = A*x*exp(-e*x) + B.
%     Permutation test confirms significance (p = 0.003).
%
%   Panel B — Per-Participant R² Distribution + Spaghetti Plot
%     Individual participant smoothed GSP curves (faint) with grand
%     median and sensitivity fit overlaid. Inset: R² histogram.
%
%   Panel C — 7-ROI Time Courses
%     Regional decomposition from anterior (prefrontal) to posterior
%     (occipital), each with sensitivity fit. Shows anterior-to-posterior
%     gradient in peak timing and fit quality.
%
%   Panel D — ROI Peak Gradient
%     Peak sensitivity time vs. anterior-posterior position, showing
%     the systematic gradient.
%
% The text also reports that the sensitivity model FAILS for the rest
% condition (R² ~ 0, confirming task specificity).
%
% INPUT:
%   EEG/globalScalpPotential_data.mat   (from preprocessGSP.m)
%   EEG/globalScalpPotential_stats.mat  (from computeGSPStats.m)
%
% OUTPUT:
%   ../../results/Figure2_Neural.png  (600 dpi)
%
% DEPENDENCIES: Statistics and Machine Learning Toolbox
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   April 2026
% =========================================================================

clear; clc; close all;

%% ========================================================================
%  1. LOAD PRECOMPUTED DATA
%  ========================================================================

dataFile  = '../../data/EEG/globalScalpPotential_data.mat';
statsFile = '../../data/EEG/globalScalpPotential_stats.mat';

if ~isfile(dataFile)
    error(['Missing: %s\n' ...
           'Run preprocessGSP.m first, or download from Zenodo.\n' ...
           'See README for instructions.'], dataFile);
end
if ~isfile(statsFile)
    error(['Missing: %s\n' ...
           'Run computeGSPStats.m first.\n' ...
           'See README for instructions.'], statsFile);
end

D = load(dataFile);
S = load(statsFile);

% Unpack frequently used variables
tTaskSm = S.tTaskSm;
nPart   = S.nPart;
nROI    = S.nROI;
roi     = S.roi;

fprintf('Loaded GSP data: N=%d participants, %d ROIs\n', nPart, nROI);
fprintf('  Grand-average R² = %.3f (perm p = %.4f)\n', S.grandR2, S.pPerm);
fprintf('  Rest R² = %.3f (task-specific control)\n', S.grandRestR2);

%% ========================================================================
%  2. STYLING
%  ========================================================================

col_task  = [0.15 0.45 0.75];    % Deep blue
col_fit   = [0.85 0.20 0.10];    % Red
col_grey  = [0.50 0.50 0.50];    % Grey

roi_colors = [
    0.90 0.29 0.21;   % Prefrontal
    0.95 0.61 0.50;   % Frontal
    0.30 0.73 0.84;   % Fronto-Central
    0.00 0.63 0.53;   % Central
    0.24 0.33 0.53;   % Centro-Parietal
    0.52 0.57 0.70;   % Parietal
    0.49 0.38 0.28;   % Occipital
];

font_sz       = 9;
font_sz_label = 10;
font_sz_title = 11;
font_sz_annot = 8;
font_sz_panel = 16;

%% ========================================================================
%  3. CREATE FIGURE — 2×2 layout
%  ========================================================================

fprintf('Creating Figure 2 ...\n');

% PNAS two-column width
fig_w = 7.09;  fig_h = 6.5;
fig = figure('Units', 'inches', 'Position', [0.5 0.5 fig_w fig_h], ...
    'Color', 'w', 'PaperUnits', 'inches', ...
    'PaperSize', [fig_w fig_h], 'PaperPosition', [0 0 fig_w fig_h]);

%% ---- Panel A: Grand-average GSP + sensitivity fit ----------------------

ax_a = axes('Position', [0.07 0.57 0.40 0.36]);
hold on;

% Smoothed inverted grand median
plot(tTaskSm, S.grandMedian_inv, '-', 'Color', col_task, 'LineWidth', 1.5);

% Sensitivity fit
plot(tTaskSm, S.grandFit, '--', 'Color', col_fit, 'LineWidth', 2);

% Peak marker
[peakVal, peakIdx] = max(S.grandFit);
peakTime = tTaskSm(peakIdx);
plot(peakTime, peakVal, 'v', 'MarkerSize', 7, ...
    'MarkerFaceColor', col_fit, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
text(peakTime + 1.5, peakVal, sprintf('%.1f s', peakTime), ...
    'FontSize', font_sz_annot, 'FontWeight', 'bold', 'Color', col_fit);

yline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
hold off;

xlim([tTaskSm(1) tTaskSm(end)]);
xlabel('Time (s)', 'FontSize', font_sz_label);
ylabel([char(8722) 'Amplitude (' char(181) 'V)'], 'FontSize', font_sz_label);
title('Grand-average GSP', 'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'FontSize', font_sz, 'Box', 'on', 'TickDir', 'out');

legend({sprintf('GSP (N=%d)', nPart), ...
        sprintf('S(x),  R^2 = %.3f', S.grandR2)}, ...
    'FontSize', font_sz_annot, 'Location', 'southeast', 'Box', 'off');

% Annotation: permutation p-value
text(0.03, 0.95, sprintf('p_{perm} = %.3f', S.pPerm), ...
    'Units', 'normalized', 'FontSize', font_sz_annot, ...
    'VerticalAlignment', 'top', 'Color', col_fit);

text(-0.15, 1.06, 'A', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

%% ---- Panel B: Spaghetti + R² inset ------------------------------------

ax_b = axes('Position', [0.57 0.57 0.40 0.36]);
hold on;

% Individual participants (faint)
for pi = 1:nPart
    plot(tTaskSm, S.partSmTC(pi,:), '-', 'Color', [col_task 0.10], 'LineWidth', 0.3);
end

% ±1 SEM band
grandMed = median(S.partSmTC, 1, 'omitnan');
grandSEM = 1.4826 * mad(S.partSmTC, 1, 1) ./ sqrt(nPart);
fill([tTaskSm, fliplr(tTaskSm)], ...
     [grandMed + grandSEM, fliplr(grandMed - grandSEM)], ...
     col_task, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

% Grand median
plot(tTaskSm, S.grandMedian_inv, '-', 'Color', col_task, 'LineWidth', 2);

% Sensitivity fit
plot(tTaskSm, S.grandFit, '--', 'Color', col_fit, 'LineWidth', 2);

yline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
hold off;

xlim([tTaskSm(1) tTaskSm(end)]);
xlabel('Time (s)', 'FontSize', font_sz_label);
ylabel([char(8722) 'Amplitude (' char(181) 'V)'], 'FontSize', font_sz_label);
title('Per-participant GSP', 'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'FontSize', font_sz, 'Box', 'on', 'TickDir', 'out');

% Inset: R² histogram
ax_inset = axes('Position', [0.62 0.77 0.14 0.13]);
hold on;
histogram(S.partR2, 15, 'FaceColor', col_task, 'EdgeColor', 'w', 'FaceAlpha', 0.7);
xline(median(S.partR2), '-', 'Color', col_fit, 'LineWidth', 1.5);
xline(0, '--', 'Color', 'k', 'LineWidth', 0.5);
hold off;
xlabel('R^2', 'FontSize', 7);
ylabel('N', 'FontSize', 7);
dz_R2 = mean(S.partR2) / std(S.partR2);
title(sprintf('d_z = %.2f', dz_R2), 'FontSize', 7, 'FontWeight', 'bold');
set(gca, 'FontSize', 7, 'Box', 'on', 'TickDir', 'out');

text(-0.15, 1.06, 'B', 'Units', 'normalized', 'Parent', ax_b, ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

%% ---- Panel C: 7-ROI time courses (stacked in one axes) ----------------

ax_c = axes('Position', [0.07 0.08 0.52 0.40]);
hold on;

legendEntries = {};
legendHandles = [];
for ri = 1:nROI
    col = roi_colors(ri,:);

    % Smoothed grand median per ROI
    h = plot(tTaskSm, S.roiGrandInv(ri,:), '-', 'Color', col, 'LineWidth', 1.2);
    legendHandles(end+1) = h; %#ok<SAGROW>
    legendEntries{end+1} = sprintf('%s (R^2=%.2f)', roi.names{ri}, S.roiGrandR2(ri)); %#ok<SAGROW>

    % Sensitivity fit (dashed, same colour)
    plot(tTaskSm, S.roiGrandFit(ri,:), '--', 'Color', col, 'LineWidth', 1.0);

    % Peak marker
    [pkVal, pkIdx] = max(S.roiGrandFit(ri,:));
    plot(tTaskSm(pkIdx), pkVal, 'v', 'MarkerSize', 5, ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
end

yline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
hold off;

xlim([tTaskSm(1) tTaskSm(end)]);
xlabel('Time (s)', 'FontSize', font_sz_label);
ylabel([char(8722) 'Amplitude (' char(181) 'V)'], 'FontSize', font_sz_label);
title('ROI time courses', 'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'FontSize', font_sz, 'Box', 'on', 'TickDir', 'out');

legend(legendHandles, legendEntries, 'FontSize', 6.5, ...
    'Location', 'southeast', 'Box', 'off', 'NumColumns', 2);

text(-0.12, 1.06, 'C', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

%% ---- Panel D: ROI peak gradient (anterior → posterior) -----------------

ax_d = axes('Position', [0.68 0.08 0.28 0.40]);
hold on;

% ROI order: 1=Prefrontal ... 7=Occipital (already anterior→posterior)
roiIdx = 1:nROI;

for ri = roiIdx
    col = roi_colors(ri,:);
    plot(ri, S.roiGrandPeak(ri), 'o', 'MarkerSize', 10, ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
end

% Trend line
p = polyfit(roiIdx(:), S.roiGrandPeak(roiIdx(:)), 1);
xline_fit = linspace(0.5, nROI+0.5, 100);
plot(xline_fit, polyval(p, xline_fit), '--', 'Color', col_grey, 'LineWidth', 1.2);

% Correlation
[rho, pval] = corr(roiIdx(:), S.roiGrandPeak(roiIdx(:)), 'type', 'Spearman');

hold off;

xlim([0.5 nROI+0.5]);
set(gca, 'XTick', 1:nROI, 'XTickLabel', ...
    {'PF','F','FC','C','CP','P','O'}, 'XTickLabelRotation', 0);
xlabel('ROI (anterior \rightarrow posterior)', 'FontSize', font_sz_label);
ylabel('Peak time (s)', 'FontSize', font_sz_label);
title('Peak gradient', 'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'FontSize', font_sz, 'Box', 'on', 'TickDir', 'out');

% Annotation
text(0.05, 0.95, sprintf('\\rho = %.2f\np = %.3f', rho, pval), ...
    'Units', 'normalized', 'FontSize', font_sz_annot, ...
    'VerticalAlignment', 'top');

text(-0.20, 1.06, 'D', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

%% ========================================================================
%  4. SAVE
%  ========================================================================

outFile = '../../results/Figure2_Neural.png';
print(fig, outFile, '-dpng', '-r600');
fprintf('  Saved: %s (600 dpi)\n', outFile);

%% ========================================================================
%  5. SUMMARY
%  ========================================================================

fprintf('\n==========================================================\n');
fprintf('  FIGURE 2 SUMMARY\n');
fprintf('==========================================================\n');
fprintf('  Panel A — Grand-Average GSP\n');
fprintf('    N participants:     %d\n', nPart);
fprintf('    Sensitivity fit:    S(x) = A*x*exp(-e*x) + B\n');
fprintf('    Grand R^2:          %.4f\n', S.grandR2);
fprintf('    Permutation p:      %.4f  (%d shifts)\n', S.pPerm, S.nPerm);
fprintf('    Grand peak:         %.1f s\n', S.grandPeak);
fprintf('    Smoothing:          %.1f s moving average\n', S.smoothWin_s);
fprintf('  --\n');
fprintf('  Panel B — Per-Participant\n');
fprintf('    Mean R^2:           %.3f (SD = %.3f)\n', mean(S.partR2), std(S.partR2));
fprintf('    d_z (R^2 > 0):      %.2f\n', mean(S.partR2)/std(S.partR2));
fprintf('  --\n');
fprintf('  Panel C — ROI Time Courses\n');
for ri = 1:nROI
    fprintf('    %-18s  R^2 = %.3f,  peak = %.1f s\n', ...
        roi.names{ri}, S.roiGrandR2(ri), S.roiGrandPeak(ri));
end
fprintf('  --\n');
fprintf('  Panel D — Peak Gradient\n');
fprintf('    Spearman rho:       %.3f,  p = %.4f\n', rho, pval);
fprintf('  --\n');
fprintf('  Rest Control\n');
fprintf('    Rest R^2:           %.4f  (no sensitivity structure)\n', S.grandRestR2);
fprintf('==========================================================\n');
fprintf('Done.\n');
