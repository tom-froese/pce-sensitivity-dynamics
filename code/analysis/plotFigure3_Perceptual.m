%% plotFigure3_Perceptual.m
% =========================================================================
% PNAS Figure 3: Perceptual Evidence for Sensitivity Dynamics
% =========================================================================
%
% Three-panel figure linking sensitivity dynamics to perceptual awareness:
%
%   Panel A — Sensitivity Schematic
%     Mathematical relationship: the exponential decay P(0) = exp(-lambda*x)
%     and its sensitivity to rate perturbations |dP(0)/dlambda| = x*exp(-lambda*x),
%     which peaks at x = 1/lambda. At lambda = e, this peak occurs at
%     x = 1/e — the optimal stopping point.
%
%   Panel B — PAS Early vs. Late
%     Within-participant comparison: proportion of PAS = 4 ("clear
%     experience") for clicks before vs. after the sensitivity peak.
%     Early clicks report clearer awareness (d_z = 0.37, p = 0.009).
%
%   Panel C — PAS Crossover + FC Peak + Per-Participant Histogram
%     The PAS 4→3 crossover time (28.4 s from logistic regression)
%     nearly coincides with the fronto-central GSP sensitivity peak
%     (27.7 s). Histogram shows per-participant crossover distribution
%     with ROI peak markers.
%
% INPUT:
%   ../../data/EEG/globalScalpPotential_stats.mat
%   ../../data/PAS/gsp_sensitivity_pas.mat
%   ../../data/PAS/gsp_pas_roi_crossover.mat
%
% OUTPUT:
%   ../../results/Figure3_Perceptual.png  (600 dpi)
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   April 2026
% =========================================================================

clear; clc; close all;

%% ========================================================================
%  1. LOAD PRECOMPUTED RESULTS
%  ========================================================================

statsFile     = '../../data/EEG/globalScalpPotential_stats.mat';
sensPasFile   = '../../data/PAS/gsp_sensitivity_pas.mat';
crossoverFile = '../../data/PAS/gsp_pas_roi_crossover.mat';

for f = {statsFile, sensPasFile, crossoverFile}
    if ~isfile(f{1})
        error('Missing: %s\nSee README for instructions.', f{1});
    end
end

S  = load(statsFile);
SP = load(sensPasFile);
XO = load(crossoverFile);

lambda = exp(1);
T = 60;

fprintf('Loaded precomputed results:\n');
fprintf('  GSP: N=%d, grand R^2=%.3f, FC peak=%.1f s\n', ...
    S.nPart, S.grandR2, S.roiGrandPeak(3));
fprintf('  PAS early/late: N=%d, d_z=%.3f, p=%.4f\n', ...
    SP.nPAS, SP.dz_prop, SP.pval_prop);
fprintf('  PAS crossover: %.1f s (trial-level), %.1f s (within-part median)\n', ...
    XO.tcross_trial, XO.grandCrossover_median);

%% ========================================================================
%  2. STYLING
%  ========================================================================

col_decay = [0.50 0.50 0.50];   % Grey for P(0)
col_sens  = [0.85 0.20 0.10];   % Red for sensitivity
col_early = [0.15 0.45 0.75];   % Blue
col_late  = [0.85 0.55 0.10];   % Orange
col_fc    = [0.30 0.73 0.84];   % Fronto-central colour (matches Fig 2)

roi_colors = [
    0.90 0.29 0.21;   % Prefrontal
    0.95 0.61 0.50;   % Frontal
    0.30 0.73 0.84;   % Fronto-Central
    0.00 0.63 0.53;   % Central
    0.24 0.33 0.53;   % Centro-Parietal
    0.52 0.57 0.70;   % Parietal
    0.49 0.38 0.28;   % Occipital
];
roiShort = {'PF','F','FC','C','CP','P','O'};

font_sz       = 9;
font_sz_label = 10;
font_sz_title = 11;
font_sz_annot = 8;
font_sz_panel = 16;

%% ========================================================================
%  3. CREATE FIGURE — 1 row on top, 2 columns on bottom
%  ========================================================================

fprintf('Creating Figure 3 ...\n');

fig_w = 7.09;  fig_h = 5.8;
fig = figure('Units', 'inches', 'Position', [0.5 0.5 fig_w fig_h], ...
    'Color', 'w', 'PaperUnits', 'inches', ...
    'PaperSize', [fig_w fig_h], 'PaperPosition', [0 0 fig_w fig_h]);

%% ---- Panel A: Sensitivity schematic (full width, top) ------------------

ax_a = axes('Position', [0.08 0.64 0.88 0.30]);
hold on;

x = linspace(0, 1.5, 300);
P0   = exp(-lambda * x);
sens = x .* exp(-lambda * x);
sens_norm = sens / max(sens);

yyaxis left;
plot(x, P0, '-', 'Color', col_decay, 'LineWidth', 2.5);
ylabel('P(0) = e^{-\lambda x}', 'FontSize', font_sz_label, 'Color', col_decay);
set(gca, 'YColor', col_decay);
ylim([-0.05 1.05]);

yyaxis right;
plot(x, sens_norm, '-', 'Color', col_sens, 'LineWidth', 2.5);
ylabel('|\partial P(0) / \partial\lambda|  (normalized)', ...
    'FontSize', font_sz_label, 'Color', col_sens);
set(gca, 'YColor', col_sens);
ylim([-0.05 1.25]);

% Mark peak
xPeak = 1/lambda;
plot(xPeak, 1.0, 'v', 'MarkerSize', 9, 'MarkerFaceColor', col_sens, ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
text(xPeak + 0.03, 1.10, 'x = 1/e', ...
    'FontSize', font_sz_annot+1, 'FontWeight', 'bold', 'Color', col_sens);

% Shade early/late regions
yl = ylim;
fill([0 xPeak xPeak 0], [yl(1) yl(1) yl(2) yl(2)], ...
    col_early, 'FaceAlpha', 0.05, 'EdgeColor', 'none');
fill([xPeak x(end) x(end) xPeak], [yl(1) yl(1) yl(2) yl(2)], ...
    col_late, 'FaceAlpha', 0.05, 'EdgeColor', 'none');

text(0.12, 0.80, {'Robust', '(low sensitivity)'}, ...
    'FontSize', font_sz_annot, 'Color', col_early, 'FontAngle', 'italic', ...
    'Units', 'normalized');
text(0.65, 0.80, {'Sensitive', '(high sensitivity)'}, ...
    'FontSize', font_sz_annot, 'Color', col_late, 'FontAngle', 'italic', ...
    'Units', 'normalized');

hold off;
xlabel('Normalized time  x = (t - \tau) / (T - \tau)', 'FontSize', font_sz_label);
title('Sensitivity of P(0) to rate perturbations', ...
    'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'TickDir', 'out', 'FontSize', font_sz, 'Box', 'off');
xlim([0 1.2]);

text(-0.06, 1.08, 'A', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

%% ---- Panel B: PAS early vs late (bottom left) -------------------------

ax_b = axes('Position', [0.08 0.10 0.38 0.42]);
hold on;

barData = [mean(SP.propPAS4_early), mean(SP.propPAS4_late)];
barSEM  = [std(SP.propPAS4_early)/sqrt(SP.nPAS), ...
           std(SP.propPAS4_late)/sqrt(SP.nPAS)];

b = bar([1 2], barData, 0.55);
b.FaceColor = 'flat';
b.CData(1,:) = col_early;
b.CData(2,:) = col_late;
b.EdgeColor  = 'k';
b.LineWidth  = 0.8;

errorbar([1 2], barData, barSEM, 'k', 'LineStyle', 'none', ...
    'LineWidth', 1.5, 'CapSize', 8);

% Individual participants (jittered)
rng(42);
jE = 0.12 * (rand(SP.nPAS,1) - 0.5);
jL = 0.12 * (rand(SP.nPAS,1) - 0.5);
scatter(1 + jE, SP.propPAS4_early, 10, col_early, 'filled', ...
    'MarkerFaceAlpha', 0.25);
scatter(2 + jL, SP.propPAS4_late, 10, col_late, 'filled', ...
    'MarkerFaceAlpha', 0.25);

% Paired lines
for i = 1:SP.nPAS
    plot([1+jE(i), 2+jL(i)], [SP.propPAS4_early(i), SP.propPAS4_late(i)], ...
        '-', 'Color', [0.5 0.5 0.5 0.10], 'LineWidth', 0.3);
end

% Significance bracket
yBracket = max(barData) + max(barSEM) + 0.04;
plot([1 1 2 2], [yBracket-0.005 yBracket yBracket yBracket-0.005], ...
    'k-', 'LineWidth', 1);
if SP.pval_prop < 0.01
    sigStr = '**';
elseif SP.pval_prop < 0.05
    sigStr = '*';
else
    sigStr = 'n.s.';
end
text(1.5, yBracket + 0.005, sigStr, 'HorizontalAlignment', 'center', ...
    'FontSize', 14, 'FontWeight', 'bold');

text(1.5, yBracket + 0.06, ...
    sprintf('p = %.3f, d_z = %.2f', SP.pval_prop, SP.dz_prop), ...
    'HorizontalAlignment', 'center', 'FontSize', font_sz_annot);

hold off;
set(gca, 'XTick', [1 2], 'XTickLabel', {'Early', 'Late'});
ylabel('Proportion PAS = 4', 'FontSize', font_sz_label);
title(sprintf('PAS clarity (N = %d)', SP.nPAS), ...
    'FontSize', font_sz_title, 'FontWeight', 'bold');
ylim([0 yBracket + 0.12]);
set(gca, 'TickDir', 'out', 'FontSize', font_sz, 'Box', 'off');

text(-0.18, 1.08, 'B', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

%% ---- Panel C: Crossover histogram + ROI peaks (bottom right) ----------

ax_c = axes('Position', [0.57 0.10 0.40 0.42]);
hold on;

validXover = ~isnan(XO.partCrossover);
validCross = XO.partCrossover(validXover);

histogram(validCross, 12, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'w', ...
    'FaceAlpha', 0.6);

yl = ylim;
yl(2) = yl(2) * 1.3;
ylim(yl);

% Within-participant median crossover
plot([XO.grandCrossover_median XO.grandCrossover_median], yl, '-', ...
    'Color', 'k', 'LineWidth', 2);
text(XO.grandCrossover_median + 0.5, yl(2)*0.92, ...
    sprintf('PAS\n%.0f s', XO.grandCrossover_median), ...
    'FontSize', font_sz_annot, 'FontWeight', 'bold', 'Color', 'k');

% Fronto-central ROI peak (highlighted)
fcPeak = S.roiGrandPeak(3);
plot([fcPeak fcPeak], yl, '--', 'Color', col_fc, 'LineWidth', 2);
text(fcPeak - 0.5, yl(2)*0.75, ...
    sprintf('FC\n%.1f s', fcPeak), ...
    'FontSize', font_sz_annot, 'FontWeight', 'bold', 'Color', col_fc, ...
    'HorizontalAlignment', 'right');

% All ROI peaks as markers at top
yMarker = yl(2) * 0.98;
for r = 1:S.nROI
    plot(S.roiGrandPeak(r), yMarker, 'v', 'MarkerSize', 7, ...
        'MarkerFaceColor', roi_colors(r,:), 'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.4);
end

% Difference annotation
diffTime = abs(fcPeak - XO.grandCrossover_median);
text(0.95, 0.50, sprintf('\\Delta = %.1f s', diffTime), ...
    'Units', 'normalized', 'FontSize', font_sz_annot+1, ...
    'HorizontalAlignment', 'right', 'FontWeight', 'bold', ...
    'Color', [0.3 0.3 0.3]);

hold off;
xlim([0 T]);
xlabel('PAS 4/3 crossover time (s)', 'FontSize', font_sz_label);
ylabel('Count', 'FontSize', font_sz_label);
title(sprintf('PAS crossover (N = %d)', XO.nValid), ...
    'FontSize', font_sz_title, 'FontWeight', 'bold');
set(gca, 'TickDir', 'out', 'FontSize', font_sz, 'Box', 'off');

text(-0.15, 1.08, 'C', 'Units', 'normalized', ...
    'FontSize', font_sz_panel, 'FontWeight', 'bold');

%% ========================================================================
%  4. SAVE
%  ========================================================================

outFile = '../../results/Figure3_Perceptual.png';
print(fig, outFile, '-dpng', '-r600');
fprintf('  Saved: %s (600 dpi)\n', outFile);

%% ========================================================================
%  5. SUMMARY
%  ========================================================================

fprintf('\n==========================================================\n');
fprintf('  FIGURE 3 SUMMARY\n');
fprintf('==========================================================\n');
fprintf('  Panel A — Sensitivity Schematic\n');
fprintf('    P(0) = exp(-lambda*x), peak at x = 1/e\n');
fprintf('    lambda = e (fixed)\n');
fprintf('  --\n');
fprintf('  Panel B — PAS Early vs. Late\n');
fprintf('    N participants:     %d\n', SP.nPAS);
fprintf('    Prop(PAS=4) early:  %.3f\n', mean(SP.propPAS4_early));
fprintf('    Prop(PAS=4) late:   %.3f\n', mean(SP.propPAS4_late));
fprintf('    t(%d) = %.2f, p = %.4f, d_z = %.3f\n', ...
        SP.nPAS-1, SP.stats_prop.tstat, SP.pval_prop, SP.dz_prop);
fprintf('  --\n');
fprintf('  Panel C — PAS Crossover + FC Peak\n');
fprintf('    PAS crossover (within-part median): %.1f s\n', XO.grandCrossover_median);
fprintf('    PAS crossover (trial-level):        %.1f s\n', XO.tcross_trial);
fprintf('    FC sensitivity peak:                %.1f s\n', fcPeak);
fprintf('    Difference:                         %.1f s\n', diffTime);
fprintf('    N valid crossovers:                 %d\n', XO.nValid);
fprintf('==========================================================\n');
fprintf('Done.\n');
