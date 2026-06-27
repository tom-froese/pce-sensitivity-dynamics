%% plotFigure3_Perceptual.m
% =========================================================================
% PNAS Figure 3: Perceptual Evidence — PAS Ratings
% =========================================================================
%
% Single-panel figure:
%
%   PAS ratings over trial time (5-s disjoint bins).  All four PAS levels
%   with SEM error bars.  Vertical shading marks the 95% CI on the
%   within-participant logistic crossover (PAS 4 vs PAS 3).
%
%   (The sensitivity schematic formerly in Panel A is now in Figure 2.)
%
% INPUT:
%   ../../results/Figure3_pas_proportions.csv
%   ../../results/Figure3_crossover_stats.csv  (from computePASCrossover.m)
%
% OUTPUT:
%   ../../results/Figure3_Perceptual.pdf  (vector PDF)
%   ../../results/Figure3_Perceptual.png  (300 dpi raster for inspection)
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   May 2026
% =========================================================================

%% ========================================================================
%  1. LOAD DATA
%  ========================================================================

scriptDir = fileparts(mfilename('fullpath'));
ROOT      = fullfile(scriptDir, '..', '..');
resDir    = fullfile(ROOT, 'results');

unsmFile  = fullfile(resDir, 'Figure3_pas_proportions.csv');
xoverFile = fullfile(resDir, 'Figure3_crossover_stats.csv');

for f = {unsmFile, xoverFile}
    if ~isfile(f{1})
        error('Missing: %s\nRun computePASCrossover.m first.', f{1});
    end
end

UNSM  = readtable(unsmFile);
XOVER = readtable(xoverFile);

T          = 60;
TAU_BOOTUP = 4;   % locked boot-up commitment (s); the shared 1/e structure is built on this
tcross    = XOVER.trial_crossover_s;
ci_lo     = XOVER.within_part_ci_lo_s;
ci_hi     = XOVER.within_part_ci_hi_s;

fprintf('Logistic crossover  = %.1f s\n', tcross);
fprintf('95%% CI: [%.1f, %.1f] s\n\n', ci_lo, ci_hi);

%% ========================================================================
%  2. AGGREGATE TO 5-s DISJOINT BINS
%  ========================================================================

binWidth = 5;
UNSM.bin5 = floor((UNSM.bin_center - 0.5) / binWidth) * binWidth + binWidth/2;

binCenters5 = unique(UNSM.bin5)';
nBin5       = length(binCenters5);

grandProp5 = nan(4, nBin5);
grandSEM5  = nan(4, nBin5);
grandN5    = nan(4, nBin5);

for lv = 1:4
    mask_lv = UNSM.pas_level == lv;
    for bi = 1:nBin5
        mask_bin = mask_lv & UNSM.bin5 == binCenters5(bi);
        if sum(mask_bin) > 0
            grandProp5(lv, bi) = mean(UNSM.grand_prop(mask_bin)) * 100;
            grandSEM5(lv, bi)  = mean(UNSM.grand_sem(mask_bin)) * 100;
            grandN5(lv, bi)    = max(UNSM.n_participants(mask_bin));
        end
    end
end

%% ========================================================================
%  3. STYLING
%  ========================================================================

col_sens  = [0.85 0.20 0.10];

col_pas = [0.65 0.65 0.65;       % PAS 1
           0.93 0.69 0.13;       % PAS 2
           0.85 0.33 0.10;       % PAS 3
           0.00 0.00 0.00];      % PAS 4
lw_pas    = [1.0 1.2 1.5 2.0];
alpha_pas = [0.6 0.8 1.0 1.0];
lw_err    = [0.6 0.6 0.8 0.8];   % thinner error bars

font_main  = 11;
font_label = 12;
font_title = 13;

%% ========================================================================
%  4. CREATE FIGURE
%  ========================================================================

fig = figure('Units', 'inches', 'Position', [0.5 0.5 8.0 3.8], 'Color', 'w', ...
    'PaperUnits', 'inches', 'PaperSize', [8.0 3.8], 'PaperPosition', [0 0 8.0 3.8]);

% Explicit position: [left, bottom, width, height]
panelLeft   = 0.10;
panelWidth  = 0.76;
panelHeight = 0.72;

%% ---- PAS ratings (5-s bins) with crossover CI -------------------------

ax_b = axes('Position', [panelLeft, 0.16, panelWidth, panelHeight]);
hold on;

% Crossover CI shading (behind everything)
fill([ci_lo ci_hi ci_hi ci_lo], [-8 -8 62 62], ...
    col_sens, 'FaceAlpha', 0.07, 'EdgeColor', 'none', ...
    'DisplayName', sprintf('95%% CI [%.0f, %.0f] s', ci_lo, ci_hi));

% Logistic crossover
xline(tcross, '-', 'Color', [0.33 0.33 0.33], 'LineWidth', 1.5, 'Alpha', 0.7, ...
    'HandleVisibility', 'off');
plot(nan, nan, '-', 'Color', [0.33 0.33 0.33], 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Crossover = %.0f s (logistic)', tcross));

% Predicted sensitivity peak: the 1/e landmark of the SHARED dynamical
% structure — the locked boot-up (TAU_BOOTUP), then 1/e decay over the
% remaining window. This is the framework's parameter-free prediction
% (lambda = e), the common anchor all four channels are tested against; it is
% NOT derived from any single channel (e.g. the click data).
tpeak = TAU_BOOTUP + (T - TAU_BOOTUP) / exp(1);   % = 4 + 56/e ~ 24.6 s -> 25 s
xline(tpeak, '--', 'Color', col_sens, 'LineWidth', 1.5, 'Alpha', 0.85, ...
    'HandleVisibility', 'off');
plot(nan, nan, '--', 'Color', col_sens, 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Sensitivity peak = %.0f s', tpeak));

% PAS traces (4 → 1 for z-ordering)
% Error bars drawn separately so they can be thinner than data lines.
pasHandles = gobjects(4,1);
errColor   = @(c) c + 0.5 * (1 - c);  % lighten colour for error bars

for lv = [4 3 2 1]
    valid = ~isnan(grandProp5(lv,:));
    xb = binCenters5(valid);
    yb = grandProp5(lv, valid);
    eb = grandSEM5(lv, valid);
    eb(isnan(eb)) = 0;

    % Thin, faint error bars (no legend entry)
    errorbar(xb, yb, eb, 'LineStyle', 'none', ...
        'Color', errColor(col_pas(lv,:)), 'LineWidth', lw_err(lv), ...
        'CapSize', 4, 'HandleVisibility', 'off');

    % Data line + markers on top
    pasHandles(lv) = plot(xb, yb, '-o', ...
        'Color', col_pas(lv,:), 'LineWidth', lw_pas(lv), ...
        'MarkerSize', 4.5, 'MarkerFaceColor', col_pas(lv,:), ...
        'DisplayName', sprintf('PAS %d', lv));

    if alpha_pas(lv) < 1
        pasHandles(lv).Color(4) = alpha_pas(lv);
        pasHandles(lv).MarkerFaceColor = col_pas(lv,:) + ...
            (1 - alpha_pas(lv)) * (1 - col_pas(lv,:));
    end
end

% Sample size annotations (unique participants per 5-s bin)
for bi = 1:nBin5
    text(binCenters5(bi), -5, sprintf('N=%d', round(grandN5(4,bi))), ...
        'HorizontalAlignment', 'center', 'FontSize', 5.5, 'Color', [0.55 0.55 0.55]);
end

hold off;
%xlim([-1 T+2]); ylim([-8 62]);
xlim([0 60]); ylim([-8 62]);
set(gca, 'XTick', 0:10:60);
xlabel('Trial time (s)', 'FontSize', font_label);
ylabel('Participant proportion (%)', 'FontSize', font_label);
title('PAS ratings over trial time (5-s disjoint bins)', ...
    'FontSize', font_title, 'FontWeight', 'bold');

% Legend: CI + lines first, then PAS 4 down to 1
lgd = legend('Location', 'southwest', 'NumColumns', 1, 'FontSize', 7.5, ...
    'Box', 'on');
set(lgd, 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);
lgd.Position(2) = lgd.Position(2) + 0.04;
set(gca, 'TickDir', 'out', 'FontSize', font_main, 'Box', 'off');

%% ========================================================================
%  5. SAVE
%  ========================================================================

outPDF = fullfile(resDir, 'Figure3_Perceptual.pdf');
outPNG = fullfile(resDir, 'Figure3_Perceptual.png');
exportgraphics(fig, outPDF, 'ContentType', 'vector', 'BackgroundColor', 'w');
exportgraphics(fig, outPNG, 'Resolution', 300, 'BackgroundColor', 'w');
fprintf('Saved: %s (vector PDF)\n', outPDF);
fprintf('Saved: %s (300 dpi PNG)\n', outPNG);

fprintf('\nDone.\n');
