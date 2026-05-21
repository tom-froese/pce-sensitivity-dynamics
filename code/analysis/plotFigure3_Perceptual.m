%% plotFigure3_Perceptual.m
% =========================================================================
% PNAS Figure 3: Perceptual Evidence — Sensitivity Dynamics and PAS
% =========================================================================
%
% Two-panel figure:
%
%   Panel A — Sensitivity of R(x) to rate perturbations (schematic)
%     R(x) = exp(-lambda*x) and dR/dlambda (normalized), trough at 1/e.
%
%   Panel B — PAS ratings over trial time (5-s disjoint bins)
%     All four PAS levels with SEM error bars.  Vertical shading marks
%     the 95% CI on the within-participant logistic crossover (PAS 4 vs
%     PAS 3).  Dashed line marks the predicted sensitivity peak t*.
%
% INPUT:
%   ../../results/pas_unsmoothed_proportions.csv
%   ../../results/pas_crossover_stats.csv  (from computePASCrossover.m)
%
% OUTPUT:
%   ../../results/Figure3_Perceptual.pdf  (vector)
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

unsmFile  = fullfile(resDir, 'pas_unsmoothed_proportions.csv');
xoverFile = fullfile(resDir, 'pas_crossover_stats.csv');

for f = {unsmFile, xoverFile}
    if ~isfile(f{1})
        error('Missing: %s\nRun computePASCrossover.m first.', f{1});
    end
end

UNSM  = readtable(unsmFile);
XOVER = readtable(xoverFile);

tau       = XOVER.tau_locked;
T         = 60;
lambda    = exp(1);
tPeak     = XOVER.tPeak;
tcross    = XOVER.trial_crossover_s;
ci_lo     = XOVER.within_part_ci_lo_s;
ci_hi     = XOVER.within_part_ci_hi_s;

fprintf('Sensitivity peak t* = %.1f s\n', tPeak);
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
            grandN5(lv, bi)    = mean(UNSM.n_participants(mask_bin));
        end
    end
end

%% ========================================================================
%  3. STYLING
%  ========================================================================

col_decay = [0.50 0.50 0.50];
col_sens  = [0.85 0.20 0.10];
col_early = [0.15 0.45 0.75];
col_late  = [0.85 0.55 0.10];

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

fig = figure('Units', 'inches', 'Position', [0.5 0.5 8.0 6.5], 'Color', 'w', ...
    'PaperUnits', 'inches', 'PaperSize', [8.0 6.5], 'PaperPosition', [0 0 8.0 6.5]);

% Explicit positions: [left, bottom, width, height] — guarantees alignment
panelLeft   = 0.10;
panelWidth  = 0.76;
panelHeight = 0.35;

%% ---- Panel A: Sensitivity schematic ------------------------------------

ax_a = axes('Position', [panelLeft, 0.60, panelWidth, panelHeight]);
hold on;

% ==================== ALIGNMENT SHIFT ====================
% Shift Panel A right so that x=0 aligns with t=3.9 s on Panel B
align_t = 3.9;          % <-- change this value if you want a different time
shift_fraction = align_t / T;   % 3.9/60 ≈ 0.065

ax_a.Position(1) = panelLeft + shift_fraction * panelWidth;  % move left edge right
ax_a.Position(3) = panelWidth * (1 - shift_fraction);        % shrink width to keep right edge aligned
% =========================================================

x = linspace(0, 1.0, 300);
Rx         = exp(-lambda * x);
dRdlam     = -x .* exp(-lambda * x);   % actual values (trough ≈ -0.135)
xTrough    = 1/lambda;
yTrough    = dRdlam(find(x >= xTrough, 1));  % -1/e^2

% Shade rejection / selection phases
yyaxis right;
ylim([-1 0]);
yl = ylim;
fill([0 xTrough xTrough 0], [yl(1) yl(1) yl(2) yl(2)], ...
    col_early, 'FaceAlpha', 0.06, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill([xTrough x(end) x(end) xTrough], [yl(1) yl(1) yl(2) yl(2)], ...
    col_late, 'FaceAlpha', 0.06, 'EdgeColor', 'none', 'HandleVisibility', 'off');

text(0.08, 0.08, 'Rejection phase', 'Units', 'normalized', ...
    'FontSize', 9, 'Color', col_early, 'FontAngle', 'italic');
text(0.58, 0.08, 'Selection phase', 'Units', 'normalized', ...
    'FontSize', 9, 'Color', col_late, 'FontAngle', 'italic');

% R(x) on left axis
yyaxis left;
plot(x, Rx, '-', 'Color', col_decay, 'LineWidth', 2.5);
ylabel('R(x) = e^{-\lambda x}', 'FontSize', font_label, 'Color', col_decay);
set(gca, 'YColor', col_decay);
ylim([-0.05 1.05]);

% Sensitivity on right axis (actual values, 0 at top, -1 at bottom)
yyaxis right;
plot(x, dRdlam, '-', 'Color', col_sens, 'LineWidth', 2.5);
ylabel('\partial R / \partial\lambda', 'FontSize', font_label, 'Color', col_sens);
set(gca, 'YColor', col_sens);
ylim([-1 0]);

% Trough marker
plot(xTrough, yTrough, '^', 'MarkerSize', 10, 'MarkerFaceColor', col_sens, ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
text(xTrough + 0.03, yTrough - 0.03, ...
    sprintf('x = 1/e,  \\partialR/\\partial\\lambda = %.3f', yTrough), ...
    'FontSize', 9, 'FontWeight', 'bold', 'Color', col_sens);

hold off;
xlim([0 1.0]);
xlabel({'Dimensionless time  x = t / T', ...
        'equivalently  \lambda x = k t,   k = \lambda / T'}, ...
        'FontSize', font_label);
title('(A)  Sensitivity of reliability R(x) to rate perturbations', ...
    'FontSize', font_title, 'FontWeight', 'bold');
set(gca, 'TickDir', 'out', 'FontSize', font_main, 'Box', 'off');

lg = legend('R(x): reliability', ...
            '\partialR/\partial\lambda: rate sensitivity', ...
            'Location', 'east');
lg.FontSize = 9; lg.Box = 'on';
lg.Color = 'w'; lg.EdgeColor = [0.7 0.7 0.7];

%% ---- Panel B: PAS ratings (5-s bins) with crossover CI ----------------

ax_b = axes('Position', [panelLeft, 0.08, panelWidth, panelHeight]);
hold on;

% Crossover CI shading (behind everything)
fill([ci_lo ci_hi ci_hi ci_lo], [-8 -8 62 62], ...
    col_sens, 'FaceAlpha', 0.07, 'EdgeColor', 'none', ...
    'DisplayName', sprintf('95%% CI [%.0f, %.0f] s', ci_lo, ci_hi));

% Predicted sensitivity peak
xline(tPeak, '--', 'Color', col_sens, 'LineWidth', 1.3, 'Alpha', 0.7, ...
    'HandleVisibility', 'off');
plot(nan, nan, '--', 'Color', col_sens, 'LineWidth', 1.3, ...
    'DisplayName', sprintf('t* = %.1f s (predicted)', tPeak));

% Logistic crossover
xline(tcross, '-', 'Color', [0.33 0.33 0.33], 'LineWidth', 1.5, 'Alpha', 0.7, ...
    'HandleVisibility', 'off');
plot(nan, nan, '-', 'Color', [0.33 0.33 0.33], 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Crossover = %.1f s (logistic)', tcross));

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

% Sample size annotations
for bi = 1:nBin5
    text(binCenters5(bi), -5, sprintf('n=%d', round(grandN5(4,bi))), ...
        'HorizontalAlignment', 'center', 'FontSize', 5.5, 'Color', [0.55 0.55 0.55]);
end

hold off;
%xlim([-1 T+2]); ylim([-8 62]);
xlim([0 60]); ylim([-8 62]);
set(gca, 'XTick', 0:10:60);
xlabel('Trial time (s)', 'FontSize', font_label);
ylabel('PAS proportion (%)', 'FontSize', font_label);
title('(B)  PAS ratings over trial time (5-s disjoint bins)', ...
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

outFile = fullfile(resDir, 'Figure3_Perceptual.pdf');
exportgraphics(fig, outFile, 'ContentType', 'vector');
fprintf('Saved: %s (vector PDF)\n', outFile);

fprintf('\nDone.\n');
