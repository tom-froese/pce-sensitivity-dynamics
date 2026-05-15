%% plotFigureS1_ClickCorrectness.m
% =========================================================================
% PNAS SI Figure S1 — Click-Target Classification and Correctness Analysis
% =========================================================================
%
% Compares the response-time distributions of correct (partner-avatar)
% and incorrect clicks in the perceptual crossing experiment. Addresses
% the question raised at the Yale ISEC 2026 talk: are wrong clicks
% left-shifted (classical decision-theoretic / SPRT prediction), or do
% they follow the same sensitivity dynamics as correct clicks
% (sensitivity-dynamics prediction)?
%
% Click targets are classified by the algorithm of Froese, Iizuka &
% Ikegami (2014, Sci. Rep.) as implemented in classifyClicks.m:
%     priority: avatar > shadow > static > none
%     look-back: 1.0 s,  range: +/- 70 units,  line length: 600
%
% THREE wrong-click contrasts are shown:
%   Shadow only      - decoy that moves with partner (same dynamics
%                      as avatar but at offset +/- 75 units).
%   Shadow + Static + None (all wrong) - the classical "error"
%                      contrast: any click that was not classified
%                      as the partner avatar.
%   Per-category KDEs of avatar vs shadow vs static vs none clicks.
%
% Panels (3 x 1):
%   A - Overlaid KDEs of avatar (correct) vs shadow (wrong, moving decoy)
%       vs static (wrong, stationary decoy) vs none (wrong, no target).
%   B - Primary contrast: avatar vs ALL wrong clicks (shadow+static+none)
%       with KS test and median/peak comparison.
%   C - Bootstrap peak-difference distributions for two contrasts:
%       (i) shadow - avatar, (ii) all-wrong - avatar.
%
% INPUT:
%   ../../data/preprocessed/ClickTimes/ClickResponseTimes.csv
%       (produced by preprocessClicks.m + classifyClicks.m)
%
% OUTPUT:
%   ../../results/FigureS1_ClickCorrectness.png  (600 dpi)
%
% AUTHOR: Embodied Cognitive Science Unit, OIST
% DATE:   April 2026
% =========================================================================

clear; clc; close all;

%% Parameters
clickFile = '../../data/preprocessed/ClickTimes/ClickResponseTimes.csv';
outFile   = '../../results/FigureS1_ClickCorrectness.png';

T_trial   = 60;
lambda    = exp(1);
n_kde     = 500;
xi_kde    = linspace(0, T_trial, n_kde);
nBoot     = 2000;
rngSeed   = 42;
bin_width = 2;
edges     = 0:bin_width:T_trial;
bin_ctr   = edges(1:end-1) + bin_width/2;
tau_grid  = 0:0.1:15;

if ~isfile(clickFile)
    error('Missing: %s\nRun preprocessClicks.m + classifyClicks.m first.', clickFile);
end

%% Load data
fprintf('Loading %s ...\n', clickFile);
data = readtable(clickFile);

if ~ismember('ClickTarget', data.Properties.VariableNames)
    error(['ClickTarget column not found. Run classifyClicks.m to augment ' ...
           'the preprocessed CSV.']);
end

valid = data.Clicked == 1 & ~isnan(data.ClickTime_s) & ...
        data.ClickTime_s >= 0 & data.ClickTime_s < T_trial;
sub = data(valid, :);
tgt = string(sub.ClickTarget);

rt_avatar = sub.ClickTime_s(tgt == "avatar");
rt_shadow = sub.ClickTime_s(tgt == "shadow");
rt_static = sub.ClickTime_s(tgt == "static");
rt_none   = sub.ClickTime_s(tgt == "none");
rt_wrong  = sub.ClickTime_s(tgt ~= "avatar");  % all wrong clicks (shadow+static+none)

n_av = length(rt_avatar);
n_sh = length(rt_shadow);
n_st = length(rt_static);
n_no = length(rt_none);
n_wr = length(rt_wrong);
n_total = n_av + n_sh + n_st + n_no;

fprintf('  Clicks (classified): %d\n', n_total);
fprintf('    avatar (correct):              %d (%.1f%%)\n', n_av, 100*n_av/n_total);
fprintf('    shadow  (wrong, moving decoy): %d (%.1f%%)\n', n_sh, 100*n_sh/n_total);
fprintf('    static  (wrong, stationary):   %d (%.1f%%)\n', n_st, 100*n_st/n_total);
fprintf('    none    (wrong, no target):    %d (%.1f%%)\n', n_no, 100*n_no/n_total);
fprintf('    ALL WRONG (sh+st+no):          %d (%.1f%%)\n', n_wr, 100*n_wr/n_total);

% Per-category RT summary
fprintf('\n  Click-time summary (s):\n');
fprintf('    %-16s  n    median   mean    std\n', 'Category');
fprintf('    %-16s  %3d   %5.2f    %5.2f   %5.2f\n', 'avatar', n_av, ...
    median(rt_avatar), mean(rt_avatar), std(rt_avatar));
fprintf('    %-16s  %3d   %5.2f    %5.2f   %5.2f\n', 'shadow', n_sh, ...
    median(rt_shadow), mean(rt_shadow), std(rt_shadow));
fprintf('    %-16s  %3d   %5.2f    %5.2f   %5.2f\n', 'static', n_st, ...
    median(rt_static), mean(rt_static), std(rt_static));
fprintf('    %-16s  %3d   %5.2f    %5.2f   %5.2f\n', 'none',   n_no, ...
    median(rt_none),   mean(rt_none),   std(rt_none));
fprintf('    %-16s  %3d   %5.2f    %5.2f   %5.2f\n', 'ALL WRONG', n_wr, ...
    median(rt_wrong),  mean(rt_wrong),  std(rt_wrong));

%% KDEs and sensitivity-model fits for each category
fprintf('\nFitting sensitivity models ...\n');

fk_av = ksdensity(rt_avatar, xi_kde, ...
    'Support', [-0.001 T_trial+0.001], 'BoundaryCorrection', 'reflection');
fk_sh = ksdensity(rt_shadow, xi_kde, ...
    'Support', [-0.001 T_trial+0.001], 'BoundaryCorrection', 'reflection');
fk_st = ksdensity(rt_static, xi_kde, ...
    'Support', [-0.001 T_trial+0.001], 'BoundaryCorrection', 'reflection');
fk_no = ksdensity(rt_none, xi_kde, ...
    'Support', [-0.001 T_trial+0.001], 'BoundaryCorrection', 'reflection');
fk_wr = ksdensity(rt_wrong, xi_kde, ...
    'Support', [-0.001 T_trial+0.001], 'BoundaryCorrection', 'reflection');

best_av = fit_sensitivity_grid(xi_kde, fk_av, T_trial, lambda, tau_grid);
best_sh = fit_sensitivity_grid(xi_kde, fk_sh, T_trial, lambda, tau_grid);
best_wr = fit_sensitivity_grid(xi_kde, fk_wr, T_trial, lambda, tau_grid);

fprintf('  avatar    (correct):   tau=%.1f s, R^2=%.3f, peak=%.2f s (KDE peak=%.2f)\n', ...
    best_av.offset, best_av.R2, best_av.peak_time, xi_kde(argmax(fk_av)));
fprintf('  shadow    (wrong):     tau=%.1f s, R^2=%.3f, peak=%.2f s (KDE peak=%.2f)\n', ...
    best_sh.offset, best_sh.R2, best_sh.peak_time, xi_kde(argmax(fk_sh)));
fprintf('  ALL WRONG (sh+st+no):  tau=%.1f s, R^2=%.3f, peak=%.2f s (KDE peak=%.2f)\n', ...
    best_wr.offset, best_wr.R2, best_wr.peak_time, xi_kde(argmax(fk_wr)));

%% Empirical KDE peaks
peak_av_emp = xi_kde(argmax(fk_av));
peak_sh_emp = xi_kde(argmax(fk_sh));
peak_wr_emp = xi_kde(argmax(fk_wr));
peak_diff_sh_av = peak_sh_emp - peak_av_emp;
peak_diff_wr_av = peak_wr_emp - peak_av_emp;

fprintf('\nEmpirical KDE peaks: avatar=%.2f, shadow=%.2f (diff %+.2f), all-wrong=%.2f (diff %+.2f) s\n', ...
    peak_av_emp, peak_sh_emp, peak_diff_sh_av, peak_wr_emp, peak_diff_wr_av);

%% Bootstrap peak-difference distributions (two contrasts)
fprintf('\nBootstrap (nBoot = %d, seed = %d) ...\n', nBoot, rngSeed);
rng(rngSeed);

xi_boot = linspace(0, T_trial, 241);   % 0.25 s resolution

boot_diff_sh_av = nan(nBoot, 1);
boot_diff_wr_av = nan(nBoot, 1);

for b = 1:nBoot
    s_av = rt_avatar(randi(n_av, n_av, 1));
    s_sh = rt_shadow(randi(n_sh, n_sh, 1));
    s_wr = rt_wrong (randi(n_wr, n_wr, 1));
    fa = ksdensity(s_av, xi_boot, ...
        'Support', [-0.001 T_trial+0.001], 'BoundaryCorrection', 'reflection');
    fs = ksdensity(s_sh, xi_boot, ...
        'Support', [-0.001 T_trial+0.001], 'BoundaryCorrection', 'reflection');
    fw = ksdensity(s_wr, xi_boot, ...
        'Support', [-0.001 T_trial+0.001], 'BoundaryCorrection', 'reflection');
    boot_diff_sh_av(b) = xi_boot(argmax(fs)) - xi_boot(argmax(fa));
    boot_diff_wr_av(b) = xi_boot(argmax(fw)) - xi_boot(argmax(fa));
end

ci_sh_av = quantile(boot_diff_sh_av, [0.025 0.975]);
ci_wr_av = quantile(boot_diff_wr_av, [0.025 0.975]);

p_sh_av = 2 * min(mean(boot_diff_sh_av <= 0), mean(boot_diff_sh_av >= 0));
p_wr_av = 2 * min(mean(boot_diff_wr_av <= 0), mean(boot_diff_wr_av >= 0));

fprintf('  shadow - avatar:   diff = %+.2f s,  95%% CI [%+.2f, %+.2f],  p = %.3f\n', ...
    peak_diff_sh_av, ci_sh_av(1), ci_sh_av(2), p_sh_av);
fprintf('  all-wrong - avatar: diff = %+.2f s,  95%% CI [%+.2f, %+.2f],  p = %.3f\n', ...
    peak_diff_wr_av, ci_wr_av(1), ci_wr_av(2), p_wr_av);

%% KS tests (two-sample)
[~, p_ks_sh, ks_sh] = kstest2(rt_avatar, rt_shadow);
[~, p_ks_wr, ks_wr] = kstest2(rt_avatar, rt_wrong);

fprintf('\nKS tests:\n');
fprintf('  avatar vs shadow:   D = %.3f, p = %.3f\n', ks_sh, p_ks_sh);
fprintf('  avatar vs ALL WRONG: D = %.3f, p = %.3f\n', ks_wr, p_ks_wr);

%% Permutation tests of median difference
rng(rngSeed);
nPerm = 5000;

% shadow - avatar
med_diff_sh = median(rt_shadow) - median(rt_avatar);
pool_sh = [rt_avatar; rt_shadow];
perm_sh = nan(nPerm, 1);
for p = 1:nPerm
    idx = randperm(length(pool_sh));
    perm_sh(p) = median(pool_sh(idx(n_av+1:end))) - median(pool_sh(idx(1:n_av)));
end
p_perm_sh = mean(abs(perm_sh) >= abs(med_diff_sh));

% all-wrong - avatar
rng(rngSeed);
med_diff_wr = median(rt_wrong) - median(rt_avatar);
pool_wr = [rt_avatar; rt_wrong];
perm_wr = nan(nPerm, 1);
for p = 1:nPerm
    idx = randperm(length(pool_wr));
    perm_wr(p) = median(pool_wr(idx(n_av+1:end))) - median(pool_wr(idx(1:n_av)));
end
p_perm_wr = mean(abs(perm_wr) >= abs(med_diff_wr));

fprintf('\nPermutation (median):\n');
fprintf('  shadow  - avatar:   diff = %+.2f s,  p = %.3f\n', med_diff_sh, p_perm_sh);
fprintf('  all-wrong - avatar: diff = %+.2f s,  p = %.3f\n', med_diff_wr, p_perm_wr);

%% ========================================================================
%  CREATE FIGURE — 3 panels
%  ========================================================================

fprintf('\nCreating SI Figure S1 ...\n');

col_av   = [0.20 0.40 0.73];    % Blue  (avatar = correct)
col_wr   = [0.80 0.15 0.15];    % Red   (all-wrong)
col_sh   = [0.85 0.32 0.10];    % Orange (shadow subcategory)
col_st   = [0.55 0.35 0.65];    % Purple (static subcategory)
col_no   = [0.45 0.45 0.45];    % Grey  (none subcategory)
col_fit  = [0.10 0.10 0.10];    % Black for fit line
col_grey = [0.50 0.50 0.50];

font_sz       = 9;
font_sz_label = 10;
font_sz_title = 11;
font_sz_annot = 8;
font_sz_panel = 14;

% PNAS double-column width
fig_w = 10.0;  fig_h = 3.2;
fig = figure('Units','inches','Position',[0.5 0.5 fig_w fig_h], ...
    'Color','w','PaperUnits','inches', ...
    'PaperSize',[fig_w fig_h],'PaperPosition',[0 0 fig_w fig_h]);

ml = 0.055;  mr = 0.02;  mb = 0.17;  mt = 0.10;
gap_ab = 0.065;
gap_bc = 0.075;
total_w = 1 - ml - mr - gap_ab - gap_bc;
pw_a = total_w * 0.37;
pw_b = total_w * 0.33;
pw_c = total_w * 0.30;
ph = 1 - mb - mt;

%% ---- Panel A: Avatar vs ALL WRONG (primary contrast) -------------------
ax_a = axes('Position',[ml, mb, pw_a, ph]);
hold on;

scale_av = n_av * bin_width;
scale_wr = n_wr * bin_width;

histogram(rt_avatar, edges, ...
    'FaceColor', col_av, 'FaceAlpha', 0.20, 'EdgeColor','w', 'LineWidth', 0.5, ...
    'HandleVisibility','off');
histogram(rt_wrong, edges, ...
    'FaceColor', col_wr, 'FaceAlpha', 0.20, 'EdgeColor','w', 'LineWidth', 0.5, ...
    'HandleVisibility','off');

plot(xi_kde, fk_av * scale_av, '-', 'Color', col_av, 'LineWidth', 1.8, ...
    'DisplayName', sprintf('Avatar (correct)  n=%d', n_av));
plot(xi_kde, fk_wr * scale_wr, '-', 'Color', col_wr, 'LineWidth', 1.8, ...
    'DisplayName', sprintf('All wrong  n=%d', n_wr));

% Sensitivity fits
t_mod = linspace(0, T_trial, 500);
y_av_mod = sensitivity_curve(t_mod, best_av.A, best_av.offset, best_av.T_eff, lambda) * scale_av;
y_wr_mod = sensitivity_curve(t_mod, best_wr.A, best_wr.offset, best_wr.T_eff, lambda) * scale_wr;
plot(t_mod, y_av_mod, '--', 'Color', col_av*0.5, 'LineWidth', 1.3, 'HandleVisibility','off');
plot(t_mod, y_wr_mod, '--', 'Color', col_wr*0.5, 'LineWidth', 1.3, 'HandleVisibility','off');

xline(peak_av_emp, ':', 'Color', col_av, 'LineWidth', 1, 'Alpha', 0.8, 'HandleVisibility','off');
xline(peak_wr_emp, ':', 'Color', col_wr, 'LineWidth', 1, 'Alpha', 0.8, 'HandleVisibility','off');

xlim([0 T_trial]);
yl_A = [0, max(histcounts(rt_avatar, edges)) * 1.4];
ylim(yl_A);
xlabel('Click response time (s)', 'FontSize', font_sz_label);
ylabel('Clicks per 2-s bin', 'FontSize', font_sz_label);
title('Correct vs ALL wrong clicks','FontSize',font_sz_title,'FontWeight','bold');
set(gca,'FontSize',font_sz,'Box','off','TickDir','out');

lg_a = legend('Location','northwest','FontSize',font_sz_annot,'Box','off');
lg_a.ItemTokenSize = [16, 8];

text(0.97, 0.04, ...
    {sprintf('Avatar:    \\tau=%.1f s, peak=%.1f s, R^2=%.2f', best_av.offset, best_av.peak_time, best_av.R2), ...
     sprintf('All wrong: \\tau=%.1f s, peak=%.1f s, R^2=%.2f', best_wr.offset, best_wr.peak_time, best_wr.R2), ...
     sprintf('Median \\Delta = %+.2f s (perm p=%.2f)', med_diff_wr, p_perm_wr), ...
     sprintf('KS: D=%.2f, p=%.2f', ks_wr, p_ks_wr)}, ...
    'Units','normalized','FontSize',font_sz_annot, ...
    'HorizontalAlignment','right','VerticalAlignment','bottom', ...
    'BackgroundColor','w','EdgeColor',[0.7 0.7 0.7],'Margin',3);

text(-0.10, 1.08, 'A', 'Units','normalized','FontSize',font_sz_panel,'FontWeight','bold');
hold off;

%% ---- Panel B: Per-category KDEs ----------------------------------------
ax_b = axes('Position',[ml + pw_a + gap_ab, mb, pw_b, ph]);
hold on;

% Scale to density (so 4 categories are visually comparable)
plot(xi_kde, fk_av, '-', 'Color', col_av, 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Avatar  n=%d  (med %.1f s)', n_av, median(rt_avatar)));
plot(xi_kde, fk_sh, '-', 'Color', col_sh, 'LineWidth', 1.6, ...
    'DisplayName', sprintf('Shadow  n=%d  (med %.1f s)', n_sh, median(rt_shadow)));
plot(xi_kde, fk_st, '-', 'Color', col_st, 'LineWidth', 1.6, ...
    'DisplayName', sprintf('Static  n=%d  (med %.1f s)', n_st, median(rt_static)));
plot(xi_kde, fk_no, '-', 'Color', col_no, 'LineWidth', 1.6, ...
    'DisplayName', sprintf('None  n=%d  (med %.1f s)', n_no, median(rt_none)));

% Medians as dotted vlines
xline(median(rt_avatar), ':', 'Color', col_av, 'Alpha', 0.7, 'HandleVisibility','off');
xline(median(rt_shadow), ':', 'Color', col_sh, 'Alpha', 0.7, 'HandleVisibility','off');
xline(median(rt_static), ':', 'Color', col_st, 'Alpha', 0.7, 'HandleVisibility','off');
xline(median(rt_none),   ':', 'Color', col_no, 'Alpha', 0.7, 'HandleVisibility','off');

xlim([0 T_trial]);
xlabel('Click response time (s)', 'FontSize', font_sz_label);
ylabel('KDE (density)', 'FontSize', font_sz_label);
title('Per-target breakdown','FontSize',font_sz_title,'FontWeight','bold');
set(gca,'FontSize',font_sz,'Box','off','TickDir','out');

lg_b = legend('Location','northwest','FontSize',font_sz_annot,'Box','off');
lg_b.ItemTokenSize = [18, 8];

text(-0.13, 1.08, 'B', 'Units','normalized','FontSize',font_sz_panel,'FontWeight','bold');
hold off;

%% ---- Panel C: Overlaid bootstrap peak-difference distributions ----------
ax_c = axes('Position',[ml + pw_a + gap_ab + pw_b + gap_bc, mb, pw_c, ph]);
hold on;

% Determine common x-range covering both distributions
xr = [min([ci_sh_av(1), ci_wr_av(1), peak_diff_sh_av, peak_diff_wr_av, 0]) - 2, ...
      max([ci_sh_av(2), ci_wr_av(2), peak_diff_sh_av, peak_diff_wr_av, 0]) + 2];
binE = linspace(xr(1), xr(2), 41);

histogram(boot_diff_sh_av, binE, ...
    'FaceColor', col_sh, 'FaceAlpha', 0.45, 'EdgeColor', 'none', ...
    'DisplayName', 'shadow - avatar');
histogram(boot_diff_wr_av, binE, ...
    'FaceColor', col_wr, 'FaceAlpha', 0.45, 'EdgeColor', 'none', ...
    'DisplayName', 'all-wrong - avatar');

xline(0, '-', 'Color', col_grey, 'LineWidth', 1, 'HandleVisibility','off');
xline(peak_diff_sh_av, '-', 'Color', col_sh, 'LineWidth', 1.8, 'HandleVisibility','off');
xline(peak_diff_wr_av, '-', 'Color', col_wr, 'LineWidth', 1.8, 'HandleVisibility','off');

% CI markers
xline(ci_sh_av(1), ':', 'Color', col_sh, 'Alpha', 0.7, 'HandleVisibility','off');
xline(ci_sh_av(2), ':', 'Color', col_sh, 'Alpha', 0.7, 'HandleVisibility','off');
xline(ci_wr_av(1), ':', 'Color', col_wr, 'Alpha', 0.7, 'HandleVisibility','off');
xline(ci_wr_av(2), ':', 'Color', col_wr, 'Alpha', 0.7, 'HandleVisibility','off');

xlim(xr);
yl = ylim;
fill([xr(1) 0 0 xr(1)], [0 0 yl(2) yl(2)], col_grey, ...
    'FaceAlpha', 0.07, 'EdgeColor','none','HandleVisibility','off');

xlabel('Peak difference: wrong - avatar (s)', 'FontSize', font_sz_label);
ylabel(sprintf('Bootstrap draws (of %d)', nBoot), 'FontSize', font_sz_label);
title('Peak-time bootstrap','FontSize',font_sz_title,'FontWeight','bold');
set(gca,'FontSize',font_sz,'Box','off','TickDir','out');

lg_c = legend('Location','northwest','FontSize',font_sz_annot,'Box','off');
lg_c.ItemTokenSize = [14, 8];

text(0.97, 0.97, ...
    {sprintf('\\Delta shadow = %+.2f s', peak_diff_sh_av), ...
     sprintf('  95%% CI [%+.2f, %+.2f]', ci_sh_av(1), ci_sh_av(2)), ...
     sprintf('\\Delta all-wrong = %+.2f s', peak_diff_wr_av), ...
     sprintf('  95%% CI [%+.2f, %+.2f]', ci_wr_av(1), ci_wr_av(2))}, ...
    'Units','normalized','FontSize',font_sz_annot, ...
    'HorizontalAlignment','right','VerticalAlignment','top', ...
    'BackgroundColor','w','EdgeColor',[0.7 0.7 0.7],'Margin',3);

text(0.98, 0.04, 'SPRT prediction', ...
    'Units','normalized','FontSize',font_sz_annot - 1, ...
    'HorizontalAlignment','right','Color',col_grey);

text(-0.15, 1.08, 'C', 'Units','normalized','FontSize',font_sz_panel,'FontWeight','bold');
hold off;

%% ========================================================================
%  SAVE
%  ========================================================================

if ~isfolder('../../results'); mkdir('../../results'); end
exportgraphics(fig, outFile, 'Resolution', 600);
fprintf('\n  Saved: %s (600 dpi)\n', outFile);

%% ========================================================================
%  SUMMARY
%  ========================================================================
fprintf('\n==========================================================\n');
fprintf('  SI FIGURE S1 SUMMARY (Click Correctness)\n');
fprintf('==========================================================\n');
fprintf('  Avatar (correct):  n = %d,  median = %.2f s,  peak KDE = %.2f s\n', ...
    n_av, median(rt_avatar), peak_av_emp);
fprintf('  Shadow (wrong):    n = %d,  median = %.2f s,  peak KDE = %.2f s\n', ...
    n_sh, median(rt_shadow), peak_sh_emp);
fprintf('  Static (wrong):    n = %d,  median = %.2f s\n', n_st, median(rt_static));
fprintf('  None   (wrong):    n = %d,  median = %.2f s\n', n_no, median(rt_none));
fprintf('  ALL WRONG:         n = %d,  median = %.2f s,  peak KDE = %.2f s\n', ...
    n_wr, median(rt_wrong), peak_wr_emp);

fprintf('\n  CONTRAST 1 -- shadow vs avatar:\n');
fprintf('    peak diff = %+.2f s,  95%% CI [%+.2f, %+.2f],  bootstrap p = %.3f\n', ...
    peak_diff_sh_av, ci_sh_av(1), ci_sh_av(2), p_sh_av);
fprintf('    KS D = %.3f, p = %.3f;  median diff (perm) %+.2f s, p = %.3f\n', ...
    ks_sh, p_ks_sh, med_diff_sh, p_perm_sh);

fprintf('\n  CONTRAST 2 -- all-wrong vs avatar (primary):\n');
fprintf('    peak diff = %+.2f s,  95%% CI [%+.2f, %+.2f],  bootstrap p = %.3f\n', ...
    peak_diff_wr_av, ci_wr_av(1), ci_wr_av(2), p_wr_av);
fprintf('    KS D = %.3f, p = %.3f;  median diff (perm) %+.2f s, p = %.3f\n', ...
    ks_wr, p_ks_wr, med_diff_wr, p_perm_wr);

fprintf('\n  Interpretation:\n');
if ci_wr_av(1) < 0 && ci_wr_av(2) > 0
    fprintf('    All-wrong CI spans 0: NO evidence that wrong clicks are left-shifted.\n');
elseif ci_wr_av(2) < 0
    fprintf('    All-wrong CI strictly < 0: wrong clicks ARE left-shifted (classical SPRT).\n');
else
    fprintf('    All-wrong CI strictly > 0: wrong clicks RIGHT-shifted (unusual).\n');
end
fprintf('==========================================================\n');
fprintf('Done.\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function idx = argmax(v)
    [~, idx] = max(v);
end

function y = sensitivity_curve(t, A, offset, T_eff, lambda)
% P(t) = A * x * exp(-lambda * x), where x = (t - offset)/T_eff for t > offset,
% 0 otherwise.
    x = max((t - offset) / T_eff, 0);
    y = A .* x .* exp(-lambda .* x);
    y(t < offset) = 0;
end

function res = fit_sensitivity_grid(xi_kde, f_kde, T, lambda, tau_grid)
% Grid search over onset lag tau; at each tau, fit A by scalar projection.
    best.R2 = -Inf;
    for k = 1:length(tau_grid)
        tau = tau_grid(k);
        T_eff = T - tau;
        if T_eff < 10; continue; end
        mask = xi_kde > tau;
        if sum(mask) < 10; continue; end
        x = max((xi_kde(mask) - tau) / T_eff, eps);
        f = f_kde(mask);
        shape = x .* exp(-lambda .* x);
        A = dot(shape, f) / dot(shape, shape);
        yfit = A .* shape;
        SS_res = sum((f - yfit).^2);
        SS_tot = sum((f - mean(f)).^2);
        R2 = max(1 - SS_res / SS_tot, 0);
        if R2 > best.R2
            best.R2 = R2;
            best.offset = tau;
            best.A = A;
            best.T_eff = T_eff;
            best.peak_time = tau + T_eff / lambda;
        end
    end
    if ~isfinite(best.R2) || best.R2 < 0
        best.R2 = 0; best.offset = 0; best.A = 0;
        best.T_eff = T; best.peak_time = T / lambda;
    end
    res = best;
end
