%% Figures: Geometric Mode Isolation — Coupling Matrix Comparison
%  Loads results from R1_geometric_mode_isolation_structural.m
%  and plots output/input coupling bar charts.
clearvars
textwidth = 15.98; % cm (IFASD 2026 template)
figwidth  = 0.47*textwidth;
figheight = 0.5*figwidth;
docFontSize = 9; % pt

script_path = mfilename('fullpath');
script_dir  = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end

% ---- Load analysis results -----------------------------------------------
data_dir = fullfile(script_dir, 'Data');
load(fullfile(data_dir, 'R1_geometric_mode_isolation_results.mat'));

%% ======================================================================
%  FIGURE 1: Output coupling matrix comparison
%  ======================================================================
fig1 = figure('Name','Output Coupling Comparison');
set(fig1, 'Units','centimeters', ...
         'Position', [7,7,figwidth,figheight], ...
         'defaultAxesFontSize', docFontSize, ...
         'defaultTextFontSize', docFontSize, ...
         'defaultLegendFontSize', docFontSize, ...
         'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex');
set(fig1, 'PaperUnits', 'centimeters');
set(fig1, 'PaperSize', [figwidth figheight]);
set(fig1, 'PaperPosition', [0 0 figwidth figheight]);

t1 = tiledlayout(fig1, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile(t1);
b1 = bar(ax1, abs(coupling_pinv'));
b1(1).FaceColor = [0.3 0.3 0.3];
b1(2).FaceColor = [0.7 0.7 0.7];
xticklabels(ax1, {'$f_1$','$f_2$','$f_3$','$f_4$','$f_5$'});
ylabel(ax1, 'Coupling magnitude', 'Interpreter', 'latex');
title(ax1, 'Pseudoinverse', 'Interpreter', 'latex', 'FontSize', docFontSize);
legend(ax1, {'$\hat{\eta}_1$','$\hat{\eta}_2$'}, 'Interpreter', 'latex', 'Location', 'northeast');
grid(ax1, 'on');
set(ax1, 'Color', 'w');

ax2 = nexttile(t1);
b2 = bar(ax2, abs(coupling_iso'));
b2(1).FaceColor = [0.3 0.3 0.3];
b2(2).FaceColor = [0.7 0.7 0.7];
xticklabels(ax2, {'$f_1$','$f_2$','$f_3$','$f_4$','$f_5$'});
ylabel(ax2, 'Coupling magnitude', 'Interpreter', 'latex');
title(ax2, 'Geometric isolation', 'Interpreter', 'latex', 'FontSize', docFontSize);
legend(ax2, {'$\hat{\eta}_1$','$\hat{\eta}_2$'}, 'Interpreter', 'latex', 'Location', 'northeast');
grid(ax2, 'on');
set(ax2, 'Color', 'w');

set(fig1, 'Color', 'w');
print(fig1, fullfile(figures_dir, 'Fig16_Output_Coupling_Pinv_vs_GeomIso'), '-dpdf', '-vector');

%% ======================================================================
%  FIGURE 2: Input coupling matrix comparison
%  ======================================================================
fig2 = figure('Name','Input Coupling Comparison');
set(fig2, 'Units','centimeters', ...
         'Position', [7,7,figwidth,figheight], ...
         'defaultAxesFontSize', docFontSize, ...
         'defaultTextFontSize', docFontSize, ...
         'defaultLegendFontSize', docFontSize, ...
         'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex');
set(fig2, 'PaperUnits', 'centimeters');
set(fig2, 'PaperSize', [figwidth figheight]);
set(fig2, 'PaperPosition', [0 0 figwidth figheight]);

t2 = tiledlayout(fig2, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile(t2);
b1 = bar(ax1, abs(coupling_ku_pinv));
b1(1).FaceColor = [0.3 0.3 0.3];
b1(2).FaceColor = [0.7 0.7 0.7];
xticklabels(ax1, {'$f_1$','$f_2$','$f_3$','$f_4$','$f_5$'});
ylabel(ax1, 'Coupling magnitude', 'Interpreter', 'latex');
title(ax1, 'Pseudoinverse', 'Interpreter', 'latex', 'FontSize', docFontSize);
legend(ax1, {'$u_1$','$u_2$'}, 'Interpreter', 'latex', 'Location', 'northeast');
grid(ax1, 'on');
set(ax1, 'Color', 'w');

ax2 = nexttile(t2);
b2 = bar(ax2, abs(coupling_ku_iso));
b2(1).FaceColor = [0.3 0.3 0.3];
b2(2).FaceColor = [0.7 0.7 0.7];
xticklabels(ax2, {'$f_1$','$f_2$','$f_3$','$f_4$','$f_5$'});
ylabel(ax2, 'Coupling magnitude', 'Interpreter', 'latex');
title(ax2, 'Geometric isolation', 'Interpreter', 'latex', 'FontSize', docFontSize);
legend(ax2, {'$u_1$','$u_2$'}, 'Interpreter', 'latex', 'Location', 'northeast');
grid(ax2, 'on');
set(ax2, 'Color', 'w');

set(fig2, 'Color', 'w');
print(fig2, fullfile(figures_dir, 'Fig17_Input_Coupling_Pinv_vs_GeomIso'), '-dpdf', '-vector');
