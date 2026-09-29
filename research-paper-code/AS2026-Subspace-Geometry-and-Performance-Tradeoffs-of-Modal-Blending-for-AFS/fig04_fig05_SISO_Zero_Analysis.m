%% fig04_fig05_SISO_Zero_Analysis.m
% Paper Figs. 4 and 5. Generates: fig04_siso_numerator_geometry   (.pdf, .emf)
%                                 fig05_siso_performance_vs_delta (.pdf, .emf)
% EMF (-dmeta) is written on Windows only; other platforms get the PDF.
%
% Fig. 4 — Numerator geometry: three stacked panels showing the SISO
% blending parameters alpha (c'b), beta (c'Jb), and the induced
% transmission zero z(Delta) as functions of the quadrature angle Delta.
% Canonical angles and the H2 Blending direction are marked. Reveals how blending angle controls zero placement and when
% RHP zeros appear.
%
% Fig. 5 — Performance diagnostics: three stacked panels plotting
% closed-loop soft goal (log scale) against Delta, the induced zero z,
% and alpha. Points are coded as margin-met (filled circles) vs
% margin-failed (crosses), with best-performing Delta and H2 Blending baseline
% highlighted. Shows the link between zero location and achievable
% closed-loop performance.
%
% Requires: Data/RHP_zeros_analysis.mat (from R04_SISO_RHP_Zeros.m)

clearvars

% Figure settings
textwidth = 0.0351*372.0;    % cm/pt * latex template pt;  % cm
figwidth  = 0.9 * textwidth;
figheight = 0.7 * figwidth;
docFontSize  = 9;
stdLineWidth = 1.2;
mkSz = 8;

% Paths
script_dir  = fileparts(mfilename('fullpath'));
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end

% Load data
data_file = fullfile(script_dir, 'Data', 'RHP_zeros_analysis.mat');
if ~exist(data_file, 'file')
    error('Data file not found. Run R04_SISO_RHP_Zeros.m first.');
end
load(data_file, 'results');

% Unpack results
sigma           = results.sigma;
omega           = results.omega;
Delta_deg_sweep = results.Delta_deg_sweep;
sweep_alpha     = results.sweep_alpha;
sweep_beta      = results.sweep_beta;
sweep_z_formula = results.sweep_z_formula;
sweep_softGoal  = results.sweep_softGoal;
sweep_hardGoal  = results.sweep_hardGoal;
h2_sweep_idx    = results.h2_sweep_idx;

% Recompute derived quantities
Delta_zero_deg = atan2d(sigma, omega);
canonical_Delta_deg = [0, Delta_zero_deg, 90, 180];

[bestSoft, bestIdx] = min(sweep_softGoal);
bestDelta = Delta_deg_sweep(bestIdx);
ok = sweep_hardGoal <= 1;

% Greyscale palette
clrBlack     = [0 0 0];
clrDarkGray  = [0.35 0.35 0.35];
clrLightGray = [0.65 0.65 0.65];
clrCanonical = [0.5 0.5 0.5];

%% Fig. 4 — Numerator Geometry and Induced Zero

fig1 = figure('Name', 'Numerator geometry and induced zero', 'Color', 'w');
set(fig1, 'Units', 'centimeters', ...
    'Position', [2, 2, figwidth, figheight], ...
    'defaultAxesFontSize', docFontSize, ...
    'defaultTextFontSize', docFontSize, ...
    'defaultLegendFontSize', docFontSize, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');
set(fig1, 'PaperUnits', 'centimeters');
set(fig1, 'PaperSize', [figwidth figheight]);
set(fig1, 'PaperPosition', [0 0 figwidth figheight]);

tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- alpha ---
nexttile;
h_alpha = plot(Delta_deg_sweep, sweep_alpha, 'k-', 'LineWidth', stdLineWidth);
hold on;
for ic = 1:4
    xline(canonical_Delta_deg(ic), '--', 'Color', clrCanonical, 'HandleVisibility', 'off');
end
h_p1 = plot(Delta_deg_sweep(h2_sweep_idx(1)), sweep_alpha(h2_sweep_idx(1)), ...
    'v', 'MarkerSize', mkSz, 'MarkerFaceColor', clrDarkGray, 'MarkerEdgeColor', 'k');
h_best = plot(bestDelta, sweep_alpha(bestIdx), 'kp', 'MarkerSize', 10, ...
    'MarkerFaceColor', clrBlack, 'MarkerEdgeColor', 'k');
grid on; grid minor;
ax = gca; ax.Box = 'on';
ylabel('$\alpha = c^\top b$');
lgd = legend([h_p1, h_best], {'$H_2$ Blending', sprintf('Best ($\\Delta = %g^\\circ$)', bestDelta)}, 'Location', 'southwest');
lgd.AutoUpdate = 'off';

% --- beta ---
nexttile;
plot(Delta_deg_sweep, sweep_beta, 'k-', 'LineWidth', stdLineWidth);
hold on;
for ic = 1:4
    xline(canonical_Delta_deg(ic), '--', 'Color', clrCanonical, 'HandleVisibility', 'off');
end
plot(Delta_deg_sweep(h2_sweep_idx(1)), sweep_beta(h2_sweep_idx(1)), ...
    'v', 'MarkerSize', mkSz, 'MarkerFaceColor', clrDarkGray, 'MarkerEdgeColor', 'k');
plot(bestDelta, sweep_beta(bestIdx), 'kp', 'MarkerSize', 10, ...
    'MarkerFaceColor', clrBlack, 'MarkerEdgeColor', 'k');
grid on; grid minor;
ax = gca; ax.Box = 'on';
ylabel('$\beta = c^\top J b$');

% --- zero ---
nexttile;
plot(Delta_deg_sweep, sweep_z_formula, 'k-', 'LineWidth', stdLineWidth);
hold on;
for ic = 1:4
    xline(canonical_Delta_deg(ic), '--', 'Color', clrCanonical, 'HandleVisibility', 'off');
end
plot(Delta_deg_sweep(h2_sweep_idx(1)), sweep_z_formula(h2_sweep_idx(1)), ...
    'v', 'MarkerSize', mkSz, 'MarkerFaceColor', clrDarkGray, 'MarkerEdgeColor', 'k');
plot(bestDelta, sweep_z_formula(bestIdx), 'kp', 'MarkerSize', 10, ...
    'MarkerFaceColor', clrBlack, 'MarkerEdgeColor', 'k');
yline(0, ':', 'Color', clrCanonical, 'HandleVisibility', 'off');
grid on; grid minor;
ax = gca; ax.Box = 'on';
xlabel('$\Delta$ [deg]');
ylabel('$z(\Delta)$');
ylim([-80, 80]);

figureName = 'fig04_siso_numerator_geometry';
print(fig1, fullfile(figures_dir, figureName), '-dpdf', '-vector');
if ispc, print(fig1, fullfile(figures_dir, figureName), '-dmeta', '-vector'); end  % EMF export is Windows-only

%% Fig. 5 — Controller Performance Diagnostics

fig2 = figure('Name', 'Closed-loop performance diagnostics', 'Color', 'w');
set(fig2, 'Units', 'centimeters', ...
    'Position', [2 + figwidth + 1, 2, figwidth, figheight], ...
    'defaultAxesFontSize', docFontSize, ...
    'defaultTextFontSize', docFontSize, ...
    'defaultLegendFontSize', docFontSize, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');
set(fig2, 'PaperUnits', 'centimeters');
set(fig2, 'PaperSize', [figwidth figheight]);
set(fig2, 'PaperPosition', [0 0 figwidth figheight]);

tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- performance vs Delta ---
nexttile;
h_ok = semilogy(Delta_deg_sweep(ok), sweep_softGoal(ok), 'ko', ...
    'MarkerSize', 4, 'MarkerFaceColor', clrBlack, 'LineWidth', 0.5);
hold on;
h_fail = semilogy(Delta_deg_sweep(~ok), sweep_softGoal(~ok), 'kx', ...
    'MarkerSize', 5, 'LineWidth', 0.8);
h_best = semilogy(bestDelta, bestSoft, 'kp', 'MarkerSize', 10, ...
    'MarkerFaceColor', clrBlack, 'MarkerEdgeColor', 'k');
xline(Delta_deg_sweep(h2_sweep_idx(1)), '--', 'Color', clrDarkGray, 'HandleVisibility', 'off');
h_p2 = semilogy(Delta_deg_sweep(h2_sweep_idx(1)), sweep_softGoal(h2_sweep_idx(1)), ...
    'v', 'MarkerSize', mkSz, 'MarkerFaceColor', clrDarkGray, 'MarkerEdgeColor', 'k');
for ic = 1:4
    xline(canonical_Delta_deg(ic), ':', 'Color', clrCanonical, 'HandleVisibility', 'off');
end
grid on; grid minor;
ax = gca; ax.Box = 'on';
xlabel('$\Delta$ [deg]');
ylabel('Soft goal');

% --- performance vs zero ---
nexttile;
h2_ok = semilogy(sweep_z_formula(ok), sweep_softGoal(ok), 'ko', ...
    'MarkerSize', 4, 'MarkerFaceColor', clrBlack, 'LineWidth', 0.5);
hold on;
h2_fail = semilogy(sweep_z_formula(~ok), sweep_softGoal(~ok), 'kx', ...
    'MarkerSize', 5, 'LineWidth', 0.8);
h2_best = semilogy(sweep_z_formula(bestIdx), bestSoft, 'kp', 'MarkerSize', 10, ...
    'MarkerFaceColor', clrBlack, 'MarkerEdgeColor', 'k');
h2_p2 = semilogy(sweep_z_formula(h2_sweep_idx(1)), sweep_softGoal(h2_sweep_idx(1)), ...
    'v', 'MarkerSize', mkSz, 'MarkerFaceColor', clrDarkGray, 'MarkerEdgeColor', 'k');
xline(0, '--k', 'HandleVisibility', 'off');
grid on; grid minor;
ax = gca; ax.Box = 'on';
xlabel('Induced zero $z$');
ylabel('Soft goal');
lgd = legend([h2_ok, h2_fail, h2_best, h2_p2], ...
    {'Margins met', 'Margins failed', ...
     sprintf('Best ($\\Delta = %g^\\circ$)', bestDelta), ...
     '$H_2$ Blending'}, ...
    'Location', 'southeast');
lgd.AutoUpdate = 'off';

% --- performance vs alpha ---
nexttile;
semilogy(sweep_alpha(ok), sweep_softGoal(ok), 'ko', ...
    'MarkerSize', 4, 'MarkerFaceColor', clrBlack, 'LineWidth', 0.5);
hold on;
semilogy(sweep_alpha(~ok), sweep_softGoal(~ok), 'kx', ...
    'MarkerSize', 5, 'LineWidth', 0.8);
semilogy(sweep_alpha(bestIdx), bestSoft, 'kp', 'MarkerSize', 10, ...
    'MarkerFaceColor', clrBlack, 'MarkerEdgeColor', 'k');
semilogy(sweep_alpha(h2_sweep_idx(1)), sweep_softGoal(h2_sweep_idx(1)), ...
    'v', 'MarkerSize', mkSz, 'MarkerFaceColor', clrDarkGray, 'MarkerEdgeColor', 'k');
grid on; grid minor;
ax = gca; ax.Box = 'on';
xlabel('$\alpha = c^\top b$');
ylabel('Soft goal');

figureName = 'fig05_siso_performance_vs_delta';
print(fig2, fullfile(figures_dir, figureName), '-dpdf', '-vector');
if ispc, print(fig2, fullfile(figures_dir, figureName), '-dmeta', '-vector'); end  % EMF export is Windows-only

fprintf('2 figures saved to %s\n', figures_dir);
