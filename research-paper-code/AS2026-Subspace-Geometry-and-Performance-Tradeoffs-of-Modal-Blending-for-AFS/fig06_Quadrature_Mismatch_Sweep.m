%% fig06_Quadrature_Mismatch_Sweep.m
% Paper Fig. 6. Generates: fig06_flutterwing_delta_sweep (.pdf, .emf)
% EMF (-dmeta) is written on Windows only; other platforms get the PDF.
%
% Linear-scale plot of max soft goal versus quadrature angle Delta on the
% full FlutterWing model (8 actuators, 8 sensors). The one-parameter family
% curve (output direction fixed, input direction rotated) is the R05 sweep.
% Points where the hard stability margins are violated appear in a grey
% shaded band with cross markers; feasible points are joined by a solid line; off-chart points are clipped and marked with
% upward arrows. The named direction pairs are drawn at their family
% coordinates with their locked R06 synthesis values, off the family curve:
% the best restricted member (star), the H2-optimal pair (triangle and
% dashed line) and the margin-infeasible MVF z = 0 pair (diamond, clipped).
%
% Requires: Data/Quadrature_Mismatch_DeltaSweep.mat (from R05_Quadrature_Mismatch.m)
%           Data/Locked_Comparison_Results.mat

clearvars

% Figure settings
textwidth = 0.0351*372.0;    % cm/pt * latex template pt;  % cm
figwidth  = 0.9 * textwidth;
figheight = 0.5 * figwidth;
docFontSize = 9;
mkSz = 8;
archive = 'Locked_Comparison_Results.mat';   % or Locked_Comparison_Results_test.mat

% Paths
script_dir  = fileparts(mfilename('fullpath'));
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end

% Load data
data_file = fullfile(script_dir, 'Data', archive);
if ~exist(data_file, 'file')
    error('Data file not found. Run R06_Locked_Comparison.m first.');
end
A = load(data_file, 'synth', 'model');

% Family curve: the R05 sweep
L = load(fullfile(script_dir, 'Data', 'Quadrature_Mismatch_DeltaSweep.mat'), 'results');
Delta_synth_deg = L.results.Delta_synth_deg;
synth_maxSoft   = L.results.synth_maxSoft;
synth_hardGoal  = L.results.synth_hardGoal;

% Named pairs: locked values at their family coordinates
caseIds = {A.model.cases.id};
famDelta = @(id) A.model.cases(strcmp(caseIds, id)).channel.familyDisplayDeltaDeg;
bestLocked = A.synth.R1bestRestricted.selected.MaxSoft;
bestDelta  = A.model.family.bestFamilyDeltaDeg;
h2Locked   = A.synth.H2Pusch.selected.MaxSoft;
h2Delta    = famDelta('H2Pusch');
mvfDelta   = famDelta('MVFz0');       % margin-infeasible: clipped off-chart

ok = synth_hardGoal <= 1;

% Greyscale palette
clrBlack     = [0 0 0];
clrDarkGray  = [0.35 0.35 0.35];
clrLightGray = [0.65 0.65 0.65];

%% Create Figure
fig = figure('Name', 'Delta sweep performance', 'Color', 'w');
set(fig, 'Units', 'centimeters', ...
    'Position', [2, 2, figwidth, figheight], ...
    'defaultAxesFontSize', docFontSize, ...
    'defaultTextFontSize', docFontSize, ...
    'defaultLegendFontSize', docFontSize, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperSize', [figwidth figheight]);
set(fig, 'PaperPosition', [0 0 figwidth figheight]);

% Clip y-axis for linear scale — hide extreme outliers
yClip = 5;
synth_clipped = min(synth_maxSoft, yClip);
clipped = synth_maxSoft > yClip;

% Failure band shading (hard margin > 1)
fail_idx = find(~ok);
if ~isempty(fail_idx)
    % Extend band half a step beyond outermost failing points
    dStep = median(diff(Delta_synth_deg));
    bandL = Delta_synth_deg(fail_idx(1)) - dStep/2;
    bandR = Delta_synth_deg(fail_idx(end)) + dStep/2;
    h_band = patch([bandL bandR bandR bandL], [0 0 yClip yClip], ...
        [0.9 0.9 0.9], 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
hold on;

% Connected line plot for all feasible family points (clipped)
h_line = plot(Delta_synth_deg(ok), synth_clipped(ok), '-o', ...
    'Color', clrBlack, 'MarkerSize', 4, 'MarkerFaceColor', clrBlack, 'LineWidth', 0.8);

% Failed family points (different marker, no connecting line to OK points)
if any(~ok)
    h_fail = plot(Delta_synth_deg(~ok), synth_clipped(~ok), 'kx', ...
        'MarkerSize', 6, 'LineWidth', 1.2);
end

% Upward arrows for clipped (off-chart) points
clipped_idx = find(clipped(:))';
for ii = clipped_idx
    plot(Delta_synth_deg(ii), yClip, 'k^', 'MarkerSize', 5, ...
        'MarkerFaceColor', clrBlack, 'HandleVisibility', 'off');
end

% Best restricted family member: locked re-synthesis value
h_best = plot(bestDelta, min(bestLocked, yClip), 'kp', 'MarkerSize', 10, ...
    'MarkerFaceColor', clrBlack, 'MarkerEdgeColor', 'k');

% H2-optimal pair: distinct marker off the family curve (locked value)
xline(h2Delta, '--', 'Color', clrDarkGray, 'HandleVisibility', 'off');
h_p = plot(h2Delta, min(h2Locked, yClip), ...
    'v', 'MarkerSize', mkSz, 'MarkerFaceColor', clrDarkGray, 'MarkerEdgeColor', 'k');

% MVF z = 0 pair: margin-infeasible, clipped off-chart
h_mvf = plot(mvfDelta, yClip, 'd', 'MarkerSize', mkSz, ...
    'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
plot(mvfDelta, yClip, 'k^', 'MarkerSize', 4, 'MarkerFaceColor', clrBlack, ...
    'HandleVisibility', 'off');

grid on; grid minor;
ax = gca; ax.Box = 'on';
xlim([0 180]);
ylim([0.9 yClip]);
xlabel('$\Delta$ [deg]');
ylabel('Max soft goal');

% Build legend dynamically to avoid empty-handle mismatch
lgd_h = h_line;
lgd_s = {'Family (margins met)'};
if any(~ok)
    lgd_h(end+1) = h_fail;
    lgd_s{end+1} = 'Family (margins violated)';
end
lgd_h(end+1) = h_best;
lgd_s{end+1} = sprintf('Best restricted ($\\Delta = %.0f^\\circ$)', bestDelta);
lgd_h(end+1) = h_p;
lgd_s{end+1} = '$H_2$-optimal pair';
lgd_h(end+1) = h_mvf;
lgd_s{end+1} = 'MVF $z=0$ pair (infeasible)';
lgd = legend(lgd_h, lgd_s, 'Location', 'northeast');
lgd.AutoUpdate = 'off';

%% Export
figureName = 'fig06_flutterwing_delta_sweep';
print(fig, fullfile(figures_dir, figureName), '-dpdf', '-vector');
if ispc, print(fig, fullfile(figures_dir, figureName), '-dmeta', '-vector'); end  % EMF export is Windows-only
fprintf('Figure saved to Figures/%s.pdf\n', figureName);
