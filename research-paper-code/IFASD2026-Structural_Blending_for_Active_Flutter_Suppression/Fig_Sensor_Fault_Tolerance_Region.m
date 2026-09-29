%% Plot: Sensor Fault Tolerance — Progressive Sensor Failures
%  Loads Data/FaultTolerance_SensorFailure.mat (from R2_Sensor_Fault_Tolerance.m)
%  Produces: Figures/Fig26_Sensor_Fault_Tolerance_Region.pdf
clearvars
textwidth = 15.98; % cm (IFASD 2026 template: 455.24pt)
figwidth  = 0.7*textwidth;
figheight = 0.5*figwidth;
docFontSize  = 9; % pt
stdLineWidth = 1.2;

script_path = mfilename('fullpath');
script_dir  = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end

% ---- Load results -------------------------------------------------------
data_dir = fullfile(script_dir, 'Data');
load(fullfile(data_dir, 'FaultTolerance_SensorFailure.mat'))

% ---- Figure: Flutter velocity under progressive sensor failures ---------
fig = figure('Name','Sensor Fault Tolerance');
set(fig, 'Units','centimeters', ...
    'Position', [5,5,figwidth,figheight], ...
    'defaultAxesFontSize', docFontSize, ...
    'defaultTextFontSize', docFontSize, ...
    'defaultLegendFontSize', docFontSize, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperSize', [figwidth figheight]);
set(fig, 'PaperPosition', [0 0 figwidth figheight]);

kv = (1:max_fail)';

% H2 envelope (min to max)
fill_x = [kv; flipud(kv)];
fill_y = [V_H2_max; flipud(V_H2_min)];
h_fill = fill(fill_x, fill_y, [0.85 0.85 0.85], 'EdgeColor', 'none');
hold on;

% H2 worst case (lower edge of envelope)
h_H2 = plot(kv, V_H2_min, 'k--s', 'LineWidth', stdLineWidth, ...
    'MarkerSize', 4, 'MarkerFaceColor', [0.6 0.6 0.6]);

% SB worst case (flat at nominal — identical for all combos)
h_SB = plot(kv, V_SB_min, 'k-o', 'LineWidth', stdLineWidth, ...
    'MarkerSize', 4, 'MarkerFaceColor', 'k');

% OL flutter reference
h_OL = yline(V_flutter_OL, 'k:', 'LineWidth', 1.0);

xlabel('Number of simultaneous sensor failures', 'Interpreter', 'latex');
ylabel('Flutter velocity [m/s]', 'Interpreter', 'latex');
legend([h_SB, h_H2, h_fill, h_OL], ...
    {'SB (reconfig.)', '$H_2$ worst case', '$H_2$ range', 'OL flutter'}, ...
    'Interpreter', 'latex', 'Location', 'best', 'FontSize', docFontSize-1);
xticks(kv);
xlim([0.5, max_fail+0.5]);
ylim([95 165]);
grid on;
set(fig, 'Color', 'w'); set(gca, 'Color', 'w');

FigureName = 'Fig26_Sensor_Fault_Tolerance_Region';
fullFigurePath = fullfile(figures_dir, FigureName);
print(fig, fullFigurePath, '-dpdf', '-vector');
% print(fig, fullFigurePath, '-dmeta', '-vector');
