%% fig03_Sigmin_Velocity_Sweep.m
% Paper Fig. 3. Generates: fig03a_sigmin_observability   (.pdf, .emf)
%                          fig03b_sigmin_controllability (.pdf, .emf)
% EMF (-dmeta) is written on Windows only; other platforms get the PDF.
%
% Two greyscale line plots of minimum singular values (sigma_min) versus
% freestream velocity V_inf, with actuator dynamics included. Each figure
% shows 3 MIMO blending sets (flutter 2D, residual 2D, flutter+residual 4D)
% using distinct line styles. Illustrates how modal observability and
% controllability vary across the velocity envelope for different subspace
% selections.
%
% Requires: Data/Modal_Obsrv_Contr_data.mat (from R01_Modal_Observability.m)

clearvars

% Figure settings
textwidth = 0.0351*372.0;    % cm/pt * latex template pt;  % cm
figwidth  = 0.47 * textwidth;
figheight = 0.7 * figwidth;
docFontSize  = 9;
stdLineWidth = 1.2;

% Paths
script_dir  = fileparts(mfilename('fullpath'));
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end

% Load data
data_file = fullfile(script_dir, 'Data', 'Modal_Obsrv_Contr_data.mat');
if ~exist(data_file, 'file')
    error('Data file not found. Run R01_Modal_Observability.m first.');
end
load(data_file, 'V_inf_arr', 'set_names', 'SigMin_C', 'SigMin_B');

%% Plot 2 sigma_min Figures
ioLabels     = {'Observability', 'Controllability'};
figureNames  = {'fig03a_sigmin_observability', 'fig03b_sigmin_controllability'};
sigmin_data  = {SigMin_C, SigMin_B};
setIdx       = 7:9;

% Greyscale line styles for MIMO (3 lines)
mimoColors = [0.0 0.0 0.0;    % flutter 2D - black
              0.35 0.35 0.35;  % residual 2D - dark grey
              0.65 0.65 0.65]; % flutter+res 4D - light grey
mimoStyles = {'-', '--', ':'};

for iIO = 1:2
    data = sigmin_data{iIO};
    figureName = figureNames{iIO};

    fig = figure('Name', figureName);
    t = tiledlayout(fig, 1, 1, 'Padding', 'none', 'TileSpacing', 'none');
    set(fig, 'Units', 'centimeters', ...
             'Position', [7, 7, figwidth, figheight], ...
             'defaultAxesFontSize', docFontSize, ...
             'defaultTextFontSize', docFontSize, ...
             'defaultLegendFontSize', docFontSize, ...
             'defaultTextInterpreter', 'latex', ...
             'defaultAxesTickLabelInterpreter', 'latex', ...
             'defaultLegendInterpreter', 'latex');
    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', [figwidth figheight]);
    set(fig, 'PaperPosition', [0 0 figwidth figheight]);
    ax = nexttile;
    hold(ax, 'on');

    legendEntries = {};
    for k = 1:length(setIdx)
        si = setIdx(k);
        plot(ax, V_inf_arr, data(si,:), mimoStyles{k}, ...
             'Color', mimoColors(k,:), 'LineWidth', stdLineWidth);
        legendEntries{end+1} = set_names{si}; %#ok<SAGROW>
    end

    xlabel('Velocity $V_\infty$ [m/s]', 'Interpreter', 'latex', 'FontSize', docFontSize);
    ylabel(['$\sigma_{\min}$ -- ', ioLabels{iIO}], ...
           'Interpreter', 'latex', 'FontSize', docFontSize);
    if iIO == 2
        legend(legendEntries, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', docFontSize);
    end

    grid on
    grid minor
    set(fig, 'Color', 'w');
    ax.TickLabelInterpreter = 'latex';
    set(ax, 'Color', 'w');
    ax.YMinorTick = 'on';

    fullFigurePath = fullfile(figures_dir, figureName);
    print(fig, fullFigurePath, '-dpdf', '-vector');
    if ispc, print(fig, fullFigurePath, '-dmeta', '-vector'); end  % EMF export is Windows-only
end

fprintf('Done. 2 figures saved to %s\n', figures_dir);
