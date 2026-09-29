%% fig07_Condition_Numbers.m
% Paper Fig. 7. Generates: fig07a_observability_condition   (.pdf, .emf)
%                          fig07b_controllability_condition (.pdf, .emf)
% EMF (-dmeta) is written on Windows only; other platforms get the PDF.
%
% Two greyscale line plots of the directional conditioning of the flutter
% and residual modes versus freestream velocity:
%   fig07a: kappa(C_S) = sigma_max/sigma_min of the modal output map
%   fig07b: kappa(U_S) = sigma_max/sigma_min of the modal input map
%
% Requires: Data/Modal_Obsrv_Contr_data.mat  (from R01_Modal_Observability.m)

clearvars

%% 1. Figure Settings
textwidth    = 0.0351*372.0;    % cm/pt * latex template pt
figwidth     = 0.47 * textwidth;   % paired side-by-side
figheight    = 0.7 * figwidth;
docFontSize  = 9;
stdLineWidth = 1.2;

%% 2. Paths
script_dir  = fileparts(mfilename('fullpath'));
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end

%% 3. Load R01 Data
data_file = fullfile(script_dir, 'Data', 'Modal_Obsrv_Contr_data.mat');
if ~exist(data_file, 'file')
    error('Data file not found. Run R01_Modal_Observability.m first.');
end
load(data_file, 'V_inf_arr', 'Cond_B', 'Cond_C');

iFlut2D = 7;   % index in R01 set_names: 'flutter 2D'
iRes2D  = 8;   % index in R01 set_names: 'residual 2D'

%% 4. Plot Both Figures
% Greyscale: flutter = black, residual = dark grey (matching Fig. 3)
colFlut = [0.0 0.0 0.0];
colRes  = [0.35 0.35 0.35];

figureNames = {'fig07a_observability_condition', ...
               'fig07b_controllability_condition'};
condData    = {Cond_C, Cond_B};
yLabels     = {'$\kappa(C_S)$', '$\kappa(U_S)$'};

for iFig = 1:2
    fig = figure('Name', figureNames{iFig});
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

    plot(ax, V_inf_arr, condData{iFig}(iFlut2D,:), '-',  'Color', colFlut, 'LineWidth', stdLineWidth);
    plot(ax, V_inf_arr, condData{iFig}(iRes2D,:),  '--', 'Color', colRes,  'LineWidth', stdLineWidth);
    yline(1, 'k:', 'LineWidth', 0.8);

    xlabel('Velocity $V_\infty$ [m/s]', 'Interpreter', 'latex', 'FontSize', docFontSize);
    ylabel(yLabels{iFig}, 'Interpreter', 'latex', 'FontSize', docFontSize);
    if iFig == 2
        legend({'flutter 2D', 'residual 2D'}, ...
               'Location', 'best', 'Interpreter', 'latex', 'FontSize', docFontSize);
    end

    grid on
    grid minor
    set(fig, 'Color', 'w');
    ax.TickLabelInterpreter = 'latex';
    set(ax, 'Color', 'w');
    ax.YMinorTick = 'on';

    fullFigurePath = fullfile(figures_dir, figureNames{iFig});
    print(fig, fullFigurePath, '-dpdf', '-vector');
    if ispc, print(fig, fullFigurePath, '-dmeta', '-vector'); end  % EMF export is Windows-only
end

fprintf('Done. 2 figures saved to %s\n', figures_dir);
