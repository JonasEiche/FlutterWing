clearvars
textwidth = 16.5730; % cm
figwidth = 0.7*textwidth; % cm
figheight = 0.5*figwidth;
docFontSize = 9; % pt
stdLineWidth = 1.2; 

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end

E1_y = -1;
E2_y = 1;

S_1 = @(y, E1_y, E2_y) -((y - E2_y).^2.*(2*y - 3*E1_y + E2_y))./(E1_y - E2_y).^3;
S_2 = @(y, E1_y, E2_y) ((y - E1_y).*(y - E2_y).^2)./(E1_y - E2_y).^2;
S_3 = @(y, E1_y, E2_y) (y - E2_y)./(E1_y - E2_y);
S_4 = @(y, E1_y, E2_y) ((y - E1_y).^2.*(2*y + E1_y - 3*E2_y))./(E1_y - E2_y).^3;
S_5 = @(y, E1_y, E2_y) ((y - E1_y).^2.*(y - E2_y))./(E1_y - E2_y).^2;
S_6 = @(y, E1_y, E2_y) -(y - E1_y)./(E1_y - E2_y);


x_=linspace(-1,1,32);

S_1_x_ = S_1(x_, E1_y, E2_y);
S_2_x_ = S_2(x_, E1_y, E2_y);
S_4_x_ = S_4(x_, E1_y, E2_y);
S_5_x_ = S_5(x_, E1_y, E2_y);


fig = figure('Name','FEM 4DoF Weight Functions');
set(fig, 'Units','centimeters', ...
         'Position', [7,7,figwidth,figheight], ...
         'defaultAxesFontSize', docFontSize, ...
         'defaultTextFontSize', docFontSize, ...
         'defaultLegendFontSize', docFontSize, ...
         'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex');
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperSize', [figwidth figheight]);        % Set PDF page size to match figure
set(fig, 'PaperPosition', [0 0 figwidth figheight]); % Position plot to fill the PDF page exactly
plot(x_,S_1_x_, 'k-', ...
    x_,S_2_x_, 'k:', ...
    x_,S_4_x_, 'k-.', ...
    x_,S_5_x_, 'k--','LineWidth',stdLineWidth)
% title('Weight Functions of 4DoF Beam Element', 'Interpreter', 'latex', 'FontSize', docFontSize)
xlabel('Normalized Length', 'Interpreter', 'latex', 'FontSize', docFontSize) 
ylabel('Normalized Amplitude', 'Interpreter', 'latex', 'FontSize', docFontSize) 
legend({'H1','H2','H3','H4'},'Location','northeast', 'Interpreter', 'latex', 'FontSize', docFontSize)
grid on
ax = fig.CurrentAxes; % Get current axes
ax.TickLabelInterpreter = 'latex'; % Set tick labels to LaTeX
ax.FontSize = docFontSize;
set(ax, 'Color', 'w'); % Set axes background to white
set(fig, 'Color', 'w'); % Set figure background to white

NyquistPlotName = 'FEM_WeightFunctions';
fullVpzmap_PlotPath = fullfile(figures_dir, NyquistPlotName);
print(fig, fullVpzmap_PlotPath, '-dpdf', '-vector');
print(fig, fullVpzmap_PlotPath, '-dmeta', '-vector');