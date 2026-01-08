%% V=60 m/s

clearvars

textwidth = 16.5730; % cm
figwidth = 0.47*textwidth; % cm
figheight = 0.7*figwidth;
docFontSize = 9; % pt
stdMarkerSize = 12;

fig = figure('Name','Pole Map 60 m/s');
t = tiledlayout(fig, 1, 1, 'Padding', 'none', 'TileSpacing', 'none');
set(fig, 'Units','centimeters', ...
         'Position', [7,7,figwidth,figheight], ...
         'defaultAxesFontSize', docFontSize, ...
         'defaultTextFontSize', docFontSize, ...
         'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex');
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperSize', [figwidth figheight]);        % Set PDF page size to match figure
set(fig, 'PaperPosition', [0 0 figwidth figheight]); % Position plot to fill the PDF page exactly
ax = nexttile;
figureName = 'PoleMap_V60';
V_inf=60;
imuIDX = [4,8];
ailIDX = 4;
G = build_G_RectWing(V_inf,imuIDX,ailIDX);
EV=pole(G);

for i = 1:length(EV)
plot(real(EV(i)), imag(EV(i)),"k.",'MarkerSize',stdMarkerSize);
hold on
end
xline(0, 'k--', 'stability limit', 'Interpreter', 'latex');
% title(['Eigenvalue Loci at $V_\infty = ',num2str(V_inf),' m/s$'], 'Interpreter', 'latex'); 
% title(' Eigenvalue Loci at $V_\infty = 60 m/s$', 'Interpreter', 'latex'); 
xlabel('Real Part $\Re (\lambda)$', 'Interpreter', 'latex'); 
ylabel('Imaginary Part $\Im (\lambda)$', 'Interpreter', 'latex'); 
grid on;
grid minor;
% axis([-25,5,-200,200]);
ax.XLim = [-25, 5];
ax.YLim = [-200, 200];
% ax = fig.CurrentAxes; % Get current axes
ax.TickLabelInterpreter = 'latex'; % Set tick labels to LaTeX
set(ax, 'Color', 'w'); % Set axes background to white
set(fig, 'Color', 'w'); % Set figure background to white

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
fullFigurePath = fullfile(figures_dir, figureName);
print(fig, fullFigurePath, '-dpdf', '-vector');
print(fig, fullFigurePath, '-dmeta', '-vector');

%% V=120 m/s

clearvars

textwidth = 16.5730; % cm
figwidth = 0.47*textwidth; % cm
figheight = 0.7*figwidth;
docFontSize = 9; % pt
stdMarkerSize = 12;

fig = figure('Name','Pole Map 120 m/s');
t = tiledlayout(fig, 1, 1, 'Padding', 'none', 'TileSpacing', 'none');
set(fig, 'Units','centimeters', ...
         'Position', [7,7,figwidth,figheight], ...
         'defaultAxesFontSize', docFontSize, ...
         'defaultTextFontSize', docFontSize, ...
         'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex');
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperSize', [figwidth figheight]);        % Set PDF page size to match figure
set(fig, 'PaperPosition', [0 0 figwidth figheight]); % Position plot to fill the PDF page exactly
ax = nexttile;
figureName = 'PoleMap_V120';
V_inf=120;
imuIDX = [4,8];
ailIDX = 4;
G = build_G_RectWing(V_inf,imuIDX,ailIDX);
EV=pole(G);

for i = 1:length(EV)
plot(real(EV(i)), imag(EV(i)),"k.",'MarkerSize',stdMarkerSize);
hold on
end
xline(0, 'k--', 'stability limit', 'Interpreter', 'latex');
% title(' Eigenvalue Loci at $V_\infty = 120 m/s$', 'Interpreter', 'latex'); 
% title(' Eigenvalue Loci at $V_\infty = 60 m/s$', 'Interpreter', 'latex'); 
xlabel('Real Part $\Re (\lambda)$', 'Interpreter', 'latex'); 
ylabel('Imaginary Part $\Im (\lambda)$', 'Interpreter', 'latex'); 
grid on;
grid minor;
% axis([-25,5,-200,200])
ax.XLim = [-25, 5];
ax.YLim = [-200, 200];
% ax = fig.CurrentAxes; % Get current axes
ax.TickLabelInterpreter = 'latex'; % Set tick labels to LaTeX
set(ax, 'Color', 'w'); % Set axes background to white
set(fig, 'Color', 'w'); % Set figure background to white

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
fullFigurePath = fullfile(figures_dir, figureName);
print(fig, fullFigurePath, '-dpdf', '-vector');
print(fig, fullFigurePath, '-dmeta', '-vector');