%% Pole maps at V=60 m/s and V=120 m/s
clearvars

textwidth = 15.98; % cm (IFASD 2026 template: 455.24pt)
figwidth = 0.47*textwidth; % cm
figheight = 0.7*figwidth;
docFontSize = 9; % pt
stdMarkerSize = 12;

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end

imuIDX = 1:8;
ailIDX = 1:8;

%% V=60 m/s
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
set(fig, 'PaperSize', [figwidth figheight]);
set(fig, 'PaperPosition', [0 0 figwidth figheight]);
ax = nexttile;
figureName = 'Fig09a_PoleMap_V60';
V_inf=60;
G = build_G_RectWing(V_inf,imuIDX,ailIDX);
EV=pole(G);

for i = 1:length(EV)
    plot(real(EV(i)), imag(EV(i)),"k.",'MarkerSize',stdMarkerSize);
    hold on
end
xline(0, 'k--', 'stability limit', 'Interpreter', 'latex');
xlabel('Real Part $\Re (\lambda)$', 'Interpreter', 'latex');
ylabel('Imaginary Part $\Im (\lambda)$', 'Interpreter', 'latex');
grid on;
grid minor;
ax.XLim = [-25, 5];
ax.YLim = [-200, 200];
ax.TickLabelInterpreter = 'latex';
set(ax, 'Color', 'w');
set(fig, 'Color', 'w');

fullFigurePath = fullfile(figures_dir, figureName);
print(fig, fullFigurePath, '-dpdf', '-vector');
% print(fig, fullFigurePath, '-dmeta', '-vector');

%% V=120 m/s
fig2 = figure('Name','Pole Map 120 m/s');
t = tiledlayout(fig2, 1, 1, 'Padding', 'none', 'TileSpacing', 'none');
set(fig2, 'Units','centimeters', ...
         'Position', [7,7,figwidth,figheight], ...
         'defaultAxesFontSize', docFontSize, ...
         'defaultTextFontSize', docFontSize, ...
         'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex');
set(fig2, 'PaperUnits', 'centimeters');
set(fig2, 'PaperSize', [figwidth figheight]);
set(fig2, 'PaperPosition', [0 0 figwidth figheight]);
ax = nexttile;
figureName = 'Fig09b_PoleMap_V120';
V_inf=120;
G = build_G_RectWing(V_inf,imuIDX,ailIDX);
EV=pole(G);

for i = 1:length(EV)
    plot(real(EV(i)), imag(EV(i)),"k.",'MarkerSize',stdMarkerSize);
    hold on
end
xline(0, 'k--', 'stability limit', 'Interpreter', 'latex');
xlabel('Real Part $\Re (\lambda)$', 'Interpreter', 'latex');
ylabel('Imaginary Part $\Im (\lambda)$', 'Interpreter', 'latex');
grid on;
grid minor;
ax.XLim = [-25, 5];
ax.YLim = [-200, 200];
ax.TickLabelInterpreter = 'latex';
set(ax, 'Color', 'w');
set(fig2, 'Color', 'w');

fullFigurePath = fullfile(figures_dir, figureName);
print(fig2, fullFigurePath, '-dpdf', '-vector');
% print(fig2, fullFigurePath, '-dmeta', '-vector');
