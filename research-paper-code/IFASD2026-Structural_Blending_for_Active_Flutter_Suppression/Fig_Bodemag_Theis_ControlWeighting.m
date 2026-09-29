% Bode magnitude of control weighting W_u(s) used in H-infinity synthesis
clearvars
textwidth = 15.98; % cm (IFASD 2026 template: 455.24pt)
figwidth = 0.7*textwidth; % cm
figheight = 0.5*figwidth;
docFontSize = 9; % pt
stdLineWidth = 1.2;

% Control weighting parameters (must match R0_afs_structural_blending_imu18_ail18_synthesis.m)
w1 = 12;
w2 = 64;
s = tf('s');
W_theis=((s+w1)*(s+w2))/((s+0.01*w1)*(0.01*s+w2));
invW_theis = ((s+0.01*w1)*(0.01*s+w2))/((s+w1)*(s+w2));

% Frequency vector (log scale)
w = logspace(0, 3, 1000);  % From 1 to 1000 rad/s

% Frequency response
[mag, ~] = bode(W_theis, w);
mag = squeeze(mag);
mag_db = 20*log10(mag);

fig = figure('Name', 'Bode Magnitude of Control Weighting');
t = tiledlayout(fig, 1, 1, 'Padding', 'none', 'TileSpacing', 'none');
set(fig, 'Units','centimeters', ...
         'Position', [7,7,figwidth,figheight], ...
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
semilogx(w, mag_db, 'k', 'LineWidth', stdLineWidth);
hold on;
xline(28, '--k', 'flutter frequency','Interpreter', 'latex','LabelVerticalAlignment', 'bottom');
ax.XScale = 'log';
ax.YLabel.Interpreter = 'latex';
ax.XLabel.Interpreter = 'latex';
ax.TickLabelInterpreter = 'latex';
ax.FontSize = docFontSize;
ax.XLabel.String = '$\omega$ [rad/s]';
ax.YLabel.String = 'Magnitude [dB]';
ax.Box = 'on';
set(ax, 'Color', 'w');

ylim([-5, 30]);
yticks(0:5:25);

grid on;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
set(fig, 'Color', 'w');
set(ax, 'Color', 'w');

figureName = 'Fig18_Bodemag_Theis_ControlWeighting';
script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
fullFigurePath = fullfile(figures_dir, figureName);
print(fig, fullFigurePath, '-dpdf', '-vector');
% print(fig, fullFigurePath, '-dmeta', '-vector');
