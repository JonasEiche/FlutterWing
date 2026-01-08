clearvars
textwidth = 16.5730; % cm
figwidth = 0.7*textwidth; % cm
figheight = 0.5*figwidth;
docFontSize = 9; % pt
stdLineWidth = 1.2; 

w1 = 12;
w2 = 64;
s = tf('s');
W_theis=((s+w1)*(s+w2))/((s+0.01*w1)*(0.01*s+w2));
invW_theis = ((s+0.01*w1)*(0.01*s+w2))/((s+w1)*(s+w2));
BUTTER=W_theis;

% Frequency vector (log scale)
w = logspace(0, 3, 1000);  % From 1 to 1000 rad/s

% Frequency response
[mag, phase] = bode(BUTTER, w);
mag = squeeze(mag);  % Remove singleton dimensions

% Convert magnitude to dB
mag_db = 20*log10(mag);




fig = figure('Name', 'Bode Magnitude of BPF');
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
set(fig, 'PaperSize', [figwidth figheight]);        % Set PDF page size to match figure
set(fig, 'PaperPosition', [0 0 figwidth figheight]); % Position plot to fill the PDF page exactly
ax = nexttile;
% Plot magnitude
semilogx(w, mag_db, 'k', 'LineWidth', stdLineWidth);  % Black line
hold on;
% Add vertical line with LaTeX label
xline(28, '--k', 'flutter frequency','Interpreter', 'latex','LabelVerticalAlignment', 'bottom');
% ax = gca;
ax.XScale = 'log';
ax.YLabel.Interpreter = 'latex';
ax.XLabel.Interpreter = 'latex';
ax.TickLabelInterpreter = 'latex';
ax.FontSize = docFontSize;
ax.XLabel.String = '$\omega$ [rad/s]';
ax.YLabel.String = 'Magnitude [dB]';
ax.Box = 'on';
set(ax, 'Color', 'w'); % Set axes background to white
% Y-axis limits and ticks
ylim([-50, 10]);
% yticks(-40:10:0);

ylim([-5, 30]);
yticks(0:5:25);

% ylim([-30, 10]);
% yticks(-20:10:0);

grid on;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
set(fig, 'Color', 'w'); % Set figure background to white
set(ax, 'Color', 'w'); % Set axes background to white

figureName = 'bodemag_BPF';
script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
fullFigurePath = fullfile(figures_dir, figureName);
print(fig, fullFigurePath, '-dpdf', '-vector');
print(fig, fullFigurePath, '-dmeta', '-vector');


