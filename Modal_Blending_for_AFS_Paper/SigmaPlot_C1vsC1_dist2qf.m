% Sigma Plot Cont1 vs Cont2

%% LOAD DATA
clearvars
textwidth = 16.5730; % cm
figwidth = 0.9*textwidth; % cm
figheight = 0.5*figwidth;


docFontSize = 9; % pt
stdLineWidth = 1.2; 


script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
data_dir = fullfile(script_dir, 'Data');
ContPath = fullfile(data_dir, 'Controller_imu1_8_ail4_advanced_ReqMarg645.mat');
load(ContPath)
Cont1 = Cont_minQF_H2;
Cont2 = Cont_minQF_MB;



%  -5-6-7-8-
% |         |
%  -1-2-3-4-
num_modes = 5;     % ACHTUNG ! must be equal to setting in defineRectWing !
V_inf=130;
imuIDX = 1:8;
ailIDX = 4;
modesIDX = [1,2];
P = build_P_RectWing(V_inf,imuIDX,ailIDX,modesIDX);

CL1 = lft(P,Cont1);
CL2 = lft(P,Cont2);

nm=length(modesIDX);
nym = length(imuIDX);
nud = length(ailIDX);

distIDX = 1:nm;
qfIDX = 1:nm;
qf_dotIDX = nm+(1:nm);
phi4_dIDX = 2*nm+(1:nud);


%% DIST-->QF  '$H_2$ Blending','Modal Blending'  --------------------------

[sv1,wout1] = sigma(CL1(qfIDX,distIDX));
[sv2,wout2] = sigma(CL2(qfIDX,distIDX),wout1);
% [sv3,wout3] = sigma(P(qfIDX,distIDX),wout1);
sv1_db = 20*log10(sv1);
sv2_db = 20*log10(sv2);
% sv3_db = 20*log10(sv3);

f_dist2qf = figure('Name',['Singular Values of Frequency Response dist2qf','   V_inf=',num2str(V_inf)]);
% t = tiledlayout(f_dist2qf, 1, 1, 'Padding', 'none', 'TileSpacing', 'none');
set(f_dist2qf,'defaultTextInterpreter','latex');
set(f_dist2qf, 'Color', 'w'); % Set figure background to white
set(f_dist2qf, 'Units','centimeters', ...
         'Position', [7,7,figwidth,figheight], ...
         'defaultAxesFontSize', docFontSize, ...
         'defaultTextFontSize', docFontSize, ...
         'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex');
set(f_dist2qf, 'PaperUnits', 'centimeters');
set(f_dist2qf, 'PaperSize', [figwidth figheight]);        % Set PDF page size to match figure
set(f_dist2qf, 'PaperPosition', [0 0 figwidth figheight]); % Position plot to fill the PDF page exactly
% ax = nexttile;
% Plot magnitude
h1=semilogx(wout1, sv1_db, 'k-', 'LineWidth', stdLineWidth);  % Black line
hold on
h2=semilogx(wout2, sv2_db, 'k--', 'LineWidth', stdLineWidth);  % Black line

% h3=semilogx(wout3, sv3_db, 'k:', 'LineWidth', 1.5);  % Black line

xline(28, '--k', 'flutter frequency','Interpreter', 'latex','LabelVerticalAlignment', 'bottom', 'FontSize', docFontSize);
grid on
ax = gca;
ax.XScale = 'log';
ax.YLabel.Interpreter = 'latex';
ax.XLabel.Interpreter = 'latex';
ax.TickLabelInterpreter = 'latex';
ax.FontSize = docFontSize;
ax.XLabel.String = '$Frequency$ [rad/s]';
ax.YLabel.String = 'Singular Values [dB]';
ax.Box = 'on';
set(ax, 'Color', 'w'); % Set axes background to white
% Y-axis limits and ticks
% ylim([-50, 10]);
% yticks(-40:10:0);
grid on;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';

lgd = legend([h1(1),h2(1)],{'$H_2$ Blending','Modal Blending'}, 'Location', 'southwest', 'FontSize', docFontSize);
% lgd = legend([h1(1),h2(1),h3(1)],{'$H_2$ Blending','Modal Blending', 'Open Loop'}, 'Location', 'southwest');
set(lgd, 'Interpreter', 'latex');

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
figureName = 'sigmaplot_dist2qf_C1vsC2';
fullFigurePath = fullfile(figures_dir, figureName);
print(f_dist2qf, fullFigurePath, '-dpdf', '-vector');
print(f_dist2qf, fullFigurePath, '-dmeta', '-vector');


%% DIST-->QF    'Modal Blending', 'Open Loop' -----------------------------

% [sv1,wout1] = sigma(CL1(qfIDX,distIDX));
[sv2,wout2] = sigma(CL2(qfIDX,distIDX),wout1);
[sv3,wout3] = sigma(P(qfIDX,distIDX),wout1);
% sv1_db = 20*log10(sv1);
sv2_db = 20*log10(sv2);
sv3_db = 20*log10(sv3);

f_dist2qf_OL = figure('Name',['Singular Values of Frequency Response dist2qf','   V_inf=',num2str(V_inf)]);
% t = tiledlayout(f_dist2qf_OL, 1, 1, 'Padding', 'none', 'TileSpacing', 'none');
set(f_dist2qf_OL,'defaultTextInterpreter','latex');
set(f_dist2qf_OL, 'Color', 'w'); % Set figure background to white
set(f_dist2qf_OL, 'Units','centimeters', ...
         'Position', [7,7,figwidth,figheight], ...
         'defaultAxesFontSize', docFontSize, ...
         'defaultTextFontSize', docFontSize, ...
         'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex');
set(f_dist2qf_OL, 'PaperUnits', 'centimeters');
set(f_dist2qf_OL, 'PaperSize', [figwidth figheight]);        % Set PDF page size to match figure
set(f_dist2qf_OL, 'PaperPosition', [0 0 figwidth figheight]); % Position plot to fill the PDF page exactly
% ax = nexttile;

h2=semilogx(wout2, sv2_db, 'k--', 'LineWidth', stdLineWidth);  % Black line
hold on
h3=semilogx(wout3, sv3_db, 'k:', 'LineWidth', stdLineWidth);  % Black line

xline(28, '--k', 'flutter frequency','Interpreter', 'latex','LabelVerticalAlignment', 'bottom', 'FontSize', docFontSize);
grid on
ax = gca;
ax.XScale = 'log';
ax.YLabel.Interpreter = 'latex';
ax.XLabel.Interpreter = 'latex';
ax.TickLabelInterpreter = 'latex';
ax.FontSize = docFontSize;
ax.XLabel.String = 'Frequency [rad/s]';
ax.YLabel.String = 'Singular Values [dB]';
ax.Box = 'on';
set(ax, 'Color', 'w'); % Set axes background to white
% Y-axis limits and ticks
% ylim([-50, 10]);
% yticks(-40:10:0);
grid on;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';


lgd = legend([h2(1),h3(1)],{'Modal Blending', 'Open Loop'}, 'Location', 'southwest', 'FontSize', docFontSize);
set(lgd, 'Interpreter', 'latex');

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
figureName = 'sigmaplot_dist2qf_C2vsOL';
fullFigurePath = fullfile(figures_dir, figureName);
print(f_dist2qf_OL, fullFigurePath, '-dpdf', '-vector');
print(f_dist2qf_OL, fullFigurePath, '-dmeta', '-vector');



%% distqf --> phi4_d ------------------------------------------------------

[sv1,wout1] = sigma(CL1(phi4_dIDX,distIDX));
[sv2,wout2] = sigma(CL2(phi4_dIDX,distIDX),wout1);
% [sv1,wout1] = sigma(CL1(phi4_dIDX,distIDX(1)));
% [sv2,wout2] = sigma(CL2(phi4_dIDX,distIDX(1)),wout1);
sv1_db = 20*log10(sv1);
sv2_db = 20*log10(sv2);

f_dist2phi4_d = figure('Name',['Singular Values of Frequency Response dist2phi4_d','   V_inf=',num2str(V_inf)]);
% t = tiledlayout(f_dist2qf, 1, 1, 'Padding', 'none', 'TileSpacing', 'none');
set(f_dist2phi4_d,'defaultTextInterpreter','latex');
set(f_dist2phi4_d, 'Color', 'w'); % Set figure background to white
set(f_dist2phi4_d, 'Units','centimeters', ...
         'Position', [7,7,figwidth,figheight], ...
         'defaultAxesFontSize', docFontSize, ...
         'defaultTextFontSize', docFontSize, ...
         'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex');
set(f_dist2phi4_d, 'PaperUnits', 'centimeters');
set(f_dist2phi4_d, 'PaperSize', [figwidth figheight]);        % Set PDF page size to match figure
set(f_dist2phi4_d, 'PaperPosition', [0 0 figwidth figheight]); % Position plot to fill the PDF page exactly
% ax = nexttile;
% Plot magnitude
semilogx(wout1, sv1_db, 'k-',wout2, sv2_db, 'k--', 'LineWidth',stdLineWidth);  % Black line
xline(28, '--k', 'flutter frequency','Interpreter', 'latex','LabelVerticalAlignment', 'bottom', 'FontSize', docFontSize);
grid on
ax = gca;
ax.XScale = 'log';
ax.YLabel.Interpreter = 'latex';
ax.XLabel.Interpreter = 'latex';
ax.TickLabelInterpreter = 'latex';
ax.FontSize = docFontSize;
ax.XLabel.String = 'Frequency [rad/s]';
ax.YLabel.String = 'Singular Values [dB]';
ax.Box = 'on';
set(ax, 'Color', 'w'); % Set axes background to white
% Y-axis limits and ticks
% ylim([-50, 10]);
% yticks(-40:10:0);
xlim([1e1 1e2])
grid on;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';

lgd = legend({'$H_2$ Blending','Modal Blending'}, 'Location', 'southeast', 'FontSize', docFontSize);
set(lgd, 'Interpreter', 'latex');

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
figureName = 'sigmaplot_dist2phi4_d_CL1vsCL2';
fullFigurePath = fullfile(figures_dir, figureName);
print(f_dist2phi4_d, fullFigurePath, '-dpdf', '-vector');
print(f_dist2phi4_d, fullFigurePath, '-dmeta', '-vector');

%% u_z_ddot --> phi4 ------------------------------------------------------
 
[sv1,wout1] = sigma(Cont1);
[sv2,wout2] = sigma(Cont2,wout1);
sv1_db = 20*log10(sv1);
sv2_db = 20*log10(sv2);

f_uz_ddot2phi4_d = figure('Name',['Singular Values of Frequency Response uz_ddot2phi4_d','   V_inf=',num2str(V_inf)]);
% t = tiledlayout(f_dist2qf, 1, 1, 'Padding', 'none', 'TileSpacing', 'none');
set(f_uz_ddot2phi4_d,'defaultTextInterpreter','latex');
set(f_uz_ddot2phi4_d, 'Color', 'w'); % Set figure background to white
set(f_uz_ddot2phi4_d, 'Units','centimeters', ...
         'Position', [7,7,figwidth,figheight], ...
         'defaultAxesFontSize', docFontSize, ...
         'defaultTextFontSize', docFontSize, ...
         'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex');
set(f_uz_ddot2phi4_d, 'PaperUnits', 'centimeters');
set(f_uz_ddot2phi4_d, 'PaperSize', [figwidth figheight]);        % Set PDF page size to match figure
set(f_uz_ddot2phi4_d, 'PaperPosition', [0 0 figwidth figheight]); % Position plot to fill the PDF page exactly
% ax = nexttile;
% Plot magnitude
semilogx(wout1, sv1_db, 'k-',wout2, sv2_db, 'k--', 'LineWidth', stdLineWidth);  % Black line
xline(28, '--k', 'flutter frequency','Interpreter', 'latex','LabelVerticalAlignment', 'bottom', 'FontSize', docFontSize);
grid on
ax = gca;
ax.XScale = 'log';
ax.YLabel.Interpreter = 'latex';
ax.XLabel.Interpreter = 'latex';
ax.TickLabelInterpreter = 'latex';
ax.FontSize = docFontSize;
ax.XLabel.String = 'Frequency [rad/s]';
ax.YLabel.String = 'Singular Values [dB]';
ax.Box = 'on';
set(ax, 'Color', 'w'); % Set axes background to white
% Y-axis limits and ticks
% ylim([-50, 10]);
% yticks(-40:10:0);
grid on;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';

lgd = legend({'$H_2$ Blending','Modal Blending'}, 'Location', 'northwest', 'FontSize', docFontSize);
set(lgd, 'Interpreter', 'latex');

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
figureName = 'sigmaplot_uz_ddot2phi4_C1vsC2';
fullFigurePath = fullfile(figures_dir, figureName);
print(f_uz_ddot2phi4_d, fullFigurePath, '-dpdf', '-vector');
print(f_uz_ddot2phi4_d, fullFigurePath, '-dmeta', '-vector');