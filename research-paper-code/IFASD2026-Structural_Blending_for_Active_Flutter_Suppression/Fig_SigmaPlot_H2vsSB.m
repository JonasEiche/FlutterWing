% Sigma Plot Cont_H2 vs Cont_SB
%
% Singular-value plots at V_inf = 130 m/s: disturbance -> modal displacement,
% disturbance -> control-surface command, and the controllers themselves.
% Produces: Figures/Fig21_SigmaPlot_dist2qf_H2vsSB.pdf
%           Figures/Fig22_SigmaPlot_dist2phi_d_H2vsSB.pdf
%           Figures/Fig23_SigmaPlot_uz_ddot2phi_d_H2vsSB.pdf
%           Figures/FigXX_SigmaPlot_dist2qf_SBvsOL.pdf (not in the paper)
% Requires: Data/controller_imu18_ail18_structural_blending.mat (from R0)
% Recomputes the frequency responses from the controller file on every run.
% The Data/SigmaPlot_H2vsSB.mat written below is a by-product for inspection
% only: it is git-ignored and read by no script.

%% LOAD DATA
clearvars
textwidth = 15.98; % cm (IFASD 2026 template: 455.24pt)
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
data_dir = fullfile(script_dir, 'Data');

ContPath = fullfile(data_dir, 'controller_imu18_ail18_structural_blending.mat');
load(ContPath, 'cont_H2', 'cont_SB')



%  -5-6-7-8-
% |         |
%  -1-2-3-4-
num_modes = 5;     % must equal num_modes hard-coded in build_G_RectWing / build_P_RectWing
V_inf=130;
imuIDX = 1:8;
ailIDX = 1:8;
modesIDX = [1,2];
P = build_P_RectWing(V_inf,imuIDX,ailIDX,modesIDX);

CL_H2 = lft(P,cont_H2);
CL_SB = lft(P,cont_SB);

nm=length(modesIDX);
nym = length(imuIDX);
nud = length(ailIDX);

distIDX = 1:nm;
qfIDX = 1:nm;
qf_dotIDX = nm+(1:nm);
phi18_dIDX = 2*nm+(1:nud);


%% COMPUTE SIGMA VALUES

% DIST-->QF  (H_2 Blending, Structural Blending, Open Loop)
[sv_dist2qf_H2,wout_dist2qf_H2] = sigma(CL_H2(qfIDX,distIDX));
[sv_dist2qf_SB,wout_dist2qf_SB] = sigma(CL_SB(qfIDX,distIDX),wout_dist2qf_H2);
[sv_dist2qf_OL,wout_dist2qf_OL] = sigma(P(qfIDX,distIDX),wout_dist2qf_H2);
sv_dist2qf_H2_db = 20*log10(sv_dist2qf_H2);
sv_dist2qf_SB_db = 20*log10(sv_dist2qf_SB);
sv_dist2qf_OL_db = 20*log10(sv_dist2qf_OL);

% DIST-->PHI_D  (largest singular value)
[sv_dist2phi_d_H2,wout_dist2phi_d_H2] = sigma(CL_H2(phi18_dIDX,distIDX));
[sv_dist2phi_d_SB,wout_dist2phi_d_SB] = sigma(CL_SB(phi18_dIDX,distIDX),wout_dist2phi_d_H2);
sv_dist2phi_d_H2_db = 20*log10(sv_dist2phi_d_H2);
sv_dist2phi_d_SB_db = 20*log10(sv_dist2phi_d_SB);

sv_dist2phi_d_H2_db = sv_dist2phi_d_H2_db(1,:);
sv_dist2phi_d_SB_db = sv_dist2phi_d_SB_db(1,:);

% CONTROLLER SIGMA (u_z_ddot --> phi_d)
[sv_Cont_H2,wout_Cont_H2] = sigma(cont_H2);
[sv_Cont_SB,wout_Cont_SB] = sigma(cont_SB,wout_Cont_H2);
sv_Cont_H2_db = 20*log10(sv_Cont_H2);
sv_Cont_SB_db = 20*log10(sv_Cont_SB);

sv_Cont_H2_db = sv_Cont_H2_db(1,:);
sv_Cont_SB_db = sv_Cont_SB_db(1,:);


%% SAVE RESULTS TO MAT FILE
sigmaDataPath = fullfile(data_dir, 'SigmaPlot_H2vsSB.mat');
save(sigmaDataPath, ...
    'sv_dist2qf_H2', 'wout_dist2qf_H2', 'sv_dist2qf_H2_db', ...
    'sv_dist2qf_SB', 'wout_dist2qf_SB', 'sv_dist2qf_SB_db', ...
    'sv_dist2qf_OL', 'wout_dist2qf_OL', 'sv_dist2qf_OL_db', ...
    'sv_dist2phi_d_H2', 'wout_dist2phi_d_H2', 'sv_dist2phi_d_H2_db', ...
    'sv_dist2phi_d_SB', 'wout_dist2phi_d_SB', 'sv_dist2phi_d_SB_db', ...
    'sv_Cont_H2', 'wout_Cont_H2', 'sv_Cont_H2_db', ...
    'sv_Cont_SB', 'wout_Cont_SB', 'sv_Cont_SB_db');


%% DIST-->QF  '$H_2$ Blending','Modal Blending'  --------------------------

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
h1=semilogx(wout_dist2qf_H2, sv_dist2qf_H2_db, 'k-', 'LineWidth', stdLineWidth);  % Black line
hold on
h2=semilogx(wout_dist2qf_SB, sv_dist2qf_SB_db, 'k--', 'LineWidth', stdLineWidth);  % Black line

% h3=semilogx(wout_dist2qf_OL, sv_dist2qf_OL_db, 'k:', 'LineWidth', 1.5);  % Black line

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
xlim([1e0 1e3]);

lgd = legend([h1(1),h2(1)],{'$H_2$ Blending','Structural Blending'}, 'Location', 'southwest', 'FontSize', docFontSize);
% lgd = legend([h1(1),h2(1),h3(1)],{'$H_2$ Blending','Modal Blending', 'Open Loop'}, 'Location', 'southwest');
set(lgd, 'Interpreter', 'latex');

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
figureName = 'Fig21_SigmaPlot_dist2qf_H2vsSB';
fullFigurePath = fullfile(figures_dir, figureName);
print(f_dist2qf, fullFigurePath, '-dpdf', '-vector');
% print(f_dist2qf, fullFigurePath, '-dmeta', '-vector');


%% DIST-->QF    'Modal Blending', 'Open Loop' -----------------------------

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

h2=semilogx(wout_dist2qf_SB, sv_dist2qf_SB_db, 'k--', 'LineWidth', stdLineWidth);  % Black line
hold on
h3=semilogx(wout_dist2qf_OL, sv_dist2qf_OL_db, 'k:', 'LineWidth', stdLineWidth);  % Black line

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
xlim([1e0 1e3]);


lgd = legend([h2(1),h3(1)],{'Structural Blending', 'Open Loop'}, 'Location', 'southwest', 'FontSize', docFontSize);
set(lgd, 'Interpreter', 'latex');

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
figureName = 'FigXX_SigmaPlot_dist2qf_SBvsOL';
fullFigurePath = fullfile(figures_dir, figureName);
print(f_dist2qf_OL, fullFigurePath, '-dpdf', '-vector');
% print(f_dist2qf_OL, fullFigurePath, '-dmeta', '-vector');



%% distqf --> phi_d ------------------------------------------------------

f_dist2phi4_d = figure('Name',['Singular Values of Frequency Response dist2phi_d','   V_inf=',num2str(V_inf)]);
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
semilogx(wout_dist2phi_d_H2, sv_dist2phi_d_H2_db, 'k-',wout_dist2phi_d_SB, sv_dist2phi_d_SB_db, 'k--', 'LineWidth',stdLineWidth);  % Black line
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
xlim([1e0 1e3])
grid on;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';

lgd = legend({'$H_2$ Blending','Structural Blending'}, 'Location', 'southeast', 'FontSize', docFontSize);
set(lgd, 'Interpreter', 'latex');

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
figureName = 'Fig22_SigmaPlot_dist2phi_d_H2vsSB';
fullFigurePath = fullfile(figures_dir, figureName);
print(f_dist2phi4_d, fullFigurePath, '-dpdf', '-vector');
% print(f_dist2phi4_d, fullFigurePath, '-dmeta', '-vector');

%% u_z_ddot --> phi4 ------------------------------------------------------

f_uz_ddot2phi_d = figure('Name',['Singular Values of Frequency Response uz_ddot2phi_d','   V_inf=',num2str(V_inf)]);
% t = tiledlayout(f_dist2qf, 1, 1, 'Padding', 'none', 'TileSpacing', 'none');
set(f_uz_ddot2phi_d,'defaultTextInterpreter','latex');
set(f_uz_ddot2phi_d, 'Color', 'w'); % Set figure background to white
set(f_uz_ddot2phi_d, 'Units','centimeters', ...
         'Position', [7,7,figwidth,figheight], ...
         'defaultAxesFontSize', docFontSize, ...
         'defaultTextFontSize', docFontSize, ...
         'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex');
set(f_uz_ddot2phi_d, 'PaperUnits', 'centimeters');
set(f_uz_ddot2phi_d, 'PaperSize', [figwidth figheight]);        % Set PDF page size to match figure
set(f_uz_ddot2phi_d, 'PaperPosition', [0 0 figwidth figheight]); % Position plot to fill the PDF page exactly
% ax = nexttile;
% Plot magnitude
semilogx(wout_Cont_H2, sv_Cont_H2_db, 'k-',wout_Cont_SB, sv_Cont_SB_db, 'k--', 'LineWidth', stdLineWidth);  % Black line
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
xlim([1e0 1e3]);

lgd = legend({'$H_2$ Blending','Structural Blending'}, 'Location', 'southwest', 'FontSize', docFontSize);
set(lgd, 'Interpreter', 'latex');

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
figureName = 'Fig23_SigmaPlot_uz_ddot2phi_d_H2vsSB';
fullFigurePath = fullfile(figures_dir, figureName);
print(f_uz_ddot2phi_d, fullFigurePath, '-dpdf', '-vector');
% print(f_uz_ddot2phi4_d, fullFigurePath, '-dmeta', '-vector');
