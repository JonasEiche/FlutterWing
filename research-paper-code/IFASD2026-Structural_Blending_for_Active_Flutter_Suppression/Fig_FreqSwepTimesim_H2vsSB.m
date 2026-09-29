% Frequency Sweep Timesim Cont_H2 vs Cont_SB
%
% Time simulation of a 0-8 Hz chirp disturbance on the two flutter modes at
% V_inf = 130 m/s for the H2-blending and structural-blending closed loops.
% Produces: Figures/Fig24_FreqSwepTimesim_H2_SB_qf1.pdf
%           Figures/Fig25_FreqSwepTimesim_H2vsSB_phi4_d.pdf
% Requires: Data/controller_imu18_ail18_structural_blending.mat (from R0)
%           and the Signal Processing Toolbox (chirp).
% Recomputes the simulation from the controller file on every run. The
% Data/FreqSwepTimesim_H2vsSB_results.mat written below is a by-product for
% inspection only: it is git-ignored and read by no script.
clearvars
textwidth = 15.98; % cm (IFASD 2026 template: 455.24pt)
figwidth = 0.9*textwidth; % cm
figheight = 0.5*figwidth;


docFontSize = 9; % pt
stdLineWidth = 0.5; 


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
% 
% distIDX = 1:nm;
% qfIDX = 1:nm;
% qf_dotIDX = nm+[1:nm];
% phi4_dIDX = 2*nm+[1:nud];

%% distf1

dist_strength_qf1=0.3;
dist_strength_qf2=0.3;
freq0 = 0; % Hz
freq1 = 8;
t = 0:0.001:16; % s
nt = length(t);
u_chirp1 = chirp(t,freq0,t(end),freq1,'linear');    %(t_,freq@t0,t1,freq@t1)
u_chirp2 = chirp(t,freq0,t(end),freq1,'linear',-75);    %(t_,freq@t0,t1,freq@t1,method,initialphase)

IN_dist_q_f_ddot = [u_chirp1*dist_strength_qf1;
                    u_chirp2*dist_strength_qf2;
                    zeros(nm-2,nt)];

% IN_dist_q_f_ddot = [u_sine1*dist_strength_qf1;
%                     u_sine2*dist_strength_qf2;
%                     zeros(nm-2,L)];

% IN_dist_q_f_ddot = [u_square'*dist_strength_qf1;
%                     u_square'*dist_strength_qf2;
%                     zeros(nm-2,L)];

IN_noise_u_z_ddot = zeros(nym,nt);

IN_CL = [IN_dist_q_f_ddot;
         IN_noise_u_z_ddot];

%% LSIM
OUT_H2 = lsim(CL_H2,IN_CL',t');
OUT_SB = lsim(CL_SB,IN_CL',t');

OUT_H2_q_f1 = OUT_H2(:,1);
OUT_H2_q_f2 = OUT_H2(:,2);
OUT_H2_phi4_d = OUT_H2(:,2*nm+4);
OUT_H2_phi18_d = OUT_H2(:,2*nm+(1:nud));

OUT_SB_q_f1 = OUT_SB(:,1);
OUT_SB_q_f2 = OUT_SB(:,2);
OUT_SB_phi4_d = OUT_SB(:,2*nm+4);
OUT_SB_phi18_d = OUT_SB(:,2*nm+(1:nud));

%% Save LSIM results
lsimResultsPath = fullfile(data_dir, 'FreqSwepTimesim_H2vsSB_results.mat');
save(lsimResultsPath, 't', 'IN_CL', 'IN_dist_q_f_ddot', 'IN_noise_u_z_ddot', ...
    'OUT_H2', 'OUT_H2_q_f1', 'OUT_H2_q_f2', 'OUT_H2_phi4_d', 'OUT_H2_phi18_d', ...
    'OUT_SB', 'OUT_SB_q_f1', 'OUT_SB_q_f2', 'OUT_SB_phi4_d', 'OUT_SB_phi18_d', ...
    'V_inf', 'imuIDX', 'ailIDX', 'modesIDX', ...
    'dist_strength_qf1', 'dist_strength_qf2', 'freq0', 'freq1');

%% (1) distf1, qf1
fig_qf1 = figure('Name',['Closed Loop Disturbance Strength qf1:',num2str(dist_strength_qf1),'  / qf2:',num2str(dist_strength_qf2),'   V_inf=',num2str(V_inf)]);
% t = tiledlayout(fig_qf1, 1, 1, 'Padding', 'none', 'TileSpacing', 'none');
set(fig_qf1,'defaultTextInterpreter','latex');
set(fig_qf1, 'Color', 'w'); % Set figure background to white
set(fig_qf1, 'Units','centimeters', ...
         'Position', [7,7,figwidth,figheight], ...
         'defaultAxesFontSize', docFontSize, ...
         'defaultTextFontSize', docFontSize, ...
         'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex');
set(fig_qf1, 'PaperUnits', 'centimeters');
set(fig_qf1, 'PaperSize', [figwidth figheight]);        % Set PDF page size to match figure
set(fig_qf1, 'PaperPosition', [0 0 figwidth figheight]); % Position plot to fill the PDF page exactly
% ax = nexttile;
plot(t,OUT_H2_q_f1, 'k-',...
     t,OUT_SB_q_f1,'k--', ...
     'LineWidth',stdLineWidth)
xline(9, ':', '4.5 Hz', 'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'right', ...
    'LabelOrientation', 'horizontal', ...
    'LineWidth',1.5, ...
    'Interpreter', 'latex');

% title('Timesimulation of Frequency Sweep on Disturbance', 'Interpreter', 'latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', docFontSize) 
ylabel('Displacement of 1st Structural Mode', 'Interpreter', 'latex', 'FontSize', docFontSize) 
legend({'Bending Disp. $H_2$ Blending', 'Bending Disp. Structural Blending'},'Location','southwest', 'Interpreter', 'latex', 'FontSize', docFontSize);
grid on
set(fig_qf1, 'Color', 'w'); % Set figure background to white
ax = fig_qf1.CurrentAxes; % Get current axes
ax.TickLabelInterpreter = 'latex'; % Set tick labels to LaTeX
ax.FontSize = docFontSize;
set(ax, 'Color', 'w'); % Set axes background to white
script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
% figureName = 'FreqSwepTimesim_CL_H2_CL_SB_q_f1.pdf';
figureName = 'Fig24_FreqSwepTimesim_H2_SB_qf1';
fullFigurePath = fullfile(figures_dir, figureName);
print(fig_qf1, fullFigurePath, '-dpdf', '-vector');
% print(fig_qf1, fullFigurePath, '-dmeta', '-vector');
 

% %% (2) distf1, qf2
% fig_qf2 = figure('Name',['Closed Loop Disturbance Strength qf1:',num2str(dist_strength_qf1),'  / qf2:',num2str(dist_strength_qf2),'   V_inf=',num2str(V_inf)]);
% % t = tiledlayout(fig_qf2, 1, 1, 'Padding', 'none', 'TileSpacing', 'none');
% set(fig_qf2,'defaultTextInterpreter','latex');
% set(fig_qf2, 'Color', 'w'); % Set figure background to white
% set(fig_qf2, 'Units','centimeters', ...
%          'Position', [7,7,figwidth,figheight], ...
%          'defaultAxesFontSize', docFontSize, ...
%          'defaultTextFontSize', docFontSize, ...
%          'defaultTextInterpreter', 'latex', ...
%          'defaultAxesTickLabelInterpreter', 'latex', ...
%          'defaultLegendInterpreter', 'latex');
% set(fig_qf2, 'PaperUnits', 'centimeters');
% set(fig_qf2, 'PaperSize', [figwidth figheight]);        % Set PDF page size to match figure
% set(fig_qf2, 'PaperPosition', [0 0 figwidth figheight]); % Position plot to fill the PDF page exactly
% % ax = nexttile;
% plot(t,OUT_H2_q_f2, 'k-',...
%      t,OUT_SB_q_f2,'k--', ...
%      t,IN_dist_q_f_ddot(1,:),'k:','LineWidth',stdLineWidth)
% 
% % title('Timesimulation of Frequency Sweep on Disturbance', 'Interpreter', 'latex', 'FontSize', docFontSize)
% xlabel('Time, s', 'Interpreter', 'latex', 'FontSize', docFontSize) 
% ylabel('Displacement of 2nd Structural Mode', 'Interpreter', 'latex', 'FontSize', docFontSize) 
% legend({'Displacement of 2nd Structural Mode $H_2$ Blending', 'Displacement of 2nd Structural Mode Modal Blending','Disturbance on 1st Structural Mode'},'Location','southeast', 'Interpreter', 'latex', 'FontSize', docFontSize);
% grid on
% set(fig_qf2, 'Color', 'w'); % Set figure background to white
% ax = fig_qf2.CurrentAxes; % Get current axes
% ax.TickLabelInterpreter = 'latex'; % Set tick labels to LaTeX
% set(ax, 'Color', 'w'); % Set axes background to white
% script_path = mfilename('fullpath');
% script_dir = fileparts(script_path);
% figures_dir = fullfile(script_dir, 'Figures');
% if ~exist(figures_dir, 'dir')
%     mkdir(figures_dir);
% end
% figureName = 'FreqSwepTimesim_CL_V130_8imu_rCflut_LPF.pdf';
% fullFigurePath = fullfile(figures_dir, figureName);
% print(fig_qf2, fullFigurePath, '-dpdf', '-vector');

%% (3) distf1, phi4d
fig_phi4_d = figure('Name',['Closed Loop Disturbance Strength qf1:',num2str(dist_strength_qf1),'  / qf2:',num2str(dist_strength_qf2),'   V_inf=',num2str(V_inf)]);
% t = tiledlayout(fig_phi4_d, 1, 1, 'Padding', 'none', 'TileSpacing', 'none');
set(fig_phi4_d,'defaultTextInterpreter','latex');
set(fig_phi4_d, 'Color', 'w'); % Set figure background to white
set(fig_phi4_d, 'Units','centimeters', ...
         'Position', [7,7,figwidth,figheight], ...
         'defaultAxesFontSize', docFontSize, ...
         'defaultTextFontSize', docFontSize, ...
         'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex');
set(fig_phi4_d, 'PaperUnits', 'centimeters');
set(fig_phi4_d, 'PaperSize', [figwidth figheight]);        % Set PDF page size to match figure
set(fig_phi4_d, 'PaperPosition', [0 0 figwidth figheight]); % Position plot to fill the PDF page exactly
% ax = nexttile;
% plot(t,OUT_H2_phi4_d, 'k-',...
%      t,OUT_SB_phi4_d,'k--', ...
%      t,IN_dist_q_f_ddot(1,:),'k:')

plot(t,OUT_H2_phi4_d, 'k-',...
     t,OUT_SB_phi4_d,'k--','LineWidth',stdLineWidth)



% % Plot all 8 channels of OUT_H2_phi18_d in blue with solid lines
% h1 = plot(t, OUT_H2_phi18_d(:,1), 'b-', 'LineWidth', stdLineWidth);
% plot(t, OUT_H2_phi18_d(:,2:end), 'b-', 'LineWidth', stdLineWidth, 'HandleVisibility', 'off');
% hold on;
% % Plot all 8 channels of OUT_SB_phi18_d in red with dashed lines
% h2 = plot(t, OUT_SB_phi18_d(:,1), 'r--', 'LineWidth', stdLineWidth);
% plot(t, OUT_SB_phi18_d(:,2:end), 'r--', 'LineWidth', stdLineWidth, 'HandleVisibility', 'off');

xline(9, ':', '4.5 Hz','LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'right', ...
    'LabelOrientation', 'horizontal', ... 
    'LineWidth',1.5, ...
    'Interpreter', 'latex');
% title('Timesimulation of Frequency Sweep on Disturbance', 'Interpreter', 'latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', docFontSize) 
ylabel('Outboard Aileron Command [rad] ', 'Interpreter', 'latex', 'FontSize', docFontSize) 
legend({'Aileron Command $H_2$ Blending','Aileron Command Structural Blending'},'Location','southwest', 'Interpreter', 'latex', 'FontSize', docFontSize);
grid on
set(fig_phi4_d, 'Color', 'w'); % Set figure background to white
ax = fig_phi4_d.CurrentAxes; % Get current axes
ax.TickLabelInterpreter = 'latex'; % Set tick labels to LaTeX
ax.FontSize = docFontSize;
set(ax, 'Color', 'w'); % Set axes background to white
script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
figureName = 'Fig25_FreqSwepTimesim_H2vsSB_phi4_d';
fullFigurePath = fullfile(figures_dir, figureName);
print(fig_phi4_d, fullFigurePath, '-dpdf', '-vector');
% print(fig_phi4_d, fullFigurePath, '-dmeta', '-vector');