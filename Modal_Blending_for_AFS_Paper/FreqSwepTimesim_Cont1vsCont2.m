% Frequency Sweep Timesim Cont1 vs Cont2
clearvars
textwidth = 16.5730; % cm
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


IN_noise_u_z_ddot = zeros(nym,nt);

IN_CL = [IN_dist_q_f_ddot;
         IN_noise_u_z_ddot];

%% LSIM
OUT1 = lsim(CL1,IN_CL',t');
OUT2 = lsim(CL2,IN_CL',t');

OUT1_q_f1 = OUT1(:,1);
OUT1_q_f2 = OUT1(:,2);
OUT1_phi4_d = OUT1(:,2*nm+1);

OUT2_q_f1 = OUT2(:,1);
OUT2_q_f2 = OUT2(:,2);
OUT2_phi4_d = OUT2(:,2*nm+1);

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
plot(t,OUT1_q_f1, 'k-',...
     t,OUT2_q_f1,'k--', ...
     'LineWidth',stdLineWidth)
xline(9, ':', '4.5 Hz', 'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'right', ...
    'LabelOrientation', 'horizontal', ...
    'LineWidth',1.5, ...
    'Interpreter', 'latex');

% title('Timesimulation of Frequency Sweep on Disturbance', 'Interpreter', 'latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', docFontSize) 
ylabel('Displacement of 1st Structural Mode', 'Interpreter', 'latex', 'FontSize', docFontSize) 
legend({'Bending Disp. $H_2$ Blending', 'Bending Disp. Modal Blending'},'Location','southwest', 'Interpreter', 'latex', 'FontSize', docFontSize);
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
figureName = 'FreqSwepTimesim_CL1_CL2_q_f1';
fullFigurePath = fullfile(figures_dir, figureName);

print(fig_qf1, fullFigurePath, '-dpdf', '-vector');
print(fig_qf1, fullFigurePath, '-dmeta', '-vector');


%% (2) distf1, phi4d
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
% plot(t,OUT1_phi4_d, 'k-',...
%      t,OUT2_phi4_d,'k--', ...
%      t,IN_dist_q_f_ddot(1,:),'k:')

plot(t,OUT1_phi4_d, 'k-',...
     t,OUT2_phi4_d,'k--','LineWidth',stdLineWidth)
xline(9, ':', '4.5 Hz','LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'right', ...
    'LabelOrientation', 'horizontal', ... 
    'LineWidth',1.5, ...
    'Interpreter', 'latex');
% title('Timesimulation of Frequency Sweep on Disturbance', 'Interpreter', 'latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', docFontSize) 
ylabel('Aileron Command [rad] ', 'Interpreter', 'latex', 'FontSize', docFontSize) 
legend({'Aileron Command $H_2$ Blending','Aileron Command Modal Blending'},'Location','southwest', 'Interpreter', 'latex', 'FontSize', docFontSize);
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
figureName = 'FreqSwepTimesim_CL1_CL2_phi4_d';
fullFigurePath = fullfile(figures_dir, figureName);

print(fig_phi4_d, fullFigurePath, '-dpdf', '-vector');
print(fig_phi4_d, fullFigurePath, '-dmeta', '-vector');