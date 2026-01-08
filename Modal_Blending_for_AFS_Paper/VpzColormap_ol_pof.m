% VpzColormap_ol_pof

clearvars
textwidth = 16.5730; % cm
figwidthVpz = 0.7*textwidth; % cm
figheightVpz = 1.0*figwidthVpz;

docFontSize = 9; % pt
stdLineWidth = 1.2; 


script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
data_dir = fullfile(script_dir, 'Data');
ContPath = fullfile(data_dir, 'Controller_imu48_ail4_pof.mat');
load(ContPath)

Cont1 = Cont_V130_2imu_POF;

%  -5-6-7-8-
% |         |
%  -1-2-3-4-
num_modes = 5;     % ACHTUNG ! must be equal to setting in defineRectWing !
imuIDX = [4,8];
ailIDX = 4;
modesIDX = [1,2];
nimu = length(imuIDX);
nail = length(ailIDX);
nmod = length(modesIDX);

V_inf=linspace(20,160,32);
nv = length(V_inf);
G = build_G_RectWing(V_inf,imuIDX,ailIDX);

CL1 = feedback(G,-Cont1);

%% Vg Plot  //  Vpzmap
[EV_OL,MS_OL] = getEigenvalueModeshape(G,num_modes);
[EV_Cont1,MS_Cont1] = getEigenvalueModeshape(CL1,num_modes);

num_modes_disp=2;
EV_Cont1_disp = EV_Cont1(1:2*num_modes_disp,:);
EV_OL_disp = EV_OL(1:2*num_modes_disp,:);

legendlist={"Open Loop", "Acceleration Feedback"};
EV_list={EV_OL_disp,EV_Cont1_disp};
vel_list={V_inf,V_inf};
vel_=vel_list; 
EV_=EV_list;

fig = figure('Name','VpzColormap OL');
set(fig,'defaultTextInterpreter','latex');
set(fig, 'Color', 'w'); % Set figure background to white
set(fig, 'Units','centimeters', ...
         'Position', [7,7,figwidthVpz,figheightVpz], ...
         'defaultAxesFontSize', docFontSize, ...
         'defaultTextFontSize', docFontSize, ...
         'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex');
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperSize', [figwidthVpz figheightVpz]);        % Set PDF page size to match figure
set(fig, 'PaperPosition', [0 0 figwidthVpz figheightVpz]); % Position plot to fill the PDF page exactly

for i=1:4
    EV_OL_disp_i=EV_OL_disp(i,:);
    
    surface([real(EV_OL_disp_i); real(EV_OL_disp_i)], [imag(EV_OL_disp_i); imag(EV_OL_disp_i)], [V_inf; V_inf], [V_inf; V_inf], ...
            'FaceColor', 'none', ...
            'EdgeColor', 'interp', ...
            'LineWidth', 3);
    
end
colormap turbo 
cb = colorbar('TickLabelInterpreter', 'latex');
cb.Label.Interpreter = 'latex';
cb.Label.String = '$V_\infty$ [m/s]';
cb.Label.FontSize = docFontSize;
cb.FontSize = docFontSize;

title(' Eigenvalue Loci', 'FontSize', docFontSize); 
xlabel('Real Part $\Re (\lambda)$', 'FontSize', docFontSize); 
ylabel('Imaginary Part $\Im (\lambda)$', 'FontSize', docFontSize); 
grid on; 
axis equal; box on
%yaxis([-20 20 -35 35])       % Set limits
ylim([-40 40])
ax = fig.CurrentAxes; % Get current axes
ax.TickLabelInterpreter = 'latex'; % Set tick labels to LaTeX
set(ax, 'Color', 'w'); % Set axes background to white
set(fig, 'Color', 'w'); % Set figure background to white

Vpzmap_PlotName = 'Vpzmap_OL_Colormap';
fullVpzmap_PlotPath = fullfile(figures_dir, Vpzmap_PlotName);
print(fig, fullVpzmap_PlotPath, '-dpdf', '-vector');
print(fig, fullVpzmap_PlotPath, '-dmeta', '-vector');

