% Vg_Plot_and_Vpzmap_CL1vsCL2vsOL

clearvars
textwidth = 16.5730; % cm
figwidthVg = 0.9*textwidth; % cm
figheightVg = 0.7*figwidthVg;

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
ContPath = fullfile(data_dir, 'Controller_imu1_8_ail4_advanced_ReqMarg645.mat');
load(ContPath)


Cont1 = Cont_minQF_H2;
Cont2 = Cont_minQF_MB;

%  -5-6-7-8-
% |         |
%  -1-2-3-4-
num_modes = 5;     % ACHTUNG ! must be equal to setting in defineRectWing !
imuIDX = 1:8;
ailIDX = 4;
modesIDX = [1,2];
nimu = length(imuIDX);
nail = length(ailIDX);
nmod = length(modesIDX);

V_inf=linspace(20,160,32);
nv = length(V_inf);
G = build_G_RectWing(V_inf,imuIDX,ailIDX);

CL1 = feedback(G,-Cont1);
CL2 = feedback(G,-Cont2);

%% Vg Plot  //  Vpzmap
[EV_OL,MS_OL] = getEigenvalueModeshape(G,num_modes);
[EV_Cont1,MS_Cont1] = getEigenvalueModeshape(CL1,num_modes);
[EV_Cont2,MS_Cont2] = getEigenvalueModeshape(CL2,num_modes);

num_modes_disp=2;
EV_Cont1_disp = EV_Cont1(1:2*num_modes_disp,:);
EV_Cont2_disp = EV_Cont2(1:2*num_modes_disp,:);
EV_OL_disp = EV_OL(1:2*num_modes_disp,:);

legendlist={"Open Loop", "$H_2$ Optimal Blending", "Modal Blending"};
EV_list={EV_OL_disp,EV_Cont1_disp,EV_Cont2_disp};
vel_list={V_inf,V_inf,V_inf};
vel_=vel_list; 
EV_=EV_list;

%% [f_Vg] = Vg_plot_mult_bw(vel_list, EV_list, legendlist);
% Define color and linestyles
black = [0 0 0];
linestyles = {'-', '--', ':', '-.'};  % add more if needed

if iscell(vel_)
    num_cases = length(vel_);
    assert(iscell(EV_) & length(EV_)==length(vel_),'Inconsistent input');
else
    num_cases=1;
    vel_={vel_};
    EV_={EV_};
end

FLUT_F = cell(1,num_cases);
FLUT_V = cell(1,num_cases);
DIV_V = cell(1,num_cases);
F = cell(1,num_cases);
D = cell(1,num_cases);

for i = 1:num_cases
    vel = vel_{i};
    EV = EV_{i};

    tol = 0.1;
    divtol = 1.0;

    crit_ind = cumsum((real(EV) > tol),2) == 1;
    flut_ind = (cumsum((real(EV) > tol),2) == 1 & imag(EV)>divtol);
    div_ind = (cumsum((real(EV) > tol),2) == 1 & imag(EV)<divtol & imag(EV)>0);

    flut_f = abs(EV(flut_ind))./(2*pi);
    VEL = repmat(vel,size(EV,1),1);
    flut_v = VEL(flut_ind);
    div_v = VEL(div_ind);

    for j = 1:length(flut_f)
        disp([num2str(flut_v(j)) ' m/s = flutter     '  num2str(flut_f(j)) ' Hz'])
    end
    for j = 1:length(div_v)
        disp([num2str(div_v(j)) ' m/s = divergence'])
    end
    if sum(sum(crit_ind)) == 0
        disp('No flutter or divergence')
    end
    f = abs(EV)./(2*pi);
    d = ( -real(EV)./abs(EV) ).*100;

    FLUT_F{i} = flut_f;
    FLUT_V{i} = flut_v;
    DIV_V{i} = div_v;
    F{i} = f;
    D{i} = d;
end

%% Plot Frequency over Speed
fig = figure('Name','Vg Plot Comparison');
set(fig,'defaultTextInterpreter','latex');
set(fig, 'Color', 'w'); % Set figure background to white
set(fig, 'Units','centimeters', ...
         'Position', [7,7,figwidthVg,figheightVg], ...
         'defaultAxesFontSize', docFontSize, ...
         'defaultTextFontSize', docFontSize, ...
         'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex');
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperSize', [figwidthVg figheightVg]);        % Set PDF page size to match figure
set(fig, 'PaperPosition', [0 0 figwidthVg figheightVg]); % Position plot to fill the PDF page exactly
subplot(211)
for i = 1:num_cases
    ls = linestyles{mod(i-1,length(linestyles))+1};  % cycle through line styles
    plot(vel_{i},F{i},'Color',black,'LineStyle',ls,'LineWidth',stdLineWidth)
    hold on
    plot(DIV_V{i},zeros(size(DIV_V{i})),'Ob', 'MarkerSize', 4, 'LineWidth', 2)  % divergence
    hold on
    plot(FLUT_V{i},FLUT_F{i},'Or', 'MarkerSize', 4, 'LineWidth', 2)             % flutter
    hold on
end
title('Frequency and Damping vs. Velocity', 'Interpreter', 'latex', 'FontSize', docFontSize); 
ylabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontSize', docFontSize);
xlim([vel(1) vel(end)]);
grid on
grid minor
ax = fig.CurrentAxes;
ax.TickLabelInterpreter = 'latex';
set(ax, 'Color', 'w');

%% Plot Damping over Speed
subplot(212)
x0 = [vel(1) vel(end)]; y0 = [0 0];

firstlines = zeros(1,num_cases);
for i = 1:num_cases
    ls = linestyles{mod(i-1,length(linestyles))+1};
    p{i} = plot(vel_{i},D{i},'Color',black,'LineStyle',ls,'LineWidth',stdLineWidth);
    hold on
    plot(DIV_V{i},-100*ones(size(DIV_V{i})),'Ob', 'MarkerSize', 4, 'LineWidth', 2)
    hold on
    plot(FLUT_V{i},zeros(size(FLUT_V{i})),'Or', 'MarkerSize', 4, 'LineWidth', 2)
    hold on
    firstlines(i) = p{i}(1);
end
plot(x0,y0,'--k')
lgd = legend(firstlines,legendlist,'Location','southwest', 'Interpreter', 'latex', 'FontSize', docFontSize);
lgd.AutoUpdate = 'off';
xlabel('Velocity (m/s)', 'Interpreter', 'latex', 'FontSize', docFontSize); 
ylabel('Damping (\%)', 'Interpreter', 'latex' ,'FontSize', docFontSize); 
axis([vel(1) vel(end) -20 25]); 
grid on
grid minor;
ax = fig.CurrentAxes;
ax.TickLabelInterpreter = 'latex';
set(ax, 'Color', 'w');

f_Vg=fig;
%% -------------------------


Vg_PlotName = 'Vg_Plot_OLvsCL1vsCL2';
fullVg_PlotPath = fullfile(figures_dir, Vg_PlotName);
print(f_Vg, fullVg_PlotPath, '-dpdf', '-vector');
print(f_Vg, fullVg_PlotPath, '-dmeta', '-vector');


%% [f_Vev] = Vev_plot_mult_bw(vel_list, EV_list, legendlist);
% Plot the "root locus" movement of the flutter eigenvalues with increasing velocity

black = [0 0 0];  % all lines in black
linestyles = {'-', '--', ':', '-.'};  % cycling through these styles

if iscell(vel_)
    num_cases = length(vel_);
    assert(iscell(EV_) & length(EV_)==length(vel_),'Inconsistent input');
else
    num_cases=1;
    vel_={vel_};
    EV_={EV_};
end

fig = figure('Name','Vpzmap Plot Comparison');
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
for i = 1:num_cases
    vel = vel_{i};
    EV = EV_{i};
    ls = linestyles{mod(i-1,length(linestyles))+1};  % select line style

    % Plot eigenvalue paths with black color and varying line styles
    p{i} = plot(transpose(real(EV)), transpose(imag(EV)), ...
        'Color', black, 'LineStyle', ls, 'LineWidth', stdLineWidth);
    hold on

    % Plot the terminal point of the trajectory
    plot(transpose(real(EV(:,end))), transpose(imag(EV(:,end))), ...
        'Color', black, 'LineStyle', 'none', 'Marker', 'hexagram');
    
    firstlines(i) = p{i}(1);  % for legend
end

xline(0, '-k');     % Add vertical line at x=0 (dashed black)
yline(0, '-k');     % Add horizontal line at y=0 (dashed black)

lgd = legend(firstlines, legendlist, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', docFontSize);
lgd.AutoUpdate = 'off';
title('Eigenvalue Loci', 'Interpreter', 'latex', 'FontSize', docFontSize); 
xlabel('Real Part $\Re (\lambda)$', 'Interpreter', 'latex', 'FontSize', docFontSize); 
ylabel('Imaginary Part $\Im (\lambda)$', 'Interpreter', 'latex', 'FontSize', docFontSize); 
grid on; 
grid minor;
axis equal;
ylim([-40, 40]);

ax = fig.CurrentAxes;
ax.TickLabelInterpreter = 'latex';
set(ax, 'Color', 'w');
f_Vev = fig;
%% -------------------------
Vpzmap_PlotName = 'Vpzmap_Plot_OLvsCL1vsCL2';
fullVpzmap_PlotPath = fullfile(figures_dir, Vpzmap_PlotName);
print(f_Vev, fullVpzmap_PlotPath, '-dpdf', '-vector');
print(f_Vev, fullVpzmap_PlotPath, '-dmeta', '-vector');
