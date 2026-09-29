% Vg_Plot_and_Vpzmap_OL_vs_H2_vs_SB
%
% V-g diagram over V_inf = 20-160 m/s for open loop, H2 blending and
% structural blending (the eigenvalue-locus map at the end is drawn but
% not printed).
% Produces: Figures/Fig20_VgPlot_OLvsH2vsSB.pdf
% Requires: Data/controller_imu18_ail18_structural_blending.mat (from R0)
% Recomputes the eigenvalue sweep from the controller file on every run. The
% Data/Vg_Vpz_OLvsH2vsSB.mat written below is a by-product for inspection
% only: it is git-ignored and read by no script.

clearvars
textwidth = 15.98; % cm (IFASD 2026 template: 455.24pt)
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
ContPath = fullfile(data_dir, 'controller_imu18_ail18_structural_blending.mat');
load(ContPath, 'cont_H2', 'cont_SB')

%  -5-6-7-8-
% |         |
%  -1-2-3-4-
num_modes = 5;     % must equal num_modes hard-coded in build_G_RectWing / build_P_RectWing
imuIDX = 1:8;
ailIDX = 1:8;
modesIDX = [1,2];
nimu = length(imuIDX);
nail = length(ailIDX);
nmod = length(modesIDX);

V_inf=linspace(20,160,32);
nv = length(V_inf);
G = build_G_RectWing(V_inf,imuIDX,ailIDX);

CL_H2 = feedback(G,-cont_H2);
CL_SB = feedback(G,-cont_SB);

%% COMPUTE EIGENVALUES / MODESHAPES
[EV_OL,MS_OL] = getEigenvalueModeshape(G,num_modes);
[EV_H2,MS_H2] = getEigenvalueModeshape(CL_H2,num_modes);
[EV_SB,MS_SB] = getEigenvalueModeshape(CL_SB,num_modes);

num_modes_disp=2;
EV_H2_disp = EV_H2(1:2*num_modes_disp,:);
EV_SB_disp = EV_SB(1:2*num_modes_disp,:);
EV_OL_disp = EV_OL(1:2*num_modes_disp,:);

legendlist={"Open Loop", "$H_2$ Optimal Blending", "Structural Blending"};
EV_list={EV_OL_disp,EV_H2_disp,EV_SB_disp};
vel_list={V_inf,V_inf,V_inf};
vel_=vel_list;
EV_=EV_list;

%% COMPUTE FLUTTER / DIVERGENCE / FREQUENCY / DAMPING PER CASE
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

% Extract per-case results into suffixed variables for saving
F_OL = F{1};        F_H2 = F{2};        F_SB = F{3};
D_OL = D{1};        D_H2 = D{2};        D_SB = D{3};
FLUT_F_OL = FLUT_F{1};  FLUT_F_H2 = FLUT_F{2};  FLUT_F_SB = FLUT_F{3};
FLUT_V_OL = FLUT_V{1};  FLUT_V_H2 = FLUT_V{2};  FLUT_V_SB = FLUT_V{3};
DIV_V_OL = DIV_V{1};    DIV_V_H2 = DIV_V{2};    DIV_V_SB = DIV_V{3};


%% SAVE RESULTS TO MAT FILE
vgDataPath = fullfile(data_dir, 'Vg_Vpz_OLvsH2vsSB.mat');
save(vgDataPath, ...
    'V_inf', ...
    'EV_OL', 'MS_OL', 'EV_OL_disp', ...
    'EV_H2', 'MS_H2', 'EV_H2_disp', ...
    'EV_SB', 'MS_SB', 'EV_SB_disp', ...
    'F_OL', 'D_OL', 'FLUT_F_OL', 'FLUT_V_OL', 'DIV_V_OL', ...
    'F_H2', 'D_H2', 'FLUT_F_H2', 'FLUT_V_H2', 'DIV_V_H2', ...
    'F_SB', 'D_SB', 'FLUT_F_SB', 'FLUT_V_SB', 'DIV_V_SB');


%% Vg Plot  //  Vpzmap
% Define color and linestyles
black = [0 0 0];
linestyles = {'-', '--', ':', '-.'};  % add more if needed

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
Vg_PlotName = 'Fig20_VgPlot_OLvsH2vsSB';
fullVg_PlotPath = fullfile(figures_dir, Vg_PlotName);
print(f_Vg, fullVg_PlotPath, '-dpdf', '-vector');
% print(f_Vg, fullVg_PlotPath, '-dmeta', '-vector');


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
Vpzmap_PlotName = 'FigXX_Vpzmap_OLvsH2vsSB';
fullVpzmap_PlotPath = fullfile(figures_dir, Vpzmap_PlotName);
% print(f_Vev, fullVpzmap_PlotPath, '-dpdf', '-vector');
% print(f_Vev, fullVpzmap_PlotPath, '-dmeta', '-vector');
