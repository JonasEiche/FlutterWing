%% fig08_Vg_Diagram.m
% Paper Fig. 8. Generates: fig08_vg_diagram (.pdf, .emf)
% EMF (-dmeta) is written on Windows only; other platforms get the PDF.
%
% V-g diagram: frequency (top) and damping (bottom) of the first two
% aeroelastic modes versus freestream velocity for open loop, orthogonal
% modal blending, and H2 optimal blending. Flutter and divergence onset
% points are marked with colored circles. Shows how each controller shifts
% the flutter boundary.
%
% The eigenvalue loci are read from the R06 archive (stage 7 of
% R06_Locked_Comparison.m: build_G_RectWing on linspace(20,160,32),
% getEigenvalueModeshape, locked MB2orth and H2Pusch controllers).
%
% Requires: Data/Locked_Comparison_Results.mat

clearvars

textwidth = 0.0351*372.0;    % cm/pt * latex template pt;  % cm
figwidthVg = 0.9*textwidth; % cm
figheightVg = 0.7*figwidthVg;

docFontSize = 9; % pt
stdLineWidth = 1.2;
archive = 'Locked_Comparison_Results.mat';   % or Locked_Comparison_Results_test.mat

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end

%% Load the eigenvalue loci of the R06 archive
data_dir = fullfile(script_dir, 'Data');
A = load(fullfile(data_dir, archive), 'vg');

V_inf    = A.vg.Vsweep;
EV_OL    = A.vg.EV_OL;
EV_Cont1 = A.vg.MB2orth.EV;   % orthogonal modal blending
EV_Cont2 = A.vg.H2Pusch.EV;   % H2-optimal blending

num_modes_disp = 2;
EV_OL_disp    = EV_OL(1:2*num_modes_disp,:);
EV_Cont1_disp = EV_Cont1(1:2*num_modes_disp,:);
EV_Cont2_disp = EV_Cont2(1:2*num_modes_disp,:);

legendlist = {"Open Loop", "Orth. Modal", "$H_2$ Blending"};
EV_list = {EV_OL_disp, EV_Cont1_disp, EV_Cont2_disp};
vel_list = {V_inf, V_inf, V_inf};
vel_ = vel_list;
EV_ = EV_list;

%% Compute frequency and damping
black = [0 0 0];
linestyles = {'-', '--', ':', '-.'};

num_cases = length(vel_);

FLUT_F = cell(1,num_cases);
FLUT_V = cell(1,num_cases);
DIV_V  = cell(1,num_cases);
F = cell(1,num_cases);
D = cell(1,num_cases);

for i = 1:num_cases
    vel = vel_{i};
    EV = EV_{i};

    tol = 0.1;
    divtol = 1.0;

    crit_ind = cumsum((real(EV) > tol),2) == 1;
    flut_ind = (cumsum((real(EV) > tol),2) == 1 & imag(EV)>divtol);
    div_ind  = (cumsum((real(EV) > tol),2) == 1 & imag(EV)<divtol & imag(EV)>0);

    flut_f = abs(EV(flut_ind))./(2*pi);
    VEL = repmat(vel,size(EV,1),1);
    flut_v = VEL(flut_ind);
    div_v  = VEL(div_ind);

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
    DIV_V{i}  = div_v;
    F{i} = f;
    D{i} = d;
end

%% Vg Plot: Frequency over Speed
fig = figure('Name','Vg Plot OL vs Orth vs H2');
set(fig,'defaultTextInterpreter','latex');
set(fig, 'Color', 'w');
set(fig, 'Units','centimeters', ...
         'Position', [7,7,figwidthVg,figheightVg], ...
         'defaultAxesFontSize', docFontSize, ...
         'defaultTextFontSize', docFontSize, ...
         'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex');
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperSize', [figwidthVg figheightVg]);
set(fig, 'PaperPosition', [0 0 figwidthVg figheightVg]);
axFreq = subplot(211);
firstlines = gobjects(1,num_cases);
for i = 1:num_cases
    ls = linestyles{mod(i-1,length(linestyles))+1};
    pFreq = plot(vel_{i},F{i},'Color',black,'LineStyle',ls,'LineWidth',stdLineWidth);
    hold on
    firstlines(i) = pFreq(1);
    plot(DIV_V{i},zeros(size(DIV_V{i})),'Ob', 'MarkerSize', 4, 'LineWidth', 2)
    hold on
    plot(FLUT_V{i},FLUT_F{i},'Or', 'MarkerSize', 4, 'LineWidth', 2)
    hold on
end
lgd = legend(axFreq,firstlines,legendlist,'Location','best', 'Interpreter', 'latex', 'FontSize', docFontSize-1.7);
lgd.AutoUpdate = 'off';
title('Frequency and Damping vs. Velocity', 'Interpreter', 'latex', 'FontSize', docFontSize);
ylabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontSize', docFontSize);
xlim([vel_{1}(1) vel_{1}(end)]);
grid on
grid minor
ax = fig.CurrentAxes;
ax.TickLabelInterpreter = 'latex';
set(ax, 'Color', 'w');

%% Damping over Speed
subplot(212)
x0 = [vel_{1}(1) vel_{1}(end)]; y0 = [0 0];

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
xlabel('Velocity (m/s)', 'Interpreter', 'latex', 'FontSize', docFontSize);
ylabel('Damping (\%)', 'Interpreter', 'latex' ,'FontSize', docFontSize);
axis([vel_{1}(1) vel_{1}(end) -20 25]);
grid on
grid minor;
ax = fig.CurrentAxes;
ax.TickLabelInterpreter = 'latex';
set(ax, 'Color', 'w');

f_Vg = fig;

%% Save Vg Plot
Vg_PlotName = 'fig08_vg_diagram';
fullVg_PlotPath = fullfile(figures_dir, Vg_PlotName);
print(f_Vg, fullVg_PlotPath, '-dpdf', '-vector');
if ispc, print(f_Vg, fullVg_PlotPath, '-dmeta', '-vector'); end  % EMF export is Windows-only
