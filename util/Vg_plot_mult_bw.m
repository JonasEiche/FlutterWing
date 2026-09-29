function [fig] = Vg_plot_mult_bw(vel_, EV_, legendlist)
%VG_PLOT_MULT_BW  V-g / V-f flutter plot, black with cycled line styles (print version).
%   fig = Vg_plot_mult_bw(V_inf, EV, legendlist)                  one case
%   fig = Vg_plot_mult_bw({V1,V2,...}, {EV1,EV2,...}, legendlist)  several cases
%   V_inf       [1 x n_vel] velocities (m/s), row vector
%   EV          [n_modes x n_vel] tracked eigenvalues (rad/s) from
%               getEigenvalueModeshape or pkmethode, pre-filtered to the modes
%               of interest (e.g. EV(1:4,:))
%   legendlist  cell of labels, one per case
%   fig         figure handle. Top: frequency |lambda|/(2*pi) in Hz; bottom:
%               damping -100*Re(lambda)/|lambda| in %. Flutter points are red
%               circles, divergence points blue circles. Colour variant: Vg_plot_mult.
%   Console output: for every mode the first grid point with Re(lambda) > tol
%   (tol = 0.1) is printed as flutter if Im(lambda) > divtol (divtol = 1 rad/s),
%   else as divergence. It is a grid point, not an interpolated crossing:
%   refine V_inf near the crossing for an accurate speed. A pole that returns
%   to the left half-plane later is not reported. The same rule is used by
%   find_flutter_velocity in R2_Sensor_Fault_Tolerance.m (IFASD2026).
%   Layout as in the paper V-g figures (research-paper-code/.../fig11_vg_diagram):
%   2x1 tiled layout with compact spacing, grid and minor grid, x tick labels only
%   on the damping tile (shared velocity axis).
%   Axis limits are fixed: x = [V_inf(1) V_inf(end)] of the last case, damping
%   axis [-20 25] %; change them with axis(...) on the returned figure.

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
    div_ind = (cumsum((real(EV) > tol),2) == 1 & imag(EV)<divtol & imag(EV)>=0);

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
if isprop(fig, 'Theme'), fig.Theme = 'light'; end   % R2025a+: never inherit a dark desktop theme
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(tl, 1)
for i = 1:num_cases
    ls = linestyles{mod(i-1,length(linestyles))+1};  % cycle through line styles
    plot(vel_{i},F{i},'Color',black,'LineStyle',ls,'LineWidth',1.2)
    hold on
    plot(DIV_V{i},zeros(size(DIV_V{i})),'Ob', 'MarkerSize', 4, 'LineWidth', 2)  % divergence
    hold on
    plot(FLUT_V{i},FLUT_F{i},'Or', 'MarkerSize', 4, 'LineWidth', 2)             % flutter
    hold on
end
title('Frequency and Damping vs. Velocity', 'Interpreter', 'latex'); 
ylabel('Frequency (Hz)', 'Interpreter', 'latex');
xlim([vel(1) vel(end)]);
grid on
grid minor

ax = fig.CurrentAxes;
ax.TickLabelInterpreter = 'latex';
set(ax, 'Color', 'w');
set(ax, 'XTickLabel', []);   % shared velocity axis: tick labels only on the damping tile

%% Plot Damping over Speed
nexttile(tl, 2)
x0 = [vel(1) vel(end)]; y0 = [0 0];

firstlines = zeros(1,num_cases);
for i = 1:num_cases
    ls = linestyles{mod(i-1,length(linestyles))+1};
    p{i} = plot(vel_{i},D{i},'Color',black,'LineStyle',ls,'LineWidth',1.2);
    hold on
    plot(DIV_V{i},-100*ones(size(DIV_V{i})),'Ob', 'MarkerSize', 4, 'LineWidth', 2)
    hold on
    plot(FLUT_V{i},zeros(size(FLUT_V{i})),'Or', 'MarkerSize', 4, 'LineWidth', 2)
    hold on
    firstlines(i) = p{i}(1);
end
plot(x0,y0,'--k')
lgd = legend(firstlines,legendlist,'Location','northwest', 'Interpreter', 'latex');
lgd.AutoUpdate = 'off';
xlabel('Velocity (m/s)', 'Interpreter', 'latex'); 
ylabel('Damping (\%)', 'Interpreter', 'latex'); 
axis([vel(1) vel(end) -20 25]); 
grid on
grid minor

ax = fig.CurrentAxes;
ax.TickLabelInterpreter = 'latex';
set(ax, 'Color', 'w');

end
