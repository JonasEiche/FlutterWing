function [fig] = Vg_plot_mult(vel_, EV_, legendlist)
%VG_PLOT_MULT  V-g / V-f flutter plot, one colour per case (screen version).
%   fig = Vg_plot_mult(V_inf, EV, legendlist)                  one case
%   fig = Vg_plot_mult({V1,V2,...}, {EV1,EV2,...}, legendlist)  several cases
%   V_inf       [1 x n_vel] velocities (m/s), row vector
%   EV          [n_modes x n_vel] tracked eigenvalues (rad/s) from
%               getEigenvalueModeshape or pkmethode, pre-filtered to the modes
%               of interest (e.g. EV(1:4,:))
%   legendlist  cell of labels, one per case (at most 6 colours defined)
%   fig         figure handle. Top: frequency |lambda|/(2*pi) in Hz; bottom:
%               damping -100*Re(lambda)/|lambda| in %. Flutter points are red
%               circles, divergence points blue circles. Black print variant used by the
%               tutorial and the paper scripts: Vg_plot_mult_bw.
%   Console output: for every mode the first grid point with Re(lambda) > tol
%   (tol = 0.1) is printed as flutter if Im(lambda) > divtol (divtol = 1 rad/s),
%   else as divergence. It is a grid point, not an interpolated crossing:
%   refine V_inf near the crossing for an accurate speed. A pole that returns
%   to the left half-plane later is not reported. The same rule is used by
%   find_flutter_velocity in R2_Sensor_Fault_Tolerance.m (IFASD2026).
%   Axis limits are fixed: x = [V_inf(1) V_inf(end)] of the last case, damping
%   axis [-20 25] %; change them with axis(...) on the returned figure.

black = [0 0 0];
darkblue = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
gelb = [0.9290 0.6940 0.1250];
lila = [0.4940 0.1840 0.5560];
gruen = [0.4660 0.6740 0.1880];
colorlist={black, darkblue, orange, lila, gruen, gelb};

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
% vel, EV
vel = vel_{i};
EV = EV_{i};
% EV   : [num_X,num_vel]           eigenvalue for every state x velocities
% Tolerance what is still considered a stable pole Re up to tol in the RHP
tol = 0.1;
% Divergence Tolerance is the Im oscillatory instability still considered to be not oscillatory ie. divergence
divtol = 1.0;
% logical array [n_eigs x n_vel]: for every eigenvalue (row) the first velocity with real(ev) > tol is marked with a 1
crit_ind = cumsum((real(EV) > tol),2) == 1;
% flut_ind / div_ind: the same first crossing, split into oscillatory (flutter) and non-oscillatory (divergence)
% a pole that later returns to the stable left half-plane is not detected because the cumsum stays 1
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
if sum(sum(crit_ind)) == 0                             % number of total instable eigenvalues
    disp('No flutter or divergence')
end
f = abs(EV)./(2*pi);
% zeta = -real(EV)./abs(EV);
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
subplot(211)
for i = 1:num_cases
plot(vel_{i},F{i},'Color', colorlist{i},'LineWidth',1.2)
hold on

plot(DIV_V{i},zeros(size(DIV_V{i})),'Ob')% divergence
hold on
plot(FLUT_V{i},FLUT_F{i},'Or', 'MarkerSize', 4, 'LineWidth', 2)   % flutter
hold on
end
title('Frequency and Damping vs. Velocity', 'Interpreter', 'latex'); 
% xlabel('Velocity (m/s)'); 
ylabel('Frequency (Hz)', 'Interpreter', 'latex');
xlim([vel(1) vel(end)]);
grid on

ax = fig.CurrentAxes; % Get current axes
ax.TickLabelInterpreter = 'latex'; % Set tick labels to LaTeX
set(ax, 'Color', 'w'); % Set axes background to white

%% Plot Damping over Speed
subplot(212)
x0 = [vel(1) vel(end)]; y0 = [0 0];
% plot(vel1,D{1},':k',vel2,D{2},'k',x0,y0,'--k')

firstlines = zeros(1,num_cases);
for i = 1:num_cases
p{i} = plot(vel_{i},D{i},'Color', colorlist{i},'LineWidth',1.2);
hold on
plot(DIV_V{i},-100*ones(size(DIV_V{i})),'Ob')
hold on
plot(FLUT_V{i},zeros(size(FLUT_V{i})),'Or')
hold on
firstlines(i) = p{i}(1);
end
plot(x0,y0,'--k')
lgd=legend(firstlines,legendlist,'Location','northwest', 'Interpreter', 'latex');
lgd.AutoUpdate = 'off';
xlabel('Velocity (m/s)', 'Interpreter', 'latex'); ylabel('Damping (\%)', 'Interpreter', 'latex'); 
axis([vel(1) vel(end) -20 25]); 
grid on

ax = fig.CurrentAxes; % Get current axes
ax.TickLabelInterpreter = 'latex'; % Set tick labels to LaTeX
set(ax, 'Color', 'w'); % Set axes background to white








