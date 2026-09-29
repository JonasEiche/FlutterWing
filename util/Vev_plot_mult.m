function [fig] = Vev_plot_mult(vel_, EV_, legendlist)
%VEV_PLOT_MULT  Eigenvalue loci over velocity (root locus), one colour per case.
%   fig = Vev_plot_mult(V_inf, EV, legendlist)                  one case
%   fig = Vev_plot_mult({V1,V2,...}, {EV1,EV2,...}, legendlist)  several cases
%   V_inf       [1 x n_vel] velocities (m/s); only checked for consistency
%   EV          [n_modes x n_vel] tracked eigenvalues (rad/s) from
%               getEigenvalueModeshape or pkmethode, pre-filtered to the modes
%               of interest. One line per row from V_inf(1) to V_inf(end),
%               the end point marked with a hexagram.
%   legendlist  cell of labels, one per case (at most 6 colours defined)
%   fig         figure handle. Black print variant: Vev_plot_mult_bw.
%   "Vev" = velocity-eigenvalue plot. The imaginary axis is fixed to
%   [-40 40] rad/s (RectWing band, about 6.4 Hz); set ylim on the returned
%   figure for other models (Goland flutter mode near 69 rad/s).
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

fig = figure('Name','Vpzmap Plot Comparison');
set(fig,'defaultTextInterpreter','latex');
for i = 1:num_cases
% vel, EV
vel = vel_{i};
EV = EV_{i};
p{i} = plot(transpose(real(EV)), transpose(imag(EV)),'Color', colorlist{i},'LineStyle', "-",'LineWidth',1.2);
firstlines(i) = p{i}(1);
hold on
plot(transpose(real(EV(:,end))), transpose(imag(EV(:,end))),'Color', colorlist{i},'LineStyle', "none",'Marker', "hexagram")
% hold on
% text(real(EV(:,end))*1.01, imag(EV(:,end))*1.01, [num2str(vel(end)),'m/s'], 'Interpreter', 'latex');
end

xline(0, '-k');     % Add vertical line at x=0 (dashed black)
yline(0, '-k');     % Add horizontal line at y=0 (dashed black)

lgd=legend(firstlines,legendlist,'Location','northwest', 'Interpreter', 'latex');
lgd.AutoUpdate = 'off';
title(' Eigenvalue Loci', 'Interpreter', 'latex'); 
xlabel('Real Part $\Re (\lambda)$', 'Interpreter', 'latex'); ylabel('Imaginary Part $\Im (\lambda)$', 'Interpreter', 'latex'); grid on; axis equal; 

ylim([-40,40])

ax = fig.CurrentAxes; % Get current axes
ax.TickLabelInterpreter = 'latex'; % Set tick labels to LaTeX
set(ax, 'Color', 'w'); % Set axes background to white

end