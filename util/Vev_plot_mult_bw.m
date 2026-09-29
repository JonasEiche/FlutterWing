function [fig] = Vev_plot_mult_bw(vel_, EV_, legendlist)
%VEV_PLOT_MULT_BW  Eigenvalue loci over velocity (root locus), black with cycled line styles.
%   fig = Vev_plot_mult_bw(V_inf, EV, legendlist)                  one case
%   fig = Vev_plot_mult_bw({V1,V2,...}, {EV1,EV2,...}, legendlist)  several cases
%   V_inf       [1 x n_vel] velocities (m/s); only checked for consistency
%   EV          [n_modes x n_vel] tracked eigenvalues (rad/s) from
%               getEigenvalueModeshape or pkmethode, pre-filtered to the modes
%               of interest. One line per row from V_inf(1) to V_inf(end),
%               the end point marked with a hexagram.
%   legendlist  cell of labels, one per case
%   fig         figure handle. Colour variant: Vev_plot_mult.
%   "Vev" = velocity-eigenvalue plot. The imaginary axis is fixed to
%   [-40 40] rad/s (RectWing band, about 6.4 Hz); set ylim on the returned
%   figure for other models (Goland flutter mode near 69 rad/s).

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
for i = 1:num_cases
    vel = vel_{i};
    EV = EV_{i};
    ls = linestyles{mod(i-1,length(linestyles))+1};  % select line style

    % Plot eigenvalue paths with black color and varying line styles
    p{i} = plot(transpose(real(EV)), transpose(imag(EV)), ...
        'Color', black, 'LineStyle', ls, 'LineWidth', 1.2);
    hold on

    % Plot the terminal point of the trajectory
    plot(transpose(real(EV(:,end))), transpose(imag(EV(:,end))), ...
        'Color', black, 'LineStyle', 'none', 'Marker', 'hexagram');
    
    firstlines(i) = p{i}(1);  % for legend
end

xline(0, '-k');     % Add vertical line at x=0 (dashed black)
yline(0, '-k');     % Add horizontal line at y=0 (dashed black)

lgd = legend(firstlines, legendlist, 'Location', 'northwest', 'Interpreter', 'latex');
lgd.AutoUpdate = 'off';
title('Eigenvalue Loci', 'Interpreter', 'latex'); 
xlabel('Real Part $\Re (\lambda)$', 'Interpreter', 'latex'); 
ylabel('Imaginary Part $\Im (\lambda)$', 'Interpreter', 'latex'); 
grid on; 
axis equal;
ylim([-40, 40]);

ax = fig.CurrentAxes;
ax.TickLabelInterpreter = 'latex';
set(ax, 'Color', 'w');

end
