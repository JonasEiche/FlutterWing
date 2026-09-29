function fig = plot_mode_shapes(E, ele, PHIgf, OMEGA, modes, opts)
%PLOT_MODE_SHAPES  Bending deflection and twist of wind-off mode shapes along the span (TUTORIAL.m section (2)).
%   fig = plot_mode_shapes(E, ele, PHIgf, OMEGA)
%   fig = plot_mode_shapes(E, ele, PHIgf, OMEGA, 1:2, 'Name', 'goland_modes')
%   E, ele   beam grid from build_E_y
%   PHIgf    mode shapes from build_PHIgf, one column per mode
%   OMEGA    eigenfrequencies (rad/s) from build_PHIgf, for the legend
%   modes    columns of PHIgf to draw (default 1:2)
%   The mode shapes are interpolated between the nodes with the beam's own shape functions
%   (beamPHI) at 200 stations from root to tip: u_z in the top tile, psi_y below.
%   Options (defaults)
%     'Name'   'mode shapes'   figure name (the tutorial export uses the file stem)
%   fig      fw_figure of the standard size, 2x1 tiled layout
%   See also build_PHIgf, beamPHI, fw_figure.

arguments
    E cell
    ele cell
    PHIgf double
    OMEGA double
    modes (1,:) double = 1:2
    opts.Name (1,:) char = 'mode shapes'
end

S = fw_style();
y_plot = linspace(E{1}{1}(2), E{end}{2}(2), 200);   % root to tip
[PHIuz, PHIpsi] = beamPHI(y_plot, E, ele);
mode_lbl = arrayfun(@(m) sprintf('mode %d, %.2f Hz', m, OMEGA(m)/(2*pi)), modes, 'UniformOutput', false);

fig = fw_figure(S.size.standard(1), S.size.standard(2), 'Name', opts.Name);
tl  = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax  = nexttile(tl); hold(ax,'on')
plot(ax, y_plot, PHIuz*PHIgf(:,modes), 'LineWidth', S.lineWidth)
xlim(ax,y_plot([1 end])); ylabel(ax,'$u_z$ (mode shape units)')
title(ax,'Wind-off mode shapes of the clamped beam')
legend(ax, mode_lbl, 'Location','northwest'); set(ax,'XTickLabel',[])
ax  = nexttile(tl); hold(ax,'on')
plot(ax, y_plot, PHIpsi*PHIgf(:,modes), 'LineWidth', S.lineWidth)
xlim(ax,y_plot([1 end])); xlabel(ax,'$y$ (m)'); ylabel(ax,'$\psi_y$ (mode shape units)')
end
