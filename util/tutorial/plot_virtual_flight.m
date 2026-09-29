function fig = plot_virtual_flight(t, z_tip, u_flap, t_on, V_on, z_thresh, opts)
%PLOT_VIRTUAL_FLIGHT  Tip deflection and flap command around the controller switch-on (TUTORIAL.m section (11)).
%   fig = plot_virtual_flight(t, z_tip, u_flap, t_on, V_on, z_thresh)
%   fig = plot_virtual_flight(t_sim, zrec, urec, t_sim(i_on), V_t(i_on), z_thresh, 'Name', 'virtual_flight')
%   t        [Nt x 1] time (s)
%   z_tip    [Nt x 2] vertical tip deflection (m) at the leading and the trailing edge
%   u_flap   [Nt x 1] flap command (rad), zero before the switch-on
%   t_on     switch-on time (s), drawn as an ink line labelled AFS ON
%   V_on     velocity at the switch-on (m/s), for the title
%   z_thresh switch-on threshold on the leading-edge deflection (m), dotted at +/- z_thresh
%   Tip deflection on top (leading edge blue, trailing edge muted), flap command in degrees below.
%   Options (defaults)
%     'Window' [-0.6 0.9]          time window around t_on (s)
%     'Name'   'virtual flight'    figure name (the tutorial export uses the file stem)
%   fig      fw_figure of the standard size, 2x1 tiled layout
%   See also build_LPV_P_Goland, lsim, lft.

arguments
    t (:,1) double
    z_tip (:,2) double
    u_flap (:,1) double
    t_on (1,1) double
    V_on (1,1) double
    z_thresh (1,1) double
    opts.Window (1,2) double = [-0.6 0.9]
    opts.Name (1,:) char = 'virtual flight'
end

S = fw_style();
t_zoom = t_on + opts.Window;                         % window around the switch-on (s)

fig = fw_figure(S.size.standard(1), S.size.standard(2), 'Name', opts.Name);
tl  = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax  = nexttile(tl); hold(ax,'on')
plot(ax, t, z_tip(:,1), '-', 'Color',S.blue,  'LineWidth',S.lineWidthThin)
plot(ax, t, z_tip(:,2), '-', 'Color',S.muted, 'LineWidth',S.lineWidthThin)
yline(ax,  z_thresh, ':', 'Color',S.muted, 'LineWidth',S.lineWidthHair);
yline(ax, -z_thresh, ':', 'Color',S.muted, 'LineWidth',S.lineWidthHair);
xline(ax, t_on, '-', 'Color',S.ink, 'LineWidth',S.lineWidth);
xlim(ax,t_zoom); ylim(ax,[-0.35 0.35]); ylabel(ax,'tip deflection (m)')
title(ax, sprintf('Virtual flight: the switch-on at %.0f m/s', V_on))
legend(ax, {'leading edge','trailing edge'}, 'Location','southeast')
ax  = nexttile(tl); hold(ax,'on')
plot(ax, t, u_flap*180/pi, '-', 'Color',S.blue, 'LineWidth',S.lineWidthThin)
xline(ax, t_on, '-', 'AFS ON', 'Color',S.ink, 'LineWidth',S.lineWidth, ...     % the flap is idle before it: room for the label
    'LabelOrientation','horizontal', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','top');
xlim(ax,t_zoom); ylim(ax,[-30 30]); xlabel(ax,'time (s)'); ylabel(ax,'flap command (deg)')
end
