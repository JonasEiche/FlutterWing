function fig = plot_dlm_panel(Pa, Ps, j, opts)
%PLOT_DLM_PANEL  Numbered doublet-lattice grid with the doublet line and control point of one panel (TUTORIAL.m section (3)).
%   fig = plot_dlm_panel(Pa, Ps, j)
%   fig = plot_dlm_panel(Pa, Ps, 8, 'Name', 'goland_panels')
%   Pa, Ps   panel corners in aerodynamic and structural coordinates from build_PaPs
%   j        panel to mark: its 1/4-chord doublet line and its 3/4-chord control point, each with
%            a hairline leader to a label off the planform
%   The wing is drawn in the aerodynamic frame by wing_scene, with the panel numbers, both frame
%   triads and the freestream arrow. The leaders are placed for a panel in the outboard half of
%   the Goland grid.
%   Options (defaults)
%     'Name'   'DLM panels'   figure name (the tutorial export uses the file stem)
%   fig      fw_figure of the standard size
%   See also wing_scene, build_PaPs, build_Qjj.

arguments
    Pa cell
    Ps cell
    j (1,1) double {mustBeInteger, mustBePositive}
    opts.Name (1,:) char = 'DLM panels'
end

S = fw_style();
[fig, ax, scene] = wing_scene(struct('Pa',{Pa},'Ps',{Ps}), 'Frame', 'aero', 'PanelNumbers', true, ...
    'Freestream', true, 'Name', opts.Name);
Pj  = Pa{j};
dl1 = 0.75*Pj{1} + 0.25*Pj{2} + scene.lift(:);      % quarter chord on the inboard panel edge ...
dl2 = 0.75*Pj{4} + 0.25*Pj{3} + scene.lift(:);      % ... and on the outboard edge, lifted towards the camera
cpj = 0.5*((0.75*Pj{2}+0.25*Pj{1}) + (0.75*Pj{3}+0.25*Pj{4})) + scene.lift(:);
plot3(ax, [dl1(1) dl2(1)], [dl1(2) dl2(2)], [dl1(3) dl2(3)], '-', 'Color',S.blue, 'LineWidth',S.lineWidth)
plot3(ax, cpj(1), cpj(2), cpj(3), 'o', 'Color',S.blue, 'MarkerFaceColor',S.blue, 'MarkerSize',S.markerSize)
tj = findobj(ax, 'Type','text', 'String',num2str(j));   % the number makes room for the two marks
tj.Position = [min(Pj{1}(1),Pj{2}(1)) + 0.85*(Pj{2}(1)-Pj{1}(1)), Pj{1}(2) + 0.18*(Pj{4}(2)-Pj{1}(2)), 0] + scene.lift;
% leaders to the labels: the doublet line ahead of the leading edge, the control point aft of the
% trailing edge, where the scene is empty (both edges rise to the right on screen, so a label
% above the leading edge runs inboard and one below the trailing edge runs outboard)
lbl1 = [scene.xLE - 0.14*scene.chord, 0.62*scene.span, scene.z0];
lbl2 = [scene.xTE + 0.14*scene.chord, 0.72*scene.span, scene.z0];
plot3(ax, [dl1(1) lbl1(1)], [dl1(2) lbl1(2)], [dl1(3) lbl1(3)], '-', 'Color',S.muted, 'LineWidth',S.lineWidthHair)
text(ax, lbl1(1), lbl1(2), lbl1(3), '1/4-chord doublet line', 'HorizontalAlignment','right', ...
    'VerticalAlignment','bottom', 'FontSize',S.fontSizeSmall, 'Color',S.ink)
plot3(ax, [cpj(1) lbl2(1)], [cpj(2) lbl2(2)], [cpj(3) lbl2(3)], '-', 'Color',S.muted, 'LineWidth',S.lineWidthHair)
text(ax, lbl2(1), lbl2(2), lbl2(3), '3/4-chord control point $j$', 'HorizontalAlignment','left', ...
    'VerticalAlignment','top', 'FontSize',S.fontSizeSmall, 'Color',S.ink)
end
