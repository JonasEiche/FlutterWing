function [fig, ax, scene] = wing_scene(Structure, opts)
%WING_SCENE  3-D view of a wing model: panel mesh, control surfaces, IMU arrows, axis overlays, frame triads.
%   fig = wing_scene(Structure)                       panel mesh in structural coordinates with the triad
%   fig = wing_scene(Structure, 'Surfaces', 'all', 'IMUs', 'all')
%                                                     control surfaces in grey, blue IMU arrows with numbers
%   fig = wing_scene(Structure, 'Axes', {'flexural axis', 0; 'mass axis', x_m})
%                                                     spanwise lines at chordwise position x (first solid,
%                                                     the others dashed)
%   fig = wing_scene(Structure, 'Frame', 'aero', 'PanelNumbers', true)
%                                                     panel indices, drawn from Structure.Pa
%   fig = wing_scene(Structure, 'Triad', 'both', 'Freestream', true)
%                                                     both frame triads at their origins, the flow arrow
%   [fig, ax, scene] = wing_scene(...)                the axes and the scene geometry, for more marks
%   Structure   from define_*_Structure_Aero: Ps (or Pa with 'Frame' 'aero') and cspanels are read;
%               'Triad' 'both' needs Pa and Ps
%   Options (defaults)
%     'Frame'        'structural'  'structural' draws Structure.Ps (x upstream, z down), 'aero' draws
%                                  Structure.Pa (x downstream, z up). The coordinates are drawn as they
%                                  are, so both frames keep their right-handed sense, and the camera
%                                  gives both the same picture: it stands aft of the wing, above it and
%                                  towards the root, physical up on top, so the trailing edge is nearest
%                                  and the tip runs to the right as in the plan view of plot_wing_layout.
%                                  Only the triad differs (z points down in the structural frame).
%     'Surfaces'     []            indices into Structure.cspanels drawn grey with ink hinge lines;
%                                  'all' draws every control surface
%     'Selected'     []            surfaces drawn blue instead of grey (the ones in the loop)
%     'IMUs'         []            IMU numbers: a blue arrow (physical up) at the hinge-line centre of
%                                  surface k with its number; 'all' for every surface
%     'IMUNumbers'   true          the number above every IMU arrow; false for a wing with one sensor
%     'PanelNumbers' false         panel index at every panel centre (8 pt)
%     'Axes'         {}            [n x 2] cell {label, x} or [n x 3] cell {label, x, colour}: spanwise
%                                  lines at chordwise position x in the drawn frame; the first is
%                                  solid, the others dashed, ink unless a colour is given; the labels
%                                  sit ahead of the leading edge around mid span, off the planform,
%                                  each with a hairline leader to its line
%     'Triad'        true          x, y, z arrows of the drawn frame at its origin: the root leading
%                                  edge (aero) or the root end of the flexural axis (structural),
%                                  with a white casing so they read over the outline and the mesh;
%                                  'both' adds the other frame (x_a + x_s = x_f c, y_a = y_s,
%                                  z_a = -z_s, build_PaPs), labelled x_a ... z_a and x_s ... z_s.
%                                  Its origin lies on the root edge a fraction of a chord from the
%                                  drawn one, where two full-size triads would overlap, so it is
%                                  drawn floated inboard of the root with a dotted leader to its
%                                  true origin, as on a drafting sheet; false draws none
%     'Freestream'   false         arrow ahead of the leading edge near the root, pointing downstream,
%                                  labelled V_inf at its tail
%     'View'         [37 26]       azimuth and elevation of the camera in the z-up sense of MATLAB's
%                                  view (the animate_wing default)
%     'Size'         fw_style().size.standard   cm, for a new figure
%     'Parent'       []            figure or axes to draw into
%     'Name'         'wing scene'  figure name
%   scene   struct: n_cam (unit vector towards the camera), up (physical up in the drawn frame),
%           lift (a small step towards the camera: add it to a point drawn on the wing plane so that
%           it wins the depth test against the panel faces), span, chord, xLE, xTE, z0 (the drawn
%           coordinate of the wing plane)
%   Corner order of every panel: 1 LE-inboard, 2 TE-inboard, 3 TE-outboard, 4 LE-outboard (build_PaPs);
%   the hinge of a trailing-edge surface is its upstream edge, that of a leading-edge surface its
%   downstream edge, as in plot_wing_layout and animate_wing. Rectangular unswept planforms.
%   Example:
%     [Structure, ~] = define_RectWing_Structure_Aero(5, 6);
%     wing_scene(Structure, 'Surfaces', 'all', 'IMUs', 'all');
%   See also plot_wing_layout, animate_wing, fw_figure, fw_style.

arguments
    Structure (1,1) struct
    opts.Frame (1,:) char {mustBeMember(opts.Frame, {'structural','aero'})} = 'structural'
    opts.Surfaces = []
    opts.Selected (1,:) double = []
    opts.IMUs = []
    opts.IMUNumbers (1,1) logical = true
    opts.PanelNumbers (1,1) logical = false
    opts.Axes cell = {}
    opts.Triad = true
    opts.Freestream (1,1) logical = false
    opts.View (1,2) double = [37 26]
    opts.Size (1,2) double {mustBePositive} = fw_style().size.standard
    opts.Parent = []
    opts.Name (1,:) char = 'wing scene'
end

S = fw_style();
if ischar(opts.Triad) || isstring(opts.Triad)
    assert(strcmpi(char(opts.Triad), 'both'), 'wing_scene:triad', 'Triad must be true, false or ''both''.');
    triads = 'both';
    assert(isfield(Structure, 'Pa') && isfield(Structure, 'Ps'), 'wing_scene:fields', ...
        'Triad ''both'' needs the fields Pa and Ps.');
    x_other = Structure.Pa{1}{1}(1) + Structure.Ps{1}{1}(1);   % x_a + x_s = x_f c: the other origin
elseif opts.Triad
    triads = 'drawn';
else
    triads = 'none';
end
if strcmp(opts.Frame, 'aero')
    assert(isfield(Structure, 'Pa'), 'wing_scene:fields', 'Structure needs the field Pa for the aero frame.');
    P    = Structure.Pa;
    up   = [0 0 1];                        % physical up in the drawn frame
    s_up = -1;                             % x increases downstream: the upstream coordinate is -x
else
    assert(isfield(Structure, 'Ps'), 'wing_scene:fields', 'Structure needs the field Ps.');
    P    = Structure.Ps;
    up   = [0 0 -1];
    s_up = 1;
end
cs = {};
if isfield(Structure, 'cspanels'), cs = Structure.cspanels; end
n_cs = numel(cs);
surfaces = local_indices(opts.Surfaces, n_cs, 'Surfaces');
imus     = local_indices(opts.IMUs, n_cs, 'IMUs');
if ~isempty(opts.Axes) && ~any(size(opts.Axes, 2) == [2 3])
    error('wing_scene:axes', 'Axes must be an [n x 2] cell {label, x} or [n x 3] cell {label, x, colour}.');
end

% ---- geometry, drawn in the coordinates of the frame --------------------
n_pan = numel(P);
V = zeros(4*n_pan, 3);
for i = 1:n_pan
    for n = 1:4
        V(4*(i-1)+n, :) = P{i}{n}(:).';
    end
end
Fa   = reshape(1:4*n_pan, 4, n_pan)';
x_up = s_up*V(:,1);                        % chordwise coordinate that increases upstream
xLE  = s_up*max(x_up);  xTE = s_up*min(x_up);
c    = max(x_up) - min(x_up);
span = max(V(:,2));
z0   = mean(V(:,3));
tol  = 1e-6*max(c, 1);

% ---- figure / axes ------------------------------------------------------
if isempty(opts.Parent)
    fig = fw_figure(opts.Size(1), opts.Size(2), 'Name', opts.Name);
    ax  = axes(fig, 'Position', [0.02 0.02 0.96 0.96]);
elseif isgraphics(opts.Parent, 'figure')
    fig = opts.Parent;
    ax  = axes(fig);
elseif isgraphics(opts.Parent) && isprop(opts.Parent, 'XLim')
    ax  = opts.Parent;
    fig = ancestor(ax, 'figure');
else
    error('wing_scene:parent', 'Parent must be [], a figure handle or an axes handle.');
end
hold(ax, 'on');

% ---- panel mesh and planform outline ------------------------------------
patch(ax, 'Faces',Fa, 'Vertices',V, 'FaceColor',S.white, ...
    'EdgeColor',S.ink, 'EdgeAlpha',0.35, 'LineWidth',S.lineWidthMesh);
plot3(ax, [xLE xTE xTE xLE xLE], [0 0 span span 0], z0*ones(1,5), '-', ...
    'Color',S.ink, 'LineWidth',S.lineWidth);

% ---- control surfaces and hinge lines ----------------------------------
hinge = zeros(n_cs, 2, 3);                 % hinge end points of every surface
for k = 1:n_cs
    blk  = cs{k};
    xblk = x_up(reshape(Fa(blk,:)', [], 1));
    isTE = abs(min(xblk) - min(x_up)) < tol;
    if isTE                                % flap: hinge on the upstream (LE-side) edge
        h1 = P{blk(1)}{1}(:).';  h2 = P{blk(end)}{4}(:).';
    else                                   % slat: hinge on the downstream (TE-side) edge
        h1 = P{blk(1)}{2}(:).';  h2 = P{blk(end)}{3}(:).';
    end
    hinge(k,1,:) = h1;  hinge(k,2,:) = h2;
    if ismember(k, surfaces)
        if ismember(k, opts.Selected), col = S.blue; else, col = S.grey; end
        patch(ax, 'Faces',Fa(blk,:), 'Vertices',V, 'FaceColor',col, ...
            'EdgeColor',S.ink, 'EdgeAlpha',0.35, 'LineWidth',S.lineWidthMesh);
        plot3(ax, [h1(1) h2(1)], [h1(2) h2(2)], [h1(3) h2(3)], '-', ...
            'Color',S.ink, 'LineWidth',S.lineWidth);
    end
end

% ---- camera: aft of the wing, above, towards the root; the same picture in both frames ----
axis(ax, 'equal');  axis(ax, 'off');
xlim(ax, sort([xTE xLE]) + [-0.2 0.2]*c);
ylim(ax, [-0.36*span, 1.1*span]);
zlim(ax, z0 + [-0.16 0.16]*span);
az = opts.View(1);  el = opts.View(2);
n_cam = [sind(az)*cosd(el), -cosd(az)*cosd(el), sind(el)];  % view(az, el) with z up
if up(3) < 0, n_cam = [-n_cam(1), n_cam(2), -n_cam(3)]; end % the same camera with z down
ax.CameraTarget        = [mean(xlim(ax)), mean(ylim(ax)), mean(zlim(ax))];
ax.CameraPosition      = ax.CameraTarget + 10*span*n_cam;
ax.CameraUpVector      = up;
ax.Projection          = 'orthographic';
ax.CameraViewAngleMode = 'auto';
lift = 0.02*span*n_cam;                    % towards the camera: invisible in orthographic projection,
                                           % but text and triad win the depth test against the faces

% ---- axis overlays ------------------------------------------------------
% Labels ahead of the leading edge around mid span, off the planform, stacked outboard to
% inboard and running inboard (the leading edge falls away to the left on screen, so a label
% that ends at its anchor never crosses it); a hairline leader runs from each line to the lower
% right corner of its label.
for r = 1:size(opts.Axes, 1)
    x_r = opts.Axes{r,2};
    if r == 1, ls = '-'; else, ls = '--'; end
    col = S.ink;
    if size(opts.Axes, 2) == 3 && ~isempty(opts.Axes{r,3}), col = opts.Axes{r,3}; end
    plot3(ax, [x_r x_r], [0 span], [z0 z0], ls, 'Color',col, 'LineWidth',S.lineWidth);
    y_r = (0.70 - 0.16*(r-1))*span;
    p   = [xLE + s_up*0.40*c, y_r, z0] + 10*lift;
    q   = [x_r, y_r, z0] + 10*lift;
    plot3(ax, [q(1) p(1)], [q(2) p(2)], [q(3) p(3)], '-', 'Color',S.muted, 'LineWidth',S.lineWidthHair);
    text(ax, p(1), p(2), p(3), opts.Axes{r,1}, 'Interpreter',S.interpreter, ...
        'FontSize',S.fontSize, 'Color',S.ink, 'HorizontalAlignment','right', 'VerticalAlignment','bottom');
end
% ---- panel numbers ------------------------------------------------------
if opts.PanelNumbers
    for i = 1:n_pan
        cpt = mean(V(Fa(i,:), :), 1) + lift;
        text(ax, cpt(1), cpt(2), cpt(3), sprintf('%d', i), 'Interpreter',S.interpreter, ...
            'FontSize',S.fontSizeSmall, 'Color',S.ink, ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle');
    end
end

% ---- IMU arrows and numbers --------------------------------------------
L = 0.10*span;
for k = imus
    base = 0.5*(squeeze(hinge(k,1,:)) + squeeze(hinge(k,2,:))).' + lift;
    local_arrow(ax, base, up, L, S.blue, 0.035, n_cam, [], true, [], 'both');
    if opts.IMUNumbers
        tip = base + up*(L + 0.035*span);
        text(ax, tip(1), tip(2), tip(3), sprintf('%d', k), 'Interpreter',S.interpreter, ...
            'FontSize',S.fontSizeSmall, 'Color',S.ink, ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle');
    end
end

% ---- freestream arrow: ahead of the leading edge near the root, pointing downstream ---------
% (outside the axis limits, hence unclipped; the label sits at its tail)
if opts.Freestream
    d_flow = -s_up*[1 0 0];
    base   = [xLE + s_up*1.0*c, 0.28*span, z0];
    local_arrow(ax, base, d_flow, 0.10*span, S.ink, 0.025, n_cam, [], false, [], 'both');
    p = base - d_flow*0.03*span;
    text(ax, p(1), p(2), p(3), '$V_\infty$', 'Interpreter',S.interpreter, ...
        'FontSize',S.fontSize, 'Color',S.ink, 'HorizontalAlignment','right', 'VerticalAlignment','middle');
end

% ---- frame triads --------------------------------------------------------
% The drawn frame has its origin at the plot origin: the root leading edge (aero) or the root
% end of the flexural axis (structural). With 'both', the other frame sits at x_other on the
% root edge with its x and z axes reversed (x_a = x_f c - x_s, z_a = -z_s), a fraction of a
% chord away, where two triads would overlap: it is drawn floated inboard of the root, with a
% dotted leader from its true origin. Labels beside the root edge (x), above the arrow (y) and
% beyond the tip (z); a label that falls on the planform gets a white background.
if ~strcmp(triads, 'none')
    Lt = 0.085*span;
    r_scr = cross(up, n_cam);  u_scr = cross(n_cam, r_scr);           % screen axes in the frame
    scr = @(p) [p*r_scr(:), p*u_scr(:)];
    planform = scr([xLE 0 z0; xTE 0 z0; xTE span z0; xLE span z0]);  % outline on screen
    origins = {[0 0 z0]};
    frames  = {eye(3)};
    labs    = {{'$x$', '$y$', '$z$'}};
    if strcmp(triads, 'both')
        if strcmp(opts.Frame, 'aero'), sub = {'a', 's'}; else, sub = {'s', 'a'}; end
        labs = {strcat('$', {'x','y','z'}, '_', sub{1}, '$'), strcat('$', {'x','y','z'}, '_', sub{2}, '$')};
        F = [x_other + s_up*0.45*c, -0.20*span, z0];                  % the other frame, floated inboard
        plot3(ax, [x_other F(1)], [0 F(2)], [z0 z0] + lift(3), ':', 'Color',S.muted, ...
            'LineWidth',S.lineWidthHair, 'Clipping','off');             % beyond the limits: no clipping
        origins{2} = F;
        frames{2}  = diag([-1 1 -1]);
    end
    back = 0.25*lift;                                                  % the casings sit behind the ink
    for f = 1:numel(origins)
        o = origins{f} + lift;
        e = frames{f};
        for i = 1:3
            local_arrow(ax, o, e(i,:), Lt, S.ink, 0.025, n_cam, S.white, f == 1, back, 'halo');
        end
        for i = 1:3
            local_arrow(ax, o, e(i,:), Lt, S.ink, 0.025, n_cam, S.white, f == 1, back, 'ink');
            tip = o + e(i,:)*Lt;
            switch i
                case 1, p = tip + [0, -0.05*span, 0];
                case 2, p = tip + 0.05*span*up;
                case 3, p = tip + 0.04*span*e(3,:);
            end
            q = scr(p);
            if inpolygon(q(1), q(2), planform(:,1), planform(:,2)), bg = S.white; else, bg = 'none'; end
            text(ax, p(1), p(2), p(3), labs{f}{i}, 'Interpreter',S.interpreter, ...
                'FontSize',S.fontSize, 'Color',S.ink, 'BackgroundColor',bg, 'Margin',1, ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle');
        end
    end
end

scene = struct('n_cam', n_cam, 'up', up, 'lift', lift, 'span', span, 'chord', c, ...
    'xLE', xLE, 'xTE', xTE, 'z0', z0);
end

% -------------------------------------------------------------------------
function idx = local_indices(v, n, name)
%LOCAL_INDICES  [] -> none, 'all' -> 1:n, numeric -> validated indices.
if isempty(v)
    idx = [];
elseif (ischar(v) || isstring(v)) && strcmpi(char(v), 'all')
    idx = 1:n;
elseif isnumeric(v)
    idx = reshape(v, 1, []);
    if any(idx < 1 | idx > n | idx ~= round(idx))
        error('wing_scene:index', '%s must be indices between 1 and %d.', name, n);
    end
else
    error('wing_scene:index', '%s must be [], ''all'' or a vector of indices.', name);
end
end

% -------------------------------------------------------------------------
function local_arrow(ax, base, dir, L, col, sw, n_cam, halo, clip, back, pass)
%LOCAL_ARROW  One filled polygon in the screen plane: a shaft of half width sw*L and a flat head
%   0.30 L long and 0.24 L wide. halo is a casing colour or []: the casing is the same polygon,
%   0.02 L wider all round, drawn the step 'back' away from the camera, so that an arm reads over an
%   ink line it runs along (the y arms lie on the flexural axis and the leading edge). pass 'halo',
%   'ink' or 'both' selects what to draw, so that a triad can lay all three casings before any ink.
%   clip false draws the arrow even outside the axis limits.
if nargin < 9  || isempty(clip), clip = true; end
if nargin < 10 || isempty(back), back = [0 0 0]; end
if nargin < 11, pass = 'both'; end
if clip, cl = 'on'; else, cl = 'off'; end
hh  = 0.30*L;                                        % head length
hw  = 0.12*L;                                        % head half width
w   = cross(n_cam, dir);                             % across the arrow, in the screen plane
if norm(w) < 1e-9, w = cross(n_cam, [0 0 1]); end   % arrow along the view direction: any width
w   = w/norm(w);
    function V = polygon(b, s_w, h_w, h_l, t_ext)
        t  = b + dir*(L + t_ext);                    % tip
        hb = t - dir*h_l;                            % head base
        V  = [b + s_w*w; hb + s_w*w; hb + h_w*w; t; hb - h_w*w; hb - s_w*w; b - s_w*w];
    end
if any(strcmp(pass, {'halo','both'})) && ~isempty(halo)
    m = 0.02*L;
    patch(ax, 'Vertices',polygon(base - back - dir*m, sw*L + m, hw + 1.3*m, hh + 1.5*m, 1.5*m), ...
        'Faces',1:7, 'FaceColor',halo, 'EdgeColor','none', 'Clipping',cl);
end
if any(strcmp(pass, {'ink','both'}))
    patch(ax, 'Vertices',polygon(base, sw*L, hw, hh, 0), 'Faces',1:7, ...
        'FaceColor',col, 'EdgeColor','none', 'Clipping',cl);
end
end
