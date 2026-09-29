function fig = plot_wing_layout(Structure, imuIDX, ailIDX, opts)
%PLOT_WING_LAYOUT  Plan view of a wing model: panel mesh, numbered control surfaces, hinge lines, IMUs.
%   fig = plot_wing_layout(Structure, imuIDX, ailIDX)
%   fig = plot_wing_layout(Structure, imuIDX, ailIDX, 'Parent', ax, 'Size', [14 6], 'Name', 'wing layout', 'Colors', fw_style())
%   Structure  from define_RectWing_Structure_Aero or define_Goland_Structure_Aero;
%              only .Ps (panel corners, structural coordinates) and .cspanels
%              (panel numbers of each control surface) are used
%   imuIDX     IMUs in use (as passed to build_G_* / build_P_*): drawn as blue
%              dots, the others as hollow ink circles
%   ailIDX     control surfaces in use: drawn blue, the others grey
%   'Parent'   [] (default) opens a new fw_figure of size 'Size' (cm) named
%              'Name'; a figure handle draws into a new axes of that figure, an
%              axes handle (e.g. from nexttile) into that axes. A parent axes
%              narrower than about 10 cm crowds the fixed-size labels.
%   'Colors'   fw_style() struct (fields ink, blue, grey, white, muted,
%              lineWidth*, fontSizeSmall are read)
%   fig        the figure drawn into
%
%   The span runs along the horizontal axis and the structural chord coordinate
%   x_s up the vertical axis. Structural x points upstream, so the leading edge
%   is on top, the trailing edge at the bottom, and the freestream arrow left of
%   the root points down. Two small triads mark the aerodynamic frame (origin
%   at the root leading edge, x_a downstream) and the structural frame (origin
%   on the flexural axis at the root, x_s toward the leading edge).
%
%   Numbering (from the define_* files and their ASCII diagrams): RectWing flaps
%   1 to 4 on the trailing edge, slats 5 to 8 on the leading edge, both root to
%   tip (ailIDX 1..4 = flap1_d..flap4_d, 5..8 = slat1_d..slat4_d); IMU k sits at
%   the hinge-line centre of surface k, so imuIDX and ailIDX share one
%   numbering. Goland: one trailing-edge flap 1 on the outer 30 % of the span.
%
%   Example:
%     [Structure, ~] = define_RectWing_Structure_Aero(5, 6);
%     plot_wing_layout(Structure, [4 8], 4);      % IMUs 4 and 8, flap 4 highlighted
%   Runtime: about 1.6 s (LaTeX text rendering), measured on R2026a.

arguments
    Structure (1,1) struct
    imuIDX    (1,:) double
    ailIDX    (1,:) double
    opts.Parent = []
    opts.Size   (1,2) double {mustBePositive} = [14 6]
    opts.Name   (1,:) char = 'wing layout'
    opts.Colors (1,1) struct = fw_style()
end

assert(isfield(Structure,'Ps') && isfield(Structure,'cspanels'), ...
    'plot_wing_layout:fields', ...
    'Structure needs the fields Ps and cspanels (from define_*_Structure_Aero).')

S  = opts.Colors;
Ps = Structure.Ps;
cs = Structure.cspanels;
n_pan = numel(Ps);
n_cs  = numel(cs);

% ---- geometry (structural coordinates: x upstream, LE = max x) ----------
xc_p = zeros(4,n_pan);      % corner x of every panel, node order 1..4
yc_p = zeros(4,n_pan);      %          y              (1 LE-in, 2 TE-in, 3 TE-out, 4 LE-out)
for i = 1:n_pan
    for n = 1:4
        xc_p(n,i) = Ps{i}{n}(1);
        yc_p(n,i) = Ps{i}{n}(2);
    end
end
xLE  = max(xc_p(:));        % leading edge  (structural x points upstream)
xTE  = min(xc_p(:));        % trailing edge
span = max(yc_p(:));        % semi span
c    = xLE - xTE;           % chord
tol  = 1e-6;

% classify each control surface by the chordwise position of its panel block
isTE = false(1,n_cs);
isLE = false(1,n_cs);
for k = 1:n_cs
    blk = cs{k};
    isTE(k) = abs(min(min(xc_p(:,blk))) - xTE) < tol;
    isLE(k) = abs(max(max(xc_p(:,blk))) - xLE) < tol;
end
isLE = isLE & ~isTE;        % a full-chord block counts as a trailing-edge surface

% ---- figure / axes ------------------------------------------------------
if isempty(opts.Parent)
    fig = fw_figure(opts.Size(1), opts.Size(2), 'Name', opts.Name);
    ax  = axes(fig);
    set(ax,'Position',[0.02 0.02 0.96 0.96]);
elseif isgraphics(opts.Parent,'figure')
    fig = opts.Parent;
    ax  = axes(fig);
elseif isgraphics(opts.Parent) && isprop(opts.Parent,'XLim')
    ax  = opts.Parent;
    fig = ancestor(ax,'figure');
else
    error('plot_wing_layout:parent','''Parent'' must be a figure or axes handle.')
end
hold(ax,'on')
set(ax,'Color',S.white)

edgeCol = S.ink;                             % panel edges: ink hairlines at low alpha

% ---- panel mesh (one patch with n_pan faces) ---------------------------
V = [yc_p(:), xc_p(:)];                     % vertices, 4 per panel
Fa = reshape(1:4*n_pan, 4, n_pan)';         % faces, corner order 1..4 of each panel
patch(ax, 'Faces',Fa, 'Vertices',V, 'FaceColor',S.white, ...
    'EdgeColor',edgeCol, 'EdgeAlpha',0.35, 'LineWidth',S.lineWidthMesh);

% ---- control-surface bands and their numbers ---------------------------
for k = 1:n_cs
    blk = cs{k};
    on  = ismember(k, ailIDX);
    if on, faceCol = S.blue; else, faceCol = S.grey; end
    patch(ax, 'Faces',Fa(blk,:), 'Vertices',V, 'FaceColor',faceCol, ...
        'EdgeColor',edgeCol, 'EdgeAlpha',0.35, 'LineWidth',S.lineWidthMesh);
    y_band = mean(mean(yc_p(:,blk)));
    x_band = mean(mean(xc_p(:,blk)));
    % the bands are too thin to hold the number: put it in the wing interior
    % next to its band (trailing-edge surfaces above, leading-edge below)
    if isLE(k), x_num = x_band - 0.22*c; else, x_num = x_band + 0.22*c; end
    text(ax, y_band, x_num, sprintf('%d',k), 'Interpreter','latex', ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
        'FontSize',S.fontSizeSmall, 'Color',S.ink);
end

% ---- hinge lines and IMU markers ---------------------------------------
for k = 1:n_cs
    blk = cs{k};
    if isLE(k)
        h1 = Ps{blk(end)}{3};   h2 = Ps{blk(1)}{2};     % slat hinge (aft edge)
    else
        h1 = Ps{blk(1)}{1};     h2 = Ps{blk(end)}{4};   % flap hinge (forward edge)
    end
    plot(ax, [h1(2) h2(2)], [h1(1) h2(1)], '-', 'Color', S.ink, 'LineWidth', S.lineWidth);
    HP = 0.5*(h1 + h2);                                 % IMU k at the hinge centre
    if ismember(k, imuIDX)
        plot(ax, HP(2), HP(1), 'o', 'MarkerSize',6, ...
            'MarkerFaceColor',S.blue, 'MarkerEdgeColor',S.blue, 'LineWidth',S.lineWidthThin);
    else
        plot(ax, HP(2), HP(1), 'o', 'MarkerSize',6, ...
            'MarkerFaceColor',S.white, 'MarkerEdgeColor',S.ink, 'LineWidth',S.lineWidthThin);
    end
end

% ---- freestream arrow left of the root (air flows LE -> TE, downward) ---
y_arr  = -0.12*span;
x_tail = xLE + 0.62*c;
x_head = xLE + 0.10*c;
hl     = 0.13*c;                                        % arrow head length
hw     = 0.045*c;                                       % arrow head half width
plot(ax, [y_arr y_arr], [x_tail x_head+0.9*hl], '-', 'Color',S.ink, 'LineWidth',S.lineWidth);
patch(ax, y_arr + [0 -hw hw], [x_head x_head+hl x_head+hl], S.ink, 'EdgeColor','none');
text(ax, y_arr - 0.025*span, 0.5*(x_tail+x_head), '$V_\infty$', 'Interpreter','latex', ...
    'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
    'FontSize',S.fontSizeSmall, 'Color',S.ink);

% ---- axis triads: aerodynamic frame at the root LE, structural at the root
% 0.08*span, but short enough that x_a (down from the LE) and x_s (up from the
% flexural axis) keep a visible gap on a wing whose flexural axis sits close to
% the leading edge
L = min(0.08*span, 0.32*max(xLE, 0.1*c));
local_arrow(ax, 0, xLE, 0, -L, S.muted);                % x_a: downstream (down)
local_arrow(ax, 0, xLE, L,  0, S.muted);                % y_a: spanwise
local_arrow(ax, 0, 0,   0,  L, S.muted);                % x_s: upstream (up)
local_arrow(ax, 0, 0,   L,  0, S.muted);                % y_s: spanwise
local_label(ax, -0.025*span,  xLE-0.55*L,   '$x_a$', 'right', 'middle', S);
local_label(ax,  L+0.02*span, xLE+0.02*c,   '$y_a$', 'left',  'bottom', S);
local_label(ax,  0.025*span,  0.55*L,       '$x_s$', 'left',  'middle', S);
local_label(ax,  L+0.02*span, 0,            '$y_s$', 'left',  'middle', S);

% ---- edge captions -----------------------------------------------------
text(ax, 0.5*span, xLE, local_edge_str('leading edge','slat',find(isLE)), ...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
    'FontSize',S.fontSizeSmall, 'Color',S.ink);
text(ax, 0.5*span, xTE, local_edge_str('trailing edge','flap',find(isTE)), ...
    'HorizontalAlignment','center', 'VerticalAlignment','top', ...
    'FontSize',S.fontSizeSmall, 'Color',S.ink);

% ---- limits: room for the arrow (left), the triads and the captions ----
grid(ax,'off')
axis(ax,'equal')
xlim(ax, [-0.30*span, 1.05*span])
ylim(ax, [xTE-0.35*c, xLE+0.78*c])
axis(ax,'off')
end

% -------------------------------------------------------------------------
function local_arrow(ax, y0, x0, dy, dx, col)
% straight arrow from (y0,x0) by (dy,dx), plain primitives only
len = hypot(dy,dx);
hl  = 0.30*len;                 % head length
hw  = 0.10*len;                 % head half width
uy  = dy/len;  ux = dx/len;     % unit vector along the arrow
py  = -ux;     px = uy;         % unit vector across it
y_b = y0 + dy - hl*uy;  x_b = x0 + dx - hl*ux;      % head base
plot(ax, [y0 y_b], [x0 x_b], '-', 'Color',col, 'LineWidth',0.8);
patch(ax, [y0+dy, y_b+hw*py, y_b-hw*py], [x0+dx, x_b+hw*px, x_b-hw*px], ...
    col, 'EdgeColor','none');
end

% -------------------------------------------------------------------------
function local_label(ax, y, x, str, hAlign, vAlign, S)
text(ax, y, x, str, 'Interpreter','latex', 'FontSize',S.fontSizeSmall, 'Color',S.muted, ...
    'HorizontalAlignment',hAlign, 'VerticalAlignment',vAlign);
end

% -------------------------------------------------------------------------
function str = local_edge_str(edge, noun, idx)
% 'trailing edge: flaps 1-4' / 'trailing edge: flap 1' / 'leading edge'
if isempty(idx)
    str = edge;
    return
end
idx = sort(idx(:))';
if isscalar(idx)
    str = sprintf('%s: %s %d', edge, noun, idx);
elseif all(diff(idx) == 1)
    str = sprintf('%s: %ss %d-%d', edge, noun, idx(1), idx(end));
else
    list = strjoin(arrayfun(@(v) sprintf('%d',v), idx, 'UniformOutput',false), ', ');
    str  = sprintf('%s: %ss %s', edge, noun, list);
end
end
