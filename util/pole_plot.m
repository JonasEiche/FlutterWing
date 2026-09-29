function fig = pole_plot(P, labels, opts)
%POLE_PLOT  Pole map in the s-plane: hollow blue open loop, larger hollow muted closed loop, coral for a pole in the right half plane.
%   fig = pole_plot(sys, 'open loop')                 poles of one LTI model (ss, tf, zpk)
%   fig = pole_plot({sys1, sys2}, {'open loop','closed loop'})
%                                                     several cases: the 1st is drawn as hollow blue
%                                                     circles, the 2nd as larger hollow muted circles
%                                                     (a pole both cases share shows both rings), the
%                                                     3rd as blue dots, the 4th as muted dots; a pole
%                                                     with Re > 0 is coral whatever its case
%   fig = pole_plot(G, 'open loop', 'Velocity', V_inf)
%                                                     ss array gridded over V_inf (build_G_*): the
%                                                     poles of every grid point as muted dots, the
%                                                     last grid point with the case marker
%   fig = pole_plot(EV, 'open loop', 'Velocity', V_inf)
%                                                     tracked eigenvalues [n_modes x n_vel] in rad/s
%                                                     (getEigenvalueModeshape, pkmethode): one thin
%                                                     blue locus per mode, coral where Re > 0, the
%                                                     case marker at V_inf(end)
%   P        LTI model, ss array, EV matrix, or a cell with one of these per case
%   labels   one label per case (char, string array or cellstr); LaTeX text, escape specials
%   Options (defaults)
%     'Velocity'  []     velocities (m/s) of an ss array or EV matrix, one per grid point
%     'Limits'    []     [re_lo re_hi im_lo im_hi] in rad/s; [] pads the poles by 8 %, makes the
%                        imaginary range symmetric and keeps some right half plane in view
%     'Size'      fw_style().size.square   cm, for a new figure
%     'Parent'    []     figure or axes to draw into; [] creates fw_figure(Size, 'Name', Name)
%     'Name'      'pole map'   name of a newly created figure
%     'Title'     ''     title (LaTeX)
%     'LegendLocation' 'best'   legend Location; the tutorial passes 'southwest', which is empty
%                        in its maps and cannot be read as a pole
%   A dashed ink line marks the imaginary axis; a pole to its right is the only coral in the
%   figure. Replaces pzmap in the tutorial (pzmap is a chart object that ignores the fw_figure
%   text defaults).
%   Example:
%     V_inf = linspace(20,160,32);  G = build_G_RectWing(V_inf, [4,8], 4);
%     pole_plot(G, 'open loop', 'Velocity', V_inf, 'Limits', [-40 10 -150 150]);
%   See also Vg_plot, fw_figure, fw_style, getEigenvalueModeshape.

arguments
    P
    labels
    opts.Velocity double = []
    opts.Limits double = []
    opts.Size (1,2) double {mustBePositive} = fw_style().size.square
    opts.Parent = []
    opts.Name (1,:) char = 'pole map'
    opts.Title (1,:) char = ''
    opts.LegendLocation (1,:) char = 'best'
end

S = fw_style();
if ~iscell(P), P = {P}; end
if ischar(labels)
    labels = {labels};
elseif isstring(labels)
    labels = cellstr(labels);
end
P = reshape(P, 1, []);
labels = reshape(labels, 1, []);
n = numel(P);
if numel(labels) ~= n
    error('pole_plot:labels', 'One label per case: got %d cases and %d labels.', n, numel(labels));
end
V = reshape(opts.Velocity, 1, []);
if ~isempty(opts.Limits) && (~isequal(size(opts.Limits), [1 4]) || opts.Limits(2) <= opts.Limits(1) ...
        || opts.Limits(4) <= opts.Limits(3))
    error('pole_plot:limits', 'Limits must be [] or [re_lo re_hi im_lo im_hi] with increasing pairs.');
end

%% -------------------------------------------------------------- poles per case
poles = cell(1, n);          % column vector (one model) or [n_poles x n_vel] (over velocity)
track = false(1, n);         % true for tracked EV rows: draw loci
for i = 1:n
    Pi = P{i};
    if isnumeric(Pi)
        if isempty(V) || size(Pi, 2) ~= numel(V)
            error('pole_plot:velocity', ['Case %d: an EV matrix needs ''Velocity'' with one entry ' ...
                'per column (got %d columns and %d velocities).'], i, size(Pi, 2), numel(V));
        end
        poles{i} = Pi;
        track(i) = true;
    elseif isa(Pi, 'DynamicSystem')
        sz = size(Pi);                                   % [ny nu n3 ...] of an LTI array
        n_models = prod([sz(3:end) 1]);
        if n_models > 1
            if isempty(V) || n_models ~= numel(V)
                error('pole_plot:velocity', ['Case %d: an ss array needs ''Velocity'' with one entry ' ...
                    'per model (got %d models and %d velocities).'], i, n_models, numel(V));
            end
            p1 = pole(Pi(:,:,1));
            M  = NaN(numel(p1), n_models);
            for k = 1:n_models
                pk = pole(Pi(:,:,k));
                M(1:numel(pk), k) = pk;
            end
            poles{i} = M;
        else
            poles{i} = pole(Pi);
        end
    else
        error('pole_plot:input', 'Case %d must be an LTI model, an ss array or a numeric EV matrix.', i);
    end
end

%% ------------------------------------------------------------------ figure
if isempty(opts.Parent)
    fig = fw_figure(opts.Size(1), opts.Size(2), 'Name', opts.Name);
    ax  = axes(fig);
elseif isgraphics(opts.Parent, 'figure')
    fig = opts.Parent;
    ax  = axes(fig);
elseif isgraphics(opts.Parent) && isprop(opts.Parent, 'XLim')
    ax  = opts.Parent;
    fig = ancestor(ax, 'figure');
else
    error('pole_plot:parent', 'Parent must be [], a figure handle or an axes handle.');
end
hold(ax, 'on');
set(ax, 'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'Layer', 'top', ...
    'Box', 'on', 'TickLabelInterpreter', S.interpreter);

% limits: pad, symmetric imaginary range, right half plane in view
all_p = cellfun(@(m) m(:), poles, 'UniformOutput', false);
all_p = vertcat(all_p{:});
all_p = all_p(isfinite(all_p));
if isempty(opts.Limits)
    re  = [min(real(all_p)) max(real(all_p))];
    d_r = max(diff(re), 1);
    re  = re + 0.08*d_r*[-1 1];
    re(2) = max(re(2), 0.12*d_r);
    im  = max(abs(imag(all_p)));
    im  = max(im*1.08, 1);
    lims = [re, -im, im];
else
    lims = opts.Limits;
end
xlim(ax, lims(1:2));
ylim(ax, lims(3:4));

% the dashed imaginary axis: stability boundary
plot(ax, [0 0], lims(3:4), '--', 'Color', S.ink, 'LineWidth', S.lineWidthHair, 'HandleVisibility', 'off');

%% -------------------------------------------------------------------- cases
% marker per case: hollow blue, larger hollow muted, blue dot, muted dot (then repeat)
mk = {{'o', 'MarkerFaceColor', 'none',  'MarkerEdgeColor', S.blue,  'LineWidth', 1.5, 'MarkerSize', S.markerSize}, ...
      {'o', 'MarkerFaceColor', 'none',  'MarkerEdgeColor', S.muted, 'LineWidth', 1.2, 'MarkerSize', S.markerSize + 3}, ...
      {'o', 'MarkerFaceColor', S.blue,  'MarkerEdgeColor', S.blue,  'LineWidth', 0.5, 'MarkerSize', S.markerSize + 1}, ...
      {'o', 'MarkerFaceColor', S.muted, 'MarkerEdgeColor', S.muted, 'LineWidth', 0.5, 'MarkerSize', S.markerSize + 1}};
h_leg = gobjects(1, n);
for i = 1:n
    st = mk{mod(i-1, numel(mk)) + 1};
    M  = poles{i};
    if track(i)                                            % loci of tracked modes
        for r = 1:size(M, 1)
            plot(ax, real(M(r,:)), imag(M(r,:)), '-', 'Color', S.blue, ...
                'LineWidth', S.lineWidthThin, 'HandleVisibility', 'off');
            u = real(M(r,:)) > 0;
            if any(u)
                lam = M(r,:); lam(~u) = NaN;               % coral over the unstable part only
                plot(ax, real(lam), imag(lam), '-', 'Color', S.coral, ...
                    'LineWidth', S.lineWidthThin, 'HandleVisibility', 'off');
            end
        end
        p_end = M(:, end);
    elseif size(M, 2) > 1                                  % ss array: dots at every grid point
        M_in = M(:, 1:end-1);
        u = real(M_in) > 0;
        plot(ax, real(M_in(~u)), imag(M_in(~u)), '.', 'Color', S.muted, 'MarkerSize', 5, ...
            'HandleVisibility', 'off');
        if any(u(:))
            plot(ax, real(M_in(u)), imag(M_in(u)), '.', 'Color', S.coral, 'MarkerSize', 5, ...
                'HandleVisibility', 'off');
        end
        p_end = M(:, end);
    else
        p_end = M(:);
    end
    p_end  = p_end(isfinite(p_end));
    stable = real(p_end) <= 0;
    if any(stable)
        h_leg(i) = plot(ax, real(p_end(stable)), imag(p_end(stable)), st{:});
    else
        h_leg(i) = plot(ax, NaN, NaN, st{:});              % legend entry only
    end
    if any(~stable)
        st_c = local_coral(st, S);
        plot(ax, real(p_end(~stable)), imag(p_end(~stable)), st_c{:}, 'HandleVisibility', 'off');
    end
end

xlabel(ax, 'Re $\lambda$ (rad/s)', 'Interpreter', S.interpreter);
ylabel(ax, 'Im $\lambda$ (rad/s)', 'Interpreter', S.interpreter);
if ~isempty(opts.Title)
    title(ax, opts.Title, 'Interpreter', S.interpreter);
end
lgd = legend(ax, h_leg, labels, 'Location', opts.LegendLocation, 'Interpreter', S.interpreter);
lgd.AutoUpdate = 'off';
end

% ---------------------------------------------------------------------------------
function st = local_coral(st, S)
%LOCAL_CORAL  the same marker with every blue or muted colour replaced by coral.
for k = 1:numel(st)
    if isnumeric(st{k}) && (isequal(st{k}, S.blue) || isequal(st{k}, S.muted))
        st{k} = S.coral;
    end
end
end
