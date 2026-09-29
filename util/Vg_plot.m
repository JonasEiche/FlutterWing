function [fig, crossings] = Vg_plot(V, EV, labels, opts)
%VG_PLOT  V-g / V-f flutter plot in the paper style, with the interpolated flutter crossing.
%   [fig, crossings] = Vg_plot(V, EV, labels)
%   [fig, crossings] = Vg_plot(V, EV, labels, Name=Value)
%   [fig, crossings] = Vg_plot({V1,V2,...}, {EV1,EV2,...}, {'open loop','closed loop'})
%   V        [1 x n_vel] velocities (m/s), or a cell {V1,V2,...} with one row vector per case
%   EV       [n x n_vel] tracked eigenvalues (rad/s) from getEigenvalueModeshape or
%            pkmethode, or a cell {EV1,EV2,...} per case. One row per mode, columns in the
%            order of V; conjugate pairs occupy two rows.
%   labels   one label per case (cellstr, string array, or a char for one case). Every text
%            is rendered by LaTeX, so escape specials in a label: 'flap\_1', '5 \%'.
%   Options (defaults)
%     'Band'          []        [f_lo f_hi] Hz: keep only the rows that oscillate at V(1)
%                               (imag > 0) with abs(EV(:,1))/(2*pi) inside the band; real
%                               poles have no frequency and are dropped. Conjugate rows
%                               (imag < 0) are always dropped, so every mode is drawn once.
%                               [] keeps all rows with imag(EV(:,1)) >= 0, real poles included.
%     'MaxDamping'    []        scalar in %: drop rows whose damping -100*Re/abs at V(1)
%                               exceeds it (with or without Band). Hides the well damped
%                               poles a controller adds to a closed loop; the shipped Goland
%                               controller has a pair at 8.4 Hz with 24 % damping.
%     'Interpolate'   true      report the linear zero crossing of Re(lambda) between grid
%                               points k-1 and k, k = first point with Re(lambda) > Tol;
%                               false reports grid point k itself, as Vg_plot_mult_bw does
%     'Tol' 0.1, 'DivTol' 1.0   (rad/s) grid rule threshold on Re; Im below DivTol at k
%                               means divergence instead of flutter
%     'DampingLimits' [-20 25]  y limits of the damping tile (%); 'FrequencyLimits' [] (Hz)
%     'Labels'        true      annotate the crossing next to its marker ('104.3 m/s, 4.52 Hz')
%     'Print'         true      one console line per case
%     'Size'          fw_style().size.standard   cm, for a newly created figure
%     'Parent'        []        figure or 2x1 tiledlayout to draw into; [] creates
%                               fw_figure(Size(1), Size(2), 'Name', Name)
%     'Name'          'V-g plot'  name of a newly created figure
%     'Colors'        []        [] = the fw_style colour order [blue; muted], with the line
%                               styles cycled after the colours as fw_figure does for a plain
%                               plot (case 1 blue, 2 muted, 3 blue dashed, 4 muted dashed, ...);
%                               an [m x 3] RGB matrix replaces the colour order
%     'Title'         'Frequency and damping over velocity'   title of the frequency tile
%     'LegendLocation' 'southwest'  legend Location in the damping tile (two cases or more)
%     'LineWidth'     fw_style().lineWidthThin   0.8 pt, the thin trace of the time histories
%   fig        figure with a 2x1 tiled layout in the style of the paper V-g figures
%              (research-paper-code/.../fig11_VgPlot_OL.m): title, frequency
%              abs(lambda)/(2*pi) on top, damping -100*Re(lambda)/abs(lambda) below, shared
%              x axis (no tick labels on the top tile), LaTeX text, grid and minor grid,
%              legend in the damping tile (two cases or more; a single case is drawn without
%              one), no box. Curves are thin: blue for the first case,
%              muted for the second, then the line styles cycle (dashed blue, dashed muted);
%              the first case is drawn last, so it stays visible where curves coincide.
%              The lowest crossing of each case is a hollow coral circle at (V,f) and
%              (V,0) for flutter, a hollow coral square at (V,0) in both tiles for divergence
%              (the paper script puts that one at -100 %, off the axis).
%   crossings  [1 x n_cases] struct, the lowest-velocity crossing of each case:
%              .label; .type 'flutter' | 'divergence' | 'none'; .V .f crossing velocity
%              (m/s) and frequency (Hz), interpolated when Interpolate is true and k > 1
%              (f = 0 for divergence); .V_grid .f_grid grid point k; .mode row of the
%              filtered EV; .index k; .interpolated logical. Numeric fields NaN for 'none'.
%   Band and MaxDamping select the rows of the plot and of the crossing search alike, so a
%   hidden mode is never reported as the flutter mechanism. The grid rule is the one of
%   Vg_plot_mult_bw: a pole that returns to the left half plane later is not reported;
%   refine V near the crossing when the interpolated value has to be accurate.
%   Example (RectWing: flutter at about 104.3 m/s, first grid point 105.8 m/s):
%     V_inf = linspace(20,160,32);
%     G     = build_G_RectWing(V_inf, [4,8], 4);
%     EV    = getEigenvalueModeshape(G, 5);
%     [fig, cr] = Vg_plot(V_inf, EV, 'open loop', 'Band', [2 10]);
%   Runtime: 0.25 s (1.9 s for the first figure of a session), measured on R2026a.
%   See also getEigenvalueModeshape, pkmethode, Vg_plot_mult_bw, fw_figure, fw_export.

arguments
    V
    EV
    labels
    opts.Band double = []
    opts.MaxDamping double = []
    opts.Interpolate (1,1) logical = true
    opts.Tol (1,1) double {mustBeReal} = 0.1
    opts.DivTol (1,1) double {mustBeReal} = 1.0
    opts.DampingLimits (1,2) double {mustBeReal} = [-20 25]
    opts.FrequencyLimits double = []
    opts.Labels (1,1) logical = true
    opts.Print (1,1) logical = true
    opts.Size (1,2) double {mustBePositive} = fw_style().size.standard
    opts.Parent = []
    opts.Name (1,:) char = 'V-g plot'
    opts.Colors double = []
    opts.LineWidth (1,1) double {mustBePositive} = fw_style().lineWidthThin
    opts.Title (1,:) char = 'Frequency and damping over velocity'
    opts.LegendLocation (1,:) char = 'southwest'
end

S = fw_style();

%% ------------------------------------------------------------------ inputs
if ~iscell(V),  V  = {V};  end
if ~iscell(EV), EV = {EV}; end
if ischar(labels)
    labels = {labels};
elseif isstring(labels)
    labels = cellstr(labels);
elseif ~iscell(labels)
    error('Vg_plot:labels', ...
        'labels must be a char, a string array or a cell array of labels, not %s.', class(labels));
end
V      = reshape(V, 1, []);
EV     = reshape(EV, 1, []);
labels = reshape(labels, 1, []);
n_cases = numel(V);
if numel(EV) ~= n_cases || numel(labels) ~= n_cases
    error('Vg_plot:caseCount', ...
        'V, EV and labels must describe the same number of cases (got %d, %d and %d).', ...
        n_cases, numel(EV), numel(labels));
end
for i = 1:n_cases
    if isstring(labels{i}) && isscalar(labels{i})
        labels{i} = char(labels{i});
    end
    if ~ischar(labels{i})
        error('Vg_plot:labels', 'labels{%d} must be a char or a scalar string.', i);
    end
    if ~isnumeric(V{i}) || ~isrow(V{i}) || ~isreal(V{i}) || isempty(V{i})
        error('Vg_plot:velocity', 'V{%d} must be a non-empty real row vector of velocities.', i);
    end
    if ~isnumeric(EV{i}) || ~ismatrix(EV{i}) || size(EV{i}, 2) ~= numel(V{i})
        error('Vg_plot:size', ...
            'EV{%d} must be an [n_modes x %d] matrix, one column per velocity (got %s).', ...
            i, numel(V{i}), mat2str(size(EV{i})));
    end
end
if ~isempty(opts.Band) && (~isequal(size(opts.Band), [1 2]) || ~all(isfinite(opts.Band)) ...
        || opts.Band(2) <= opts.Band(1))
    error('Vg_plot:band', 'Band must be [] or a 1x2 vector [f_lo f_hi] in Hz with f_hi > f_lo.');
end
if ~isempty(opts.MaxDamping) && ~isscalar(opts.MaxDamping)
    error('Vg_plot:maxDamping', 'MaxDamping must be [] or a scalar damping in %%.');
end
if ~isempty(opts.FrequencyLimits) && (~isequal(size(opts.FrequencyLimits), [1 2]) ...
        || opts.FrequencyLimits(2) <= opts.FrequencyLimits(1))
    error('Vg_plot:frequencyLimits', 'FrequencyLimits must be [] or an increasing 1x2 vector.');
end
if opts.DampingLimits(2) <= opts.DampingLimits(1)
    error('Vg_plot:dampingLimits', 'DampingLimits must be increasing.');
end
if ~isempty(opts.Colors) && size(opts.Colors, 2) ~= 3
    error('Vg_plot:colors', 'Colors must be [] or an [m x 3] matrix of RGB rows.');
end

%% -------------------------------------------- filter, damping and crossing
F = cell(1, n_cases);            % frequency  [n_kept x n_vel], Hz
D = cell(1, n_cases);            % damping    [n_kept x n_vel], %
crossings = repmat(struct('label', '', 'type', 'none', 'V', NaN, 'f', NaN, ...
    'V_grid', NaN, 'f_grid', NaN, 'mode', NaN, 'index', NaN, 'interpolated', false), 1, n_cases);

for i = 1:n_cases
    vel = V{i};
    ev  = EV{i};

    % keep one row per mode: drop the conjugates, then the band at V(1)
    f_0  = abs(ev(:,1))/(2*pi);
    keep = imag(ev(:,1)) >= 0;
    if ~isempty(opts.Band)
        % a band selects oscillatory modes: real poles (imag == 0) have no frequency
        keep = imag(ev(:,1)) > 0 & f_0 >= opts.Band(1) & f_0 <= opts.Band(2);
    end
    if ~isempty(opts.MaxDamping)
        % drop what is already well damped at V(1), e.g. the poles of a controller
        keep = keep & -100*real(ev(:,1))./abs(ev(:,1)) <= opts.MaxDamping;
    end
    ev = ev(keep,:);
    F{i} = abs(ev)/(2*pi);
    D{i} = -100*real(ev)./abs(ev);
    if isempty(ev)
        warning('Vg_plot:emptyBand', 'Case "%s": no eigenvalue row left after filtering.', labels{i});
    end

    % first grid point past zero damping of every kept row
    % (same as cumsum(real(ev) > Tol, 2) == 1 in Vg_plot_mult_bw)
    row = NaN; k = NaN; v_first = Inf;
    for r = 1:size(ev, 1)
        k_r = find(real(ev(r,:)) > opts.Tol, 1);
        if ~isempty(k_r) && vel(k_r) < v_first      % ties keep the lower row index
            v_first = vel(k_r);
            row     = r;
            k       = k_r;
        end
    end

    cr = crossings(i);
    cr.label = labels{i};
    if ~isnan(row)
        lam = ev(row,:);
        cr.V_grid = vel(k);
        cr.f_grid = abs(lam(k))/(2*pi);
        cr.mode   = row;
        cr.index  = k;
        cr.V      = cr.V_grid;
        cr.f      = cr.f_grid;
        if imag(lam(k)) > opts.DivTol
            cr.type = 'flutter';
        else
            cr.type = 'divergence';
        end
        if opts.Interpolate && k > 1
            re = real(lam);
            w  = (0 - re(k-1))/(re(k) - re(k-1));   % re(k-1) <= Tol < re(k), so re(k) > re(k-1)
            f_k = abs(lam((k-1):k))/(2*pi);
            cr.V = vel(k-1) + w*(vel(k) - vel(k-1));
            cr.f = f_k(1)   + w*(f_k(2) - f_k(1));
            cr.interpolated = true;
        end
        if strcmp(cr.type, 'divergence')
            cr.f = 0;
        end
    end
    crossings(i) = cr;

    if opts.Print
        switch cr.type
            case 'flutter'
                if opts.Interpolate
                    fprintf('%s: flutter at %.1f m/s, %.2f Hz (first grid point past zero damping %.1f m/s)\n', ...
                        cr.label, cr.V, cr.f, cr.V_grid);
                else
                    disp([num2str(cr.V_grid) ' m/s = flutter     ' num2str(cr.f_grid) ' Hz'])
                end
            case 'divergence'
                if opts.Interpolate
                    fprintf('%s: divergence at %.1f m/s (first grid point %.1f m/s)\n', ...
                        cr.label, cr.V, cr.V_grid);
                else
                    disp([num2str(cr.V_grid) ' m/s = divergence'])
                end
            otherwise
                if opts.Interpolate
                    fprintf('%s: no flutter or divergence up to %g m/s\n', cr.label, vel(end));
                else
                    disp('No flutter or divergence')
                end
        end
    end
end

%% ------------------------------------------------------------------ figure
if isempty(opts.Parent)
    fig = fw_figure(opts.Size(1), opts.Size(2), 'Name', opts.Name);
    tl  = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
elseif isa(opts.Parent, 'matlab.graphics.layout.TiledChartLayout')
    tl  = opts.Parent;                                  % used as given, create it 2 x 1
    fig = ancestor(tl, 'figure');
elseif isgraphics(opts.Parent)
    fig = ancestor(opts.Parent, 'figure');
    tl  = tiledlayout(opts.Parent, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
else
    error('Vg_plot:parent', 'Parent must be [], a figure handle or a tiledlayout handle.');
end
if isa(opts.Parent, 'matlab.ui.Figure')
    fig = opts.Parent;
end

ax_f = nexttile(tl, 1);            % frequency over velocity
ax_d = nexttile(tl, 2);            % damping over velocity
for ax = [ax_f ax_d]
    hold(ax, 'on');
    % set instead of grid(ax,'minor'), which toggles and would undo the fw_figure default
    set(ax, 'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'on', 'YMinorGrid', 'on');
end

x_lim = [min(cellfun(@min, V)) max(cellfun(@max, V))];
if x_lim(2) <= x_lim(1)
    x_lim = x_lim(1) + [-1 1];
end
x0 = x_lim; y0 = [0 0];
plot(ax_d, x0, y0, '--', 'Color', S.ink, 'LineWidth', S.lineWidthHair, 'HandleVisibility', 'off');   % zero damping line

% colour order [blue; muted] (or the colours the caller passed), line styles after the colours
cols = opts.Colors;
if isempty(cols), cols = S.cases; end
n_col   = size(cols, 1);
h_first = gobjects(1, n_cases);
for i = n_cases:-1:1                                       % first case last: on top where curves coincide
    style = {'Color', cols(mod(i-1, n_col) + 1, :), ...
             'LineStyle', S.lineStyles{mod(floor((i-1)/n_col), numel(S.lineStyles)) + 1}, ...
             'LineWidth', opts.LineWidth};
    if isempty(F{i})
        h_first(i) = plot(ax_d, NaN, NaN, style{:});       % keep the legend entry
        continue
    end
    plot(ax_f, V{i}, F{i}.', style{:});
    h = plot(ax_d, V{i}, D{i}.', style{:});
    h_first(i) = h(1);                                     % the legend lives in ax_d
end

% crossing markers on top of the curves, kept out of the legend
for i = 1:n_cases
    cr = crossings(i);
    if strcmp(cr.type, 'none')
        continue
    end
    if strcmp(cr.type, 'flutter')
        marker = {S.markerFlutter, 'MarkerEdgeColor', S.coral};
        y_f    = cr.f;
    else
        marker = {S.markerDivergence, 'MarkerEdgeColor', S.coral};
        y_f    = 0;                                        % paper: -100 % in ax_d, off the axis
    end
    marker = [marker, {'MarkerFaceColor', 'none', 'MarkerSize', S.markerSize, ...
        'LineWidth', S.lineWidthMarker, 'HandleVisibility', 'off'}];                                    %#ok<AGROW>
    plot(ax_f, cr.V, y_f, marker{:});
    plot(ax_d, cr.V, 0,   marker{:});
end

title(ax_f, opts.Title, 'Interpreter', S.interpreter);
ylabel(ax_f, 'Frequency (Hz)', 'Interpreter', S.interpreter);
ylabel(ax_d, 'Damping (\%)',   'Interpreter', S.interpreter);
xlabel(ax_d, '$V_\infty$ (m/s)', 'Interpreter', S.interpreter);
xlim(ax_f, x_lim);
xlim(ax_d, x_lim);
ylim(ax_d, opts.DampingLimits);
if ~isempty(opts.FrequencyLimits)
    ylim(ax_f, opts.FrequencyLimits);
end
set(ax_f, 'XTickLabel', []);                               % shared x axis, ticks below
if n_cases > 1                                             % a single case needs no legend
    lgd = legend(h_first, labels, 'Location', opts.LegendLocation, 'Interpreter', S.interpreter);
    lgd.AutoUpdate = 'off';
else
    lgd = [];
end

%% ------------------------------------------------------- crossing annotations
% below the zero damping line next to the marker: no kept mode crosses before cr.V, so
% that quadrant is free of curves. The 'best' legend ignores text, hence the candidate
% slots: take the first that stays inside the axes and misses the legend box.
if opts.Labels && any(~strcmp({crossings.type}, 'none'))
    lgd_px = zeros(0,4);
    if ~isempty(lgd)
        try
            lgd_px = getpixelposition(lgd, true);
        catch
        end
    end
    ax_px = getpixelposition(ax_d, true);
    y_l   = ylim(ax_d);
    d_x   = 0.02*diff(x_lim);
    d_y   = diff(opts.DampingLimits);
    n_ann = 0;
    for i = 1:n_cases
        cr = crossings(i);
        if strcmp(cr.type, 'none')
            continue
        end
        if strcmp(cr.type, 'flutter')
            str = sprintf('%.1f m/s, %.2f Hz', cr.V, cr.f);
        else
            str = sprintf('%.1f m/s, divergence', cr.V);
        end
        if cr.V - d_x > x_lim(1) + 0.3*diff(x_lim)          % left of the marker if it fits
            sides = {'right', 'left'};
        else
            sides = {'left', 'right'};
        end
        h_txt = text(ax_d, cr.V, 0, str, 'Color', S.ink, 'FontSize', S.fontSizeSmall, ...
            'VerticalAlignment', 'top', 'Interpreter', S.interpreter);
        placed  = false;
        n_ann   = n_ann + 1;
        fallback = {};                                      % first candidate inside the axes
        for slot = (n_ann-1) + (0:3)
            for i_s = 1:2
                pos = [i_x(cr.V, d_x, sides{i_s}), -(0.05 + slot*0.11)*d_y, 0];
                set(h_txt, 'Position', pos, 'HorizontalAlignment', sides{i_s});
                ext = get(h_txt, 'Extent');
                if ext(1) < x_lim(1) || ext(1) + ext(3) > x_lim(2) || ext(2) < y_l(1)
                    continue                                % would leave the axes
                end
                if isempty(fallback)
                    fallback = {pos, sides{i_s}};
                end
                if isempty(lgd_px) || ~i_overlaps(ext, ax_px, x_lim, y_l, lgd_px)
                    placed = true;
                    n_ann  = slot + 1;
                    break
                end
            end
            if placed
                break
            end
        end
        if ~placed && ~isempty(fallback)                    % every slot hits the legend
            set(h_txt, 'Position', fallback{1}, 'HorizontalAlignment', fallback{2});
        end
    end
end

end

% ---------------------------------------------------------------------------------
function x = i_x(v, d_x, side)
%I_X  x position of an annotation left ('right' aligned) or right ('left' aligned) of v.
if strcmp(side, 'right')
    x = v - d_x;
else
    x = v + d_x;
end
end

function tf = i_overlaps(ext, ax_px, x_l, y_l, box_px)
%I_OVERLAPS  true if the text extent ext (data units) overlaps box_px (figure pixels).
r  = [ax_px(1) + (ext(1) - x_l(1))/diff(x_l)*ax_px(3), ...
      ax_px(2) + (ext(2) - y_l(1))/diff(y_l)*ax_px(4), ...
      ext(3)/diff(x_l)*ax_px(3), ext(4)/diff(y_l)*ax_px(4)];
tf = r(1) < box_px(1) + box_px(3) && r(1) + r(3) > box_px(1) && ...
     r(2) < box_px(2) + box_px(4) && r(2) + r(4) > box_px(2);
end
