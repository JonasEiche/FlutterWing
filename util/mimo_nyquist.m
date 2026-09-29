function mimo_nyquist(wi,Qjj,H_hat_cell,legends,rows,columns,varargin)
%MIMO_NYQUIST  Nyquist-style comparison of AIC samples Qjj(k) with RFA fits H_hat(k), one tile per entry.
%   mimo_nyquist(k_red,Qjj,H_hat_cell,legends,rows,columns)
%   mimo_nyquist(k_red,Qjj,H_hat_cell,legends,rows,columns,'Parent',fig,'Colors',[navy; muted],'LegendLocation','southeast')
%   k_red       [1 x nk]        reduced frequencies (only min/max are used, in
%                               the figure name)
%   Qjj         [ny x nu x nk]  AIC samples, drawn as hollow muted circles in
%                               the complex plane
%   H_hat_cell  cell of [ny x nu x nk] fits (or a single array), drawn as lines
%               in the colour order with one line style per fit
%               (fw_style().lineStyles, cycled)
%   legends     cell of labels: either one per fit, or 1+numel(H_hat_cell)
%               labels with the first one for the Qjj samples; LaTeX text
%   rows, columns  panel indices to show; one tile per (row, column) pair,
%               titled $Q_{jj}(p,q)$, equal axis scaling
%   'Parent'    [] (default) opens fw_figure(fw_style().size.standard) named
%               Name; a figure handle draws the tiles into that figure
%   'Name'      '' (default): name of a newly created figure; '' names it after
%               the k_red range
%   'Colors'    [] (default) cycles fw_style().cases over the fits (blue,
%               muted); a cell array of colours or an [n x 3] RGB matrix
%               replaces it. The line styles cycle after the colours, as in fw_figure
%               and Vg_plot (fit 3 is blue dashed), and the first fit is drawn last so
%               that it stays visible where the fits coincide.
%   'LegendLocation' 'best' (default): legend Location in the first tile; 'best' can
%               land on a curve after axis equal
%   'LegendTile' '' (default): the legend lives in the first tile; 'north' | 'south' |
%               'east' | 'west' moves it into that slot of the tiled layout, outside the
%               tiles (horizontal for north and south), which is what the tutorial does
%   Example (TUTORIAL.m section (5)):
%     mimo_nyquist(k_red,Qjj,{H_rog},{'DLM samples','Roger fit'},[1 25],[1 25])
%   See also rogersRFA_magW, evalRFA, fw_figure, fw_style.

    if ~iscell(H_hat_cell)
        H_hat_cell={H_hat_cell};
    end

    parser = inputParser;
    parser.FunctionName = 'mimo_nyquist';
    addParameter(parser,'Parent',[]);
    addParameter(parser,'Name','');
    addParameter(parser,'Colors',[]);
    addParameter(parser,'LegendLocation','best');
    addParameter(parser,'LegendTile','');
    parse(parser,varargin{:});
    parent    = parser.Results.Parent;
    fig_name  = char(parser.Results.Name);
    colors_in = parser.Results.Colors;
    lgd_loc   = parser.Results.LegendLocation;
    lgd_tile  = char(parser.Results.LegendTile);
    if ~isempty(lgd_tile)
        assert(any(strcmp(lgd_tile, {'north','south','east','west'})), 'mimo_nyquist:legendTile', ...
            'LegendTile must be '''', ''north'', ''south'', ''east'' or ''west''.')
    end
    S = fw_style();

    if isempty(fig_name)
        fig_name = sprintf('Frequency response H(ik), %g <= k_red <= %g',min(wi), max(wi));
    end
    if isempty(parent)
        hFig = fw_figure(S.size.standard(1), S.size.standard(2), 'Name', fig_name);
    else
        assert(isgraphics(parent,'figure'),'mimo_nyquist:parent','''Parent'' must be a figure handle.')
        hFig = parent;                               % draw into the given figure, keep its defaults
    end

    if isempty(colors_in)
        colors = num2cell(S.cases, 2)';
    elseif iscell(colors_in)
        colors = colors_in;
    else
        assert(size(colors_in,2) == 3,'mimo_nyquist:colors','''Colors'' must be a cell array or an [n x 3] RGB matrix.')
        colors = num2cell(colors_in,2)';
    end

    np = length(rows);
    nq = length(columns);
    tl = tiledlayout(hFig, np, nq, 'TileSpacing','compact', 'Padding','compact');
    for ip = 1:np
        p = rows(ip);
        for iq = 1:nq
            q  = columns(iq);
            ax = nexttile(tl);
            hold(ax, 'on');
            set(ax, 'XGrid','on', 'YGrid','on', 'XMinorGrid','on', 'YMinorGrid','on', ...
                'Box','on', 'TickLabelInterpreter',S.interpreter);
            z = squeeze(Qjj(p,q,:));
            hDot = plot(ax, real(z), imag(z), 'o', 'MarkerSize',S.markerSize, ...
                'MarkerEdgeColor',S.muted, 'MarkerFaceColor',S.white, 'LineWidth',S.lineWidthThin);

            hCurve = gobjects(numel(H_hat_cell),1);
            for h = numel(H_hat_cell):-1:1                   % first fit last: on top where they coincide
                H_hat = H_hat_cell{h};
                zh = squeeze(H_hat(p,q,:));
                hCurve(h) = plot(ax, real(zh), imag(zh), ...
                    'LineStyle', S.lineStyles{mod(floor((h-1)/numel(colors)),numel(S.lineStyles))+1}, ...
                    'Color', colors{mod(h-1,numel(colors))+1}, 'LineWidth', S.lineWidth);
            end
            axis(ax, 'equal')
            xlabel(ax, 'Re', 'Interpreter',S.interpreter)
            ylabel(ax, 'Im', 'Interpreter',S.interpreter)
            title(ax, sprintf('$Q_{jj}(%d,%d)$', p, q), 'Interpreter',S.interpreter)
            if ip==1 && iq==1
                if numel(legends) == numel(H_hat_cell)
                    lgd = legend(ax, hCurve, legends);          % labels for the fits only
                elseif numel(legends) == numel(H_hat_cell)+1
                    lgd = legend(ax, [hDot; hCurve], legends);  % first label for the Qjj samples
                else
                    lgd = legend(ax, legends);
                end
                set(lgd, 'Interpreter',S.interpreter, 'Location',lgd_loc);
                if ~isempty(lgd_tile)                       % out of the tiles, into the layout
                    lgd.Layout.Tile = lgd_tile;
                    if any(strcmp(lgd_tile, {'north','south'})), lgd.Orientation = 'horizontal'; end
                end
            end
        end
    end
end
