function fig = plot_aic_entries(k_red, Qjj, entries, opts)
%PLOT_AIC_ENTRIES  Entries of the AIC matrix Qjj(k) in the complex plane over the reduced frequency (TUTORIAL.m section (3)).
%   fig = plot_aic_entries(k_red, Qjj, entries)
%   fig = plot_aic_entries(k_red, Qjj, [1 1; 1 25], 'Name', 'aic_nyquist')
%   k_red    [1 x nk] reduced frequencies of the samples
%   Qjj      [n x n x nk] AIC samples from build_Qjj
%   entries  [m x 2] (row, column) pairs, one tile per pair, side by side
%   Each tile shows the samples as hollow muted circles joined by a thin line, with equal axis
%   scaling, and labels the samples nearest to k = 0, 0.1, 0.5 and max(k_red). The label
%   positions are tuned for the entries (1,1) and (1,25) of the Goland grid; other entries put
%   every label above right of its sample.
%   Options (defaults)
%     'Name'   'AIC entries'   figure name (the tutorial export uses the file stem)
%   fig      fw_figure of the standard size
%   See also build_Qjj, mimo_nyquist.

arguments
    k_red (1,:) double
    Qjj double
    entries (:,2) double
    opts.Name (1,:) char = 'AIC entries'
end

S = fw_style();
[~, i_klabel] = min(abs(k_red(:) - [0 0.1 0.5 max(k_red)]), [], 1);   % 1 7 11 14 for the tutorial grid
place_tuned = {[1 1],  {'center','top',    0, -1;          % Q(1,1):  k = 0   below the sample
                        'right','middle', -1,  0;          %          k = 0.1 left of it
                        'left','middle',   1,  0;          %          k = 0.5 right of it
                        'left','bottom',   1,  0.4};       %          k = 1.1 above right
               [1 25], {'right','bottom', -1,  0.5;        % Q(1,25): the curve hooks back, so k = 0 goes up left
                        'center','top',    0, -1;
                        'left','middle',   1,  0;
                        'left','bottom',   1,  0.4}};
place_default = repmat({'left','bottom', 1, 0.4}, numel(i_klabel), 1);

fig = fw_figure(S.size.standard(1), S.size.standard(2), 'Name', opts.Name);
tl  = tiledlayout(fig,1,size(entries,1),'TileSpacing','compact','Padding','compact');
for t = 1:size(entries,1)
    place = place_default;
    for p = 1:size(place_tuned,1)
        if isequal(entries(t,:), place_tuned{p,1}), place = place_tuned{p,2}; end
    end
    ax = nexttile(tl); hold(ax,'on')
    q_k = squeeze(Qjj(entries(t,1),entries(t,2),:));
    plot(ax, real(q_k), imag(q_k), '-o', 'Color',S.muted, 'LineWidth',S.lineWidthThin, ...
        'MarkerEdgeColor',S.muted, 'MarkerFaceColor',S.white, 'MarkerSize',S.markerSize)
    axis(ax,'equal')
    xl = xlim(ax) + 0.16*diff(xlim(ax))*[-1 1];     % room around the curve for the four labels
    yl = ylim(ax) + 0.16*diff(ylim(ax))*[-1 1];
    xlim(ax,xl); ylim(ax,yl)
    for m = 1:numel(i_klabel)
        text(ax, real(q_k(i_klabel(m))) + place{m,3}*0.035*diff(xl), ...
                 imag(q_k(i_klabel(m))) + place{m,4}*0.035*diff(yl), ...
            sprintf('$k = %g$', k_red(i_klabel(m))), ...
            'HorizontalAlignment',place{m,1}, 'VerticalAlignment',place{m,2}, ...
            'FontSize',S.fontSizeSmall, 'Color',S.ink)
    end
    xlabel(ax,'Re'); ylabel(ax,'Im')
    title(ax, sprintf('$Q_{jj}(%d,%d)$', entries(t,1), entries(t,2)))
end
end
