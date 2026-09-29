function fig = plot_actuator_bode(G_acts, labels, f_band, opts)
%PLOT_ACTUATOR_BODE  Bode plot of actuator models from demand to deflection, with the flutter band shaded (TUTORIAL.m section (8)).
%   fig = plot_actuator_bode(G_acts, labels, f_band)
%   fig = plot_actuator_bode({G_act{1}, G_act_slow}, {'32 Hz','10 Hz'}, OMEGA(1:2)/(2*pi), 'Name', 'actuator_bode')
%   G_acts   cell of actuator models from define_PT2Actuator; channel (1,1), demanded to actual
%            deflection, is drawn
%   labels   one legend label per model; the title reads 'PT2 actuator: <labels> bandwidth'
%   f_band   [f_lo f_hi] Hz, shaded grey in both tiles: the wind-off frequencies of the two modes
%            that flutter
%   Magnitude in dB on top, unwrapped phase in degrees below, 1 to 100 Hz on a log axis; the
%   first model in blue, the second muted.
%   Options (defaults)
%     'Name'   'actuator Bode'   figure name (the tutorial export uses the file stem)
%   fig      fw_figure of the standard size, 2x1 tiled layout
%   See also define_PT2Actuator, freqresp.

arguments
    G_acts cell
    labels cell
    f_band (1,2) double
    opts.Name (1,:) char = 'actuator Bode'
end

S = fw_style();
f_act = logspace(0,2,300);                           % Hz
n = numel(G_acts);
h = zeros(numel(f_act), n);
for i = 1:n
    h(:,i) = squeeze(freqresp(G_acts{i}(1,1), 2*pi*f_act));
end
col = S.cases(mod(0:n-1, size(S.cases,1)) + 1, :);

fig = fw_figure(S.size.standard(1), S.size.standard(2), 'Name', opts.Name);
tl  = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax  = nexttile(tl); hold(ax,'on')
patch(ax, f_band([1 2 2 1]), [-200 -200 200 200], S.grey, 'EdgeColor','none', ...
    'FaceAlpha',0.5, 'HandleVisibility','off');
for i = 1:n
    plot(ax, f_act, 20*log10(abs(h(:,i))), '-', 'Color',col(i,:), 'LineWidth',S.lineWidth)
end
set(ax,'XScale','log'); xlim(ax,[f_act(1) f_act(end)]); ylim(ax,[-40 10])
ylabel(ax,'magnitude (dB)'); set(ax,'XTickLabel',[])
text(ax, sqrt(f_band(1)*f_band(2)), 6, 'flutter band', 'HorizontalAlignment','center', ...
    'VerticalAlignment','top', 'FontSize',S.fontSizeSmall, 'Color',S.muted)
legend(ax, labels, 'Location','southwest')
title(ax, sprintf('PT2 actuator: %s bandwidth', strjoin(labels, ' and ')))
ax  = nexttile(tl); hold(ax,'on')
patch(ax, f_band([1 2 2 1]), [-400 -400 100 100], S.grey, 'EdgeColor','none', ...
    'FaceAlpha',0.5, 'HandleVisibility','off');
for i = 1:n
    plot(ax, f_act, unwrap(angle(h(:,i)))*180/pi, '-', 'Color',col(i,:), 'LineWidth',S.lineWidth)
end
set(ax,'XScale','log'); xlim(ax,[f_act(1) f_act(end)]); ylim(ax,[-190 10])
xlabel(ax,'$f$ (Hz)'); ylabel(ax,'phase (deg)')
end
