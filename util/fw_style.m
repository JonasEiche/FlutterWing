function S = fw_style()
%FW_STYLE  FlutterWing visual identity: palette, type, line and size tokens (pure data, no graphics).
%   S = fw_style()
%   The apple.com neutrals, one data colour and one rare accent; docs/STYLE.md has the rules.
%     S.ink    #1D1D1F  text, axes, outlines, hinge lines, grid (at alpha), zero and threshold lines
%     S.blue   #0066CC  data: curves, IMU arrows, selected sensors and surfaces, control surfaces
%                       while the controller acts, closed-loop poles, coupling quantities
%     S.coral  #E35336  instability, rarely: the flutter or divergence crossing, a pole in the right
%                       half plane, the velocity readout past the flutter boundary; never a fill
%     S.white  #FFFFFF  figure and animation background (so that a figure blends into the page),
%                       plotting area, legend, wing body
%     S.paper  #F5F5F7  the social preview card only
%     S.grey   #D2D2D7  control surfaces at rest, panel faces
%     S.muted  #6E6E73  secondary text, second series, sample points
%     S.cases  [blue; muted]  axes colour order; fw_figure cycles S.lineStyles after it
%   Rules of thumb: ink is the only text colour, muted for secondary text; coral text only from
%   S.fontSizeDisplay up; a colour appears only with its meaning, and most figures show no coral.
%   Grid: ink at S.gridAlpha (major) and S.minorGridAlpha (minor, solid).
%   Lines in pt: S.lineWidthMesh 0.4 (panel edges), S.lineWidthHair 0.5 (zero and threshold lines),
%     S.lineWidthThin 0.8, S.lineWidth 1.2 (curves), S.lineWidthArrow 1.8, S.lineWidthMarker 2 with
%     S.markerSize 4; S.lineStyles {'-','--',':','-.'}; crossings S.markerFlutter 'o' and
%     S.markerDivergence 's', both hollow coral.
%   Type: every text object goes through the LaTeX interpreter (S.interpreter), i.e. Computer Modern,
%     in three sizes: S.fontSizeSmall 8 (in-figure numbers, annotations, captions), S.fontSize 9
%     (everything else in a figure: ticks, labels, legend, title in normal weight, footnotes, xline
%     labels), S.fontSizeDisplay 14 (animation label and readout). No bold anywhere: hierarchy
%     comes from size and from ink against muted. Math is $V_\infty$, escape % _ &.
%   Sizes in cm for fw_figure: S.size.standard [14 9.8] (V-g, scenes, RFA fit), S.size.wide [14 6]
%     (time histories, planform), S.size.square [9.8 9.8] (pole maps); S.dpi 200 for fw_export PNGs.
%   S.lighting: patch lighting of the foil skin in animate_wing (Ambient, Diffuse, Specular,
%     SpecularExponent).
%   Example:
%     S = fw_style();  fig = fw_figure(S.size.wide(1), S.size.wide(2));
%     plot(t, z, 'Color', S.blue);  xline(t_on, '-', 'Color', S.ink, 'LineWidth', S.lineWidth)
%   See also fw_figure, fw_export, Vg_plot, pole_plot, wing_scene, plot_wing_layout, animate_wing.

S.ink   = [29 29 31]/255;
S.blue  = [0 102 204]/255;
S.coral = [227 83 54]/255;
S.white = [1 1 1];
S.paper = [245 245 247]/255;
S.grey  = [210 210 215]/255;
S.muted = [110 110 115]/255;
S.cases = [S.blue; S.muted];

S.gridAlpha      = 0.15;
S.minorGridAlpha = 0.07;

S.lineStyles      = {'-', '--', ':', '-.'};
S.lineWidthMesh   = 0.4;
S.lineWidthHair   = 0.5;
S.lineWidthThin   = 0.8;
S.lineWidth       = 1.2;
S.lineWidthArrow  = 1.8;
S.lineWidthMarker = 2;
S.markerSize      = 4;
S.markerFlutter   = 'o';
S.markerDivergence = 's';

S.interpreter     = 'latex';
S.fontSizeSmall   = 8;
S.fontSize        = 9;
S.fontSizeDisplay = 14;

S.size.standard = [14, 9.8];
S.size.wide     = [14, 6];
S.size.square   = [9.8, 9.8];
S.dpi           = 200;

S.lighting = struct('Ambient',0.50, 'Diffuse',0.50, 'Specular',0.12, 'SpecularExponent',12);
end
