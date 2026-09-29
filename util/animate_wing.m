function [frames, fig, info] = animate_wing(Structure, q_f, opts)
%ANIMATE_WING  Render a wing motion in 3-D: live figure, GIF or RGB frames, with control surfaces, IMU arrows and the AFS label.
%   [frames, fig, info] = animate_wing(Structure, q_f)
%   [frames, fig, info] = animate_wing(Structure, q_f, Name=Value)
%   One frame per column of q_f: the modal displacement field
%   w(x,y) = build_IMU(x,y,E,ele)*PHIgf*q_f is evaluated at every render vertex,
%   each control surface is rotated rigidly about its own deforming hinge line,
%   and blue arrows show the measured IMU accelerations. The wing is drawn
%   either as the flat DLM panel mesh ('flat') or as a lofted NACA 4-digit skin
%   ('foil'). All kinematics run in structural coordinates, where z points DOWN;
%   the render flips z, so physical up is screen up and a positive control
%   surface angle shows as an edge hanging down.
%
%   Structure  model dictionary from define_Goland_Structure_Aero(num_modes,
%              num_poles) (1 trailing-edge flap on the outer 3 of 10 spanwise
%              panels, 1 IMU) or define_RectWing_Structure_Aero(...) (4 TE flaps
%              + 4 LE slats, 8 IMUs). Fields used: Ps (panel corners in
%              structural coordinates, {panel}{corner} = [x;y;z], corner order
%              LE-inboard, TE-inboard, TE-outboard, LE-outboard), cspanels
%              (panel index list per control surface), E and ele (FEM beam
%              elements, for build_IMU) and PHIgf. Rectangular unswept
%              planforms only (asserted).
%   q_f        [num_modes x Nt] real modal displacements (m), one column per
%              frame, or a complex [num_modes x 1] eigenvector, which switches
%              to the eigenmode loop (see Cycles below).
%
%   Motion options (default in brackets)
%     PHIgf      [Structure.PHIgf]  mode shapes [n_dof x num_modes]
%     Surfaces   [zeros(n_cs,Nt)]   control-surface angles (rad), one row per
%                                   entry of Structure.cspanels; model sign
%                                   convention: positive = edge down, for
%                                   trailing-edge and leading-edge surfaces
%     IMUPos     []                 [n_imu x 2] IMU (x,y) in structural
%                                   coordinates; default: the midpoint of every
%                                   control-surface hinge axis, i.e. the IMU
%                                   positions of the define_* models
%     IMUAccel   []                 [n_imu x Nt] measured accelerations, any
%                                   unit; arrow length =
%                                   0.10*span*min(|a|/ArrowRef, 1.3)
%     ArrowRef   []                 default: max|IMUAccel| over the frames
%                                   before the first Active frame (over all
%                                   frames if none is active)
%     Arrows     []                 default: true when IMUAccel is given
%     Active     [false(1,Nt)]      controller ON per frame: sets the label text
%                                   and colour, and paints ALL control-surface
%                                   bodies Colors.blue (Colors.grey when off)
%     LabelOff   ['active flutter suppression (AFS) controller OFF']
%     LabelOn    ['active flutter suppression (AFS) controller ON']
%                                   top-left textbox [0.015 0.85 0.70 0.12],
%                                   Colors.fontSizeDisplay, LaTeX, regular weight,
%                                   Colors.muted when off and Colors.ink when
%                                   on; drawn only when Active or Labels was given
%     Labels     [{}]               1 x Nt cellstr, per-frame label override
%     Readout    [{}]               1 x Nt cellstr, top-right textbox
%                                   [0.60 0.85 0.385 0.12], right aligned,
%                                   Colors.fontSizeDisplay, LaTeX: write math as
%                                   '$V_\infty = 108.0$ m/s' and escape % _ &
%     ReadoutColor [Colors.ink]      colour of the readout text
%     ReadoutHighlight [false(1,Nt)] per frame: Colors.coral instead of ReadoutColor
%     Footnote   ['']               bottom-right textbox [0.60 0.02 0.385 0.06],
%                                   Colors.fontSize, Colors.muted, LaTeX
%   Appearance
%     Skin       ['flat']           'flat' (DLM panel mesh) or 'foil' (lofted skin)
%     NACA       ['2408']           4-digit code MPTT of the foil skin: camber
%                                   M/100 at P/10 chord, thickness TT/100.
%                                   '0012' is symmetric (zero camber line)
%     FlapChord  []                 drawn chord fraction of the TE surfaces
%                                   (foil: where the cut is; flat: ignored, the
%                                   strip is the panel block); [] = the model
%                                   hinge chord fraction from the panel corners
%     SlatChord  []                 the same for the LE surfaces
%     ShowMesh   [true]             flat skin: draw the panel edges (Colors.edge,
%                                   EdgeAlpha 0.32, LineWidth 0.4)
%     TipAmplitude [0.07]           scale the peak |deflection| over all frames
%                                   and render vertices to this fraction of the
%                                   semispan (one global scale, info.tipScale);
%                                   [] = physical (scale 1)
%     SurfaceAmplitude [deg2rad(30)] scalar: one common scale that brings the
%                                   peak |angle| over all surfaces and frames to
%                                   this value; [TE LE]: one scale per edge
%                                   group; [] = physical. Never per surface;
%                                   see info.surfaceScale
%     View       [[37 26]]          view(az,el)
%     Light      ['right']          camlight argument
%     Colors     [fw_style()]       colour and lighting struct
%     Width      [900]              frame width (px) before cropping
%     Aspect     [410/900]          frame height/width
%     CropPad    [10]               px of background kept around the content
%   Output
%     Output     ['figure']         'figure', 'gif', 'frames', or a cellstr of
%                                   several. 'figure': visible figure, the frames
%                                   are played Repeat times, one every FrameDelay
%                                   s, the figure stays open with a Replay
%                                   button and is returned,
%                                   frames = {}. 'gif': frames captured, cropped
%                                   and written to GifFile, figure closed,
%                                   fig = []. 'frames': the cropped frames are
%                                   returned as {Nt x 1} of uint8 [h x w x 3].
%                                   'figure' plus 'gif'/'frames' captures from
%                                   the visible figure
%     GifFile    ['']               target file, required for Output 'gif'
%     FrameDelay [0.05]             (s) frame period of the GIF and of the
%                                   figure playback
%     Loop       [Inf]              GIF LoopCount
%     PaletteFrames []              frames sampled for the ONE 256-colour GIF
%                                   palette (default: [1, round(Nt/3),
%                                   round(2*Nt/3), Nt] plus the first Active
%                                   frame), so the GIF has no palette flicker
%     Repeat     [1]                figure mode only: how often the loop is played
%     Verbose    [true]             print the one-line summary
%   Eigenmode loop (complex q_f [num_modes x 1])
%     Cycles     [[3 1]]            growth cycles, fade cycles
%     FramesPerCycle [14]           Nt = sum(Cycles)*FramesPerCycle
%     StartAmplitude [0.25]         envelope value at the first and last frame.
%                                   q_f is phase-aligned so that the render
%                                   vertex of largest |deflection| is at phase 0
%                                   and deflected UP (screen up) there;
%                                   Surfaces, IMUAccel and Active default to none
%
%   OUTPUT
%     frames  {Nt x 1} cropped uint8 [h x w x 3] frames; {} unless Output
%             contains 'gif' or 'frames'
%     fig     figure handle in figure mode, [] otherwise
%     info    .Nt, .cropbox [r1 r2 c1 c2] ([] without a capture), .size_px [h w]
%             of the cropped frames (the uncropped render size without a
%             capture), .gifBytes (NaN if no GIF was written), .tipScale,
%             .surfaceScale [1 x n_cs], .geom (xLE, xTE, c, span, hinges: struct
%             array with .xh .y1 .y2 .isTE, imu [n_imu x 2]), .secondsPerFrame,
%             .edgeZ [n_cs x Nt] display z of each control surface's moving edge,
%             averaged over its edge vertices (negative = edge below the wing
%             plane, i.e. down: the sign check of the drawing)
%
%   Example - open loop to closed loop from a simulate_afs_switch run:
%     [Structure, ~] = define_RectWing_Structure_Aero(5, 6);
%     G   = build_G_RectWing(108, 1:8, 1:8);
%     load RectWing_Cont_imu18_ail18.mat              % RectWing_Cont: the shipped static gain (R01)
%     sim = simulate_afs_switch(G, RectWing_Cont);
%     animate_wing(Structure, sim.q_f, Surfaces=sim.delta, IMUAccel=sim.acc, ...
%         Active=sim.active, Skin='foil', NACA='2408', FlapChord=0.25, ...
%         SlatChord=0.20, TipAmplitude=0.07, SurfaceAmplitude=deg2rad([33 25]), ...
%         Output='gif', GifFile='my_afs_demo.gif');
%
%   Example - flutter eigenmode loop (56 frames) in a live figure:
%     [Structure, Aero] = define_RectWing_Structure_Aero(5, 6);
%     A = build_ABCD_G(108, Structure, Aero);
%     [v, e] = eig(A);  e = diag(e);
%     cand = find(real(e) > 0 & imag(e) > 3);
%     [~, imax] = max(real(e(cand)));
%     animate_wing(Structure, v(1:5, cand(imax)), TipAmplitude=0.13, Repeat=3);
%
%   Runtime (R2026a): foil skin about 0.25 s per frame (168 frames in 43 s), flat
%   skin 0.12 s per frame; figure mode adds the screen drawing.
%
%   See also fw_style, build_IMU, define_RectWing_Structure_Aero,
%   define_Goland_Structure_Aero, simulate_afs_switch.

arguments
    Structure (1,1) struct
    q_f {mustBeNumeric, mustBeNonempty}
    opts.PHIgf {mustBeNumeric} = []
    opts.Surfaces {mustBeNumeric} = []
    opts.IMUPos {mustBeNumeric} = []
    opts.IMUAccel {mustBeNumeric} = []
    opts.ArrowRef {mustBeNumeric} = []
    opts.Arrows {mustBeNumericOrLogical} = []
    opts.Active {mustBeNumericOrLogical} = []
    opts.LabelOff {mustBeTextScalar} = 'active flutter suppression (AFS) controller OFF'
    opts.LabelOn {mustBeTextScalar} = 'active flutter suppression (AFS) controller ON'
    opts.Labels = {}
    opts.Readout = {}
    opts.ReadoutColor {mustBeNumeric} = []
    opts.ReadoutHighlight {mustBeNumericOrLogical} = []
    opts.Footnote {mustBeTextScalar} = ''
    opts.Skin {mustBeMember(opts.Skin, {'flat','foil'})} = 'flat'
    opts.NACA {mustBeTextScalar} = '2408'
    opts.FlapChord {mustBeNumeric} = []
    opts.SlatChord {mustBeNumeric} = []
    opts.ShowMesh (1,1) logical = true
    opts.TipAmplitude {mustBeNumeric} = 0.07
    opts.SurfaceAmplitude {mustBeNumeric} = deg2rad(30)
    opts.View (1,2) double = [37 26]
    opts.Light {mustBeTextScalar} = 'right'
    opts.Colors (1,1) struct = fw_style()
    opts.Width (1,1) double {mustBePositive} = 900
    opts.Aspect (1,1) double {mustBePositive} = 410/900
    opts.CropPad (1,1) double {mustBeNonnegative} = 10
    opts.Output = 'figure'
    opts.GifFile {mustBeTextScalar} = ''
    opts.FrameDelay (1,1) double {mustBePositive} = 0.05
    opts.Loop (1,1) double {mustBeNonnegative} = Inf
    opts.PaletteFrames {mustBeNumeric} = []
    opts.Repeat (1,1) double {mustBePositive} = 1
    opts.Verbose (1,1) logical = true
    opts.Cycles (1,2) double {mustBeNonnegative} = [3 1]
    opts.FramesPerCycle (1,1) double {mustBePositive} = 14
    opts.StartAmplitude (1,1) double {mustBePositive} = 0.25
end

t_start = tic;
S = opts.Colors;

% ---- 1. requested outputs ------------------------------------------------
outs = opts.Output;
if ischar(outs) || isstring(outs), outs = cellstr(outs); end
if ~iscellstr(outs) || isempty(outs)
    error('animate_wing:Output', ...
        'Output must be a char or a cellstr of ''figure'', ''gif'', ''frames''.');
end
for i = 1:numel(outs)
    if ~any(strcmp(outs{i}, {'figure','gif','frames'}))
        error('animate_wing:Output', 'unknown Output ''%s''.', outs{i});
    end
end
wantFig = any(strcmp('figure', outs));
wantGif = any(strcmp('gif',    outs));
wantFrm = any(strcmp('frames', outs));
gifFile = char(opts.GifFile);
if wantGif && isempty(gifFile)
    error('animate_wing:GifFile', 'Output ''gif'' requires GifFile.');
end
doCapture = wantGif || wantFrm;

% ---- 2. planform geometry and control surfaces ---------------------------
geom  = wing_geometry(Structure);
span  = geom.span;
xLE   = geom.xLE;
xTE   = geom.xTE;
cch   = geom.c;
n_cs  = numel(geom.hinges);
isTEs = false(1, n_cs);
for j = 1:n_cs, isTEs(j) = geom.hinges(j).isTE; end

if isempty(opts.PHIgf)
    if ~isfield(Structure, 'PHIgf')
        error('animate_wing:PHIgf', 'Structure has no PHIgf field; pass PHIgf=...');
    end
    PHIgf = Structure.PHIgf;
else
    PHIgf = opts.PHIgf;
end
num_modes = size(PHIgf, 2);
if size(q_f,1) ~= num_modes
    error('animate_wing:qf', 'q_f has %d rows, PHIgf has %d columns.', ...
        size(q_f,1), num_modes);
end

% ---- 3. render bodies ----------------------------------------------------
isFoil = strcmp(char(opts.Skin), 'foil');
g_x = 0.005;                                  % half chord gap at the cuts
g_y = 0.0015*span;                            % spanwise gap between bodies

flapChord = opts.FlapChord;
if isempty(flapChord) && any(isTEs)
    xh = zeros(1,0);
    for j = find(isTEs), xh(end+1) = geom.hinges(j).xh; end %#ok<AGROW>
    flapChord = max((xh - xTE)/cch);
end
slatChord = opts.SlatChord;
if isempty(slatChord) && any(~isTEs)
    xh = zeros(1,0);
    for j = find(~isTEs), xh(end+1) = geom.hinges(j).xh; end %#ok<AGROW>
    slatChord = max((xLE - xh)/cch);
end

[zupf, zlof] = naca_funs(char(opts.NACA));
if isFoil
    [B, hingeY, hingeX] = build_bodies_foil(geom, isTEs, flapChord, slatChord, ...
        zupf, zlof, g_x, g_y);
    zoff = 0.004*span;
else
    [B, hingeY, hingeX] = build_bodies_flat(geom);
    zoff = 0.003*span;
end
nB = numel(B);

% global render-vertex arrays (ring vertices only; cap centroids are per body)
nvert = 0;
for b = 1:nB
    B(b).idx = nvert + (1:numel(B(b).x));
    nvert = nvert + numel(B(b).x);
end
xAll = zeros(nvert,1);  yAll = zeros(nvert,1);  a0All = zeros(nvert,1);
for b = 1:nB
    xAll(B(b).idx)  = B(b).x;
    yAll(B(b).idx)  = B(b).y;
    a0All(B(b).idx) = B(b).a0;                % upward thickness offset
end
bodyOfCs = zeros(1, n_cs);
for b = 1:nB
    if B(b).cs > 0, bodyOfCs(B(b).cs) = b; end
end

% moving-edge vertices per control surface (the sign check of the drawing)
edgeIdx = cell(1, n_cs);
for j = 1:n_cs
    b = bodyOfCs(j);
    if b == 0, continue, end
    xb = B(b).x;
    if isTEs(j), xe = min(xb); else, xe = max(xb); end
    edgeIdx{j} = B(b).idx(abs(xb - xe) < 1e-9*max(1,cch) + 1e-12);
end

% ---- 4. deflection maps --------------------------------------------------
Mmesh = deflection_map(xAll, yAll, Structure, PHIgf, span);
Mh = cell(1, n_cs);
for j = 1:n_cs
    Mh{j} = deflection_map([hingeX(j); hingeX(j)], hingeY(j,:)', ...
        Structure, PHIgf, span);
end

if isempty(opts.IMUPos)
    imuXY = geom.imu;
else
    imuXY = opts.IMUPos;
    if size(imuXY,2) ~= 2
        error('animate_wing:IMUPos', 'IMUPos must be [n_imu x 2].');
    end
end
n_imu = size(imuXY,1);
Mimu   = zeros(n_imu, num_modes);
a_imu  = zeros(n_imu, 1);
csOfImu = zeros(n_imu, 1);
if n_imu > 0
    Mimu = deflection_map(imuXY(:,1), imuXY(:,2), Structure, PHIgf, span);
    if isFoil
        xi = min(max((xLE - imuXY(:,1))/cch, 0), 1);
        a_imu = zupf(xi)*cch;                 % arrows sit on the upper surface
    end
    for i = 1:n_imu                           % IMU rides with its own surface
        for j = 1:n_cs
            b = bodyOfCs(j);
            if b == 0, continue, end
            tolg = g_y + 1e-9;
            if imuXY(i,2) >= min(B(b).y) - tolg && imuXY(i,2) <= max(B(b).y) + tolg ...
                    && imuXY(i,1) >= min(B(b).x) - 1e-9 && imuXY(i,1) <= max(B(b).x) + 1e-9
                csOfImu(i) = j;
                break
            end
        end
    end
end

% ---- 5. frame series -----------------------------------------------------
if size(q_f,2) == 1 && ~isreal(q_f)
    fpc = opts.FramesPerCycle;
    ng  = round(opts.Cycles(1)*fpc);
    nfa = round(opts.Cycles(2)*fpc);
    Nt  = ng + nfa;
    if Nt < 1
        error('animate_wing:Cycles', 'Cycles*FramesPerCycle must give at least one frame.');
    end
    A0  = opts.StartAmplitude;
    zc  = Mmesh*q_f;                          % phase align on the largest vertex
    [~, imx] = max(abs(zc));
    q0  = -q_f * exp(-1i*angle(zc(imx)));      % minus: peak deflection physically UP (z_s < 0), screen up
    env = zeros(1, Nt);
    if ng > 0,  env(1:ng)     = A0.^(1 - (0:ng-1)/ng);                     end
    if nfa > 0, env(ng+1:Nt)  = A0 + (1-A0)*0.5*(1 + cos(pi*(1:nfa)/nfa)); end
    q_all = real(q0 * (env .* exp(1i*2*pi*(0:Nt-1)/fpc)));
elseif ~isreal(q_f)
    error('animate_wing:qf', ...
        'a complex q_f must be a single column (the eigenmode loop).');
else
    q_all = q_f;
    Nt    = size(q_all,2);
end

% one global deflection scale
W0 = Mmesh*q_all;
pk = max(abs(W0(:)));
if isempty(opts.TipAmplitude) || pk == 0
    tipScale = 1;
else
    tipScale = opts.TipAmplitude*span/pk;
end
q_all = tipScale*q_all;
Wall  = tipScale*W0;
Amax  = max(max(abs(Wall(:))), 0.05*span);
Wimu  = Mimu*q_all;

% control-surface angles and their scaling
Sur = opts.Surfaces;
if isempty(Sur)
    Sur = zeros(n_cs, Nt);
else
    if isvector(Sur) && n_cs == 1, Sur = reshape(Sur, 1, []); end
    if size(Sur,1) ~= n_cs || size(Sur,2) ~= Nt
        error('animate_wing:Surfaces', 'Surfaces must be [%d x %d].', n_cs, Nt);
    end
end
surfaceScale = ones(1, n_cs);
sa = opts.SurfaceAmplitude;
if ~isempty(sa) && n_cs > 0
    if isscalar(sa)
        pks = max(abs(Sur(:)));
        if pks > 0, surfaceScale(:) = sa/pks; end
    elseif numel(sa) == 2
        grp = {isTEs, ~isTEs};
        for g = 1:2
            rows = grp{g};
            if ~any(rows), continue, end
            pks = max(max(abs(Sur(rows,:))));
            if ~isempty(pks) && pks > 0, surfaceScale(rows) = sa(g)/pks; end
        end
    else
        error('animate_wing:SurfaceAmplitude', ...
            'SurfaceAmplitude must be a scalar, a [TE LE] pair or [].');
    end
end
delta = surfaceScale(:).*Sur;
if isempty(delta), delta = zeros(0, Nt); end

% ---- 6. label, readout, arrows ------------------------------------------
activeGiven = ~isempty(opts.Active);
if activeGiven
    Active = logical(reshape(opts.Active, 1, []));
    if numel(Active) ~= Nt
        error('animate_wing:Active', 'Active must have %d entries.', Nt);
    end
else
    Active = false(1, Nt);
end
iOn = find(Active, 1);

labels = opts.Labels;
if ischar(labels) || isstring(labels), labels = cellstr(labels); end
if ~isempty(labels) && numel(labels) ~= Nt
    error('animate_wing:Labels', 'Labels must be a 1 x %d cellstr.', Nt);
end
labelOn = activeGiven || ~isempty(labels);

readout = opts.Readout;
if ischar(readout) || isstring(readout), readout = cellstr(readout); end
if ~isempty(readout) && numel(readout) ~= Nt
    error('animate_wing:Readout', 'Readout must be a 1 x %d cellstr.', Nt);
end
rdCol = opts.ReadoutColor;
if isempty(rdCol), rdCol = S.ink; end
if isempty(opts.ReadoutHighlight)
    rdHi = false(1, Nt);
else
    rdHi = logical(reshape(opts.ReadoutHighlight, 1, []));
    if numel(rdHi) ~= Nt
        error('animate_wing:ReadoutHighlight', 'ReadoutHighlight must have %d entries.', Nt);
    end
end

acc = opts.IMUAccel;
if isempty(opts.Arrows)
    drawArrows = ~isempty(acc);
else
    drawArrows = logical(opts.Arrows);
end
Lr = zeros(n_imu, Nt);
if drawArrows
    if isempty(acc)
        error('animate_wing:Arrows', 'Arrows=true requires IMUAccel.');
    end
    if n_imu == 0
        error('animate_wing:Arrows', 'no IMU positions to put arrows on.');
    end
    if isvector(acc) && n_imu == 1, acc = reshape(acc, 1, []); end
    if size(acc,1) ~= n_imu || size(acc,2) ~= Nt
        error('animate_wing:IMUAccel', 'IMUAccel must be [%d x %d].', n_imu, Nt);
    end
    aref = opts.ArrowRef;
    if isempty(aref)
        if isempty(iOn) || iOn == 1
            aref = max(abs(acc(:)));
        else
            aref = max(max(abs(acc(:,1:iOn-1))));
        end
    end
    if isempty(aref) || ~(aref > 0), aref = 1; end
    Lr = min(abs(acc)/aref, 1.3);
end

% ---- 7. figure -----------------------------------------------------------
W_px = round(opts.Width);
H_px = round(opts.Width*opts.Aspect);
if mod(H_px,2) ~= 0, H_px = H_px + 1; end
vis = {};
if ~wantFig, vis = {'Visible','off'}; end   % an explicit 'on' pops a live script's figure out of the Live Editor
fig = figure(vis{:}, 'Units','pixels', 'Position',[50 50 W_px H_px], ...
    'Color',S.white, 'InvertHardcopy','off'); %#ok<INVHCRM> keeps Color on print
                                             % up to R2025a; no effect after
if isprop(fig, 'Theme'), fig.Theme = 'light'; end   % never inherit a dark desktop theme
set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 W_px H_px]/96, ...
    'PaperPositionMode','manual');      % print -r192 gives exactly 2*Width px
ax = axes(fig, 'Position',[-0.06 -0.16 1.12 1.32], 'Clipping','off');
hold(ax, 'on');

hPat = gobjects(1, nB);
hLn  = cell(1, nB);
for b = 1:nB
    if B(b).cs > 0, col = S.grey; else, col = S.white; end
    if isFoil || ~opts.ShowMesh
        eArgs = {'EdgeColor','none'};
    else
        eArgs = {'EdgeColor',S.ink, 'EdgeAlpha',0.35, 'LineWidth',S.lineWidthMesh};
    end
    hPat(b) = patch(ax, 'Vertices',zeros(numel(B(b).x)+2*B(b).hasCaps,3), ...
        'Faces',B(b).faces, 'FaceColor',col, eArgs{:}, ...
        'FaceLighting','gouraud', 'BackFaceLighting','reverselit', ...
        'AmbientStrength',S.lighting.Ambient, 'DiffuseStrength',S.lighting.Diffuse, ...
        'SpecularStrength',S.lighting.Specular, ...
        'SpecularExponent',S.lighting.SpecularExponent);
    hLn{b} = gobjects(1, numel(B(b).lines));
    for i = 1:numel(B(b).lines)
        hLn{b}(i) = line(ax, nan, nan, nan, 'Color',S.ink, 'LineWidth',B(b).lw);
    end
end

arrL = 0.10*span;
wh   = 0.11*arrL*[cosd(37), sind(37)];
hh   = 0.30*arrL;
hShaft = gobjects(0);  hHead = gobjects(0);
if drawArrows
    hShaft = line(ax, nan, nan, nan, 'Color',S.blue, 'LineWidth',S.lineWidthArrow);
    hHead  = patch(ax, 'Vertices',zeros(3*n_imu,3), ...
        'Faces',reshape(1:3*n_imu, 3, n_imu)', 'FaceColor',S.blue, ...
        'EdgeColor','none', 'FaceLighting','none');
end

hLbl = gobjects(0);  hRd = gobjects(0);
if labelOn
    hLbl = annotation(fig, 'textbox', [0.015 0.85 0.70 0.12], 'String','', ...
        'FontSize',S.fontSizeDisplay, 'EdgeColor','none', ...
        'Color',S.muted, 'FitBoxToText','off', 'VerticalAlignment','top', ...
        'Interpreter',S.interpreter);
end
if ~isempty(readout)
    hRd = annotation(fig, 'textbox', [0.60 0.85 0.385 0.12], 'String','', ...
        'FontSize',S.fontSizeDisplay, 'EdgeColor','none', 'Color',rdCol, ...
        'FitBoxToText','off', 'HorizontalAlignment','right', ...
        'VerticalAlignment','top', 'Interpreter',S.interpreter);
end
if ~isempty(char(opts.Footnote))
    annotation(fig, 'textbox', [0.60 0.02 0.385 0.06], ...
        'String',char(opts.Footnote), 'FontSize',S.fontSize, 'EdgeColor','none', ...
        'Color',S.muted, 'FitBoxToText','off', 'HorizontalAlignment','right', ...
        'VerticalAlignment','bottom', 'Interpreter',S.interpreter);
end

axis(ax, 'off');  daspect(ax, [1 1 1]);
xlim(ax, [xTE-0.1, xLE+0.1]);
ylim(ax, [-0.1, span+0.15]);
zlim(ax, [-1.9, 1.9]*Amax);
ax.Projection = 'orthographic';
view(ax, opts.View(1), opts.View(2));
camlight(ax, char(opts.Light));

% ---- 8. frames -----------------------------------------------------------
edgeZ   = nan(n_cs, Nt);
prevAct = [];
playing = false;
if wantFig
    play_figure(round(opts.Repeat));
end
frames = {};
info.cropbox = [];
info.size_px = [H_px W_px];
info.gifBytes = NaN;
if doCapture
    tmppng = [tempname '.png'];
    frames = cell(Nt,1);
    for k = 1:Nt
        set_frame(k);
        frames{k} = grab_frame(fig, tmppng);
    end
    if exist(tmppng, 'file'), delete(tmppng); end
    if abs(size(frames{1},2) - W_px) > 2
        warning('animate_wing:frameSize', ...
            'captured frame is %d px wide, expected %d px.', size(frames{1},2), W_px);
    end
    box = union_crop(frames, S.white, round(opts.CropPad));
    for k = 1:Nt
        frames{k} = frames{k}(box(1):box(2), box(3):box(4), :);
    end
    info.cropbox = box;
    info.size_px = [size(frames{1},1), size(frames{1},2)];
end
if wantFig
    % added after the capture so the button never appears in a GIF frame
    uicontrol(fig, 'Style','pushbutton', 'String','Replay', 'Units','pixels', ...
        'Position',[10 10 70 24], 'Callback',@(~,~) play_figure(1));
else
    close(fig);
    fig = [];
end

% ---- 9. GIF --------------------------------------------------------------
if wantGif
    pf = opts.PaletteFrames;
    if isempty(pf)
        pf = [1, round(Nt/3), round(2*Nt/3), Nt, iOn];
    end
    pf = unique(min(max(round(pf(:)'), 1), Nt));
    info.gifBytes = write_gif(frames, gifFile, opts.FrameDelay, opts.Loop, pf);
end
if ~wantFrm
    frames = {};
end

% ---- 10. info and report ------------------------------------------------
t_el = toc(t_start);
info.Nt = Nt;
info.tipScale = tipScale;
info.surfaceScale = surfaceScale;
info.geom = struct('xLE',xLE, 'xTE',xTE, 'c',cch, 'span',span, ...
    'hinges',geom.hinges, 'imu',imuXY);
info.secondsPerFrame = t_el/max(Nt,1);
info.edgeZ = edgeZ;
if opts.Verbose
    fprintf(['animate_wing: %d frames, %d x %d px, tip scale %.3g, ' ...
        'surface scale [%s], %.1f s (%.2f s/frame)\n'], Nt, info.size_px(2), ...
        info.size_px(1), tipScale, strtrim(sprintf('%.3g ', surfaceScale)), ...
        t_el, info.secondsPerFrame);
    if wantGif
        fprintf('%s written, %.0f kB\n', gifFile, info.gifBytes/1024);
    end
end

% =========================================================================
    function set_frame(k)
        % structural coordinates: z points down, thickness offsets are up
        Vs = [xAll, yAll, Wall(:,k) - a0All];
        if n_imu > 0
            pinS = [imuXY, Wimu(:,k) - a_imu];
        else
            pinS = zeros(0,3);
        end
        for jj = 1:n_cs
            bj = bodyOfCs(jj);
            ii = find(csOfImu == jj);
            if bj == 0 && isempty(ii), continue, end
            zh  = Mh{jj}*q_all(:,k);
            P0  = [hingeX(jj), hingeY(jj,1), zh(1)];
            aax = [0, hingeY(jj,2)-hingeY(jj,1), zh(2)-zh(1)];
            aax = aax/norm(aax);
            if bj > 0
                Vs(B(bj).idx,:) = rotate_about_hinge(Vs(B(bj).idx,:), P0, aax, delta(jj,k));
            end
            if ~isempty(ii)
                pinS(ii,:) = rotate_about_hinge(pinS(ii,:), P0, aax, delta(jj,k));
            end
        end
        Vd = [Vs(:,1), Vs(:,2), -Vs(:,3)];          % display: z up

        setCol = isempty(prevAct) || prevAct ~= Active(k);
        for bb = 1:nB
            vb = Vd(B(bb).idx,:);
            if B(bb).hasCaps
                na = B(bb).Na;
                vb = [vb; mean(vb(1:na,:),1); mean(vb(end-na+1:end,:),1)]; %#ok<AGROW>
            end
            set(hPat(bb), 'Vertices', vb);
            if setCol && B(bb).cs > 0
                if Active(k), set(hPat(bb), 'FaceColor', S.blue);
                else,         set(hPat(bb), 'FaceColor', S.grey);
                end
            end
            for i2 = 1:numel(B(bb).lines)
                li = B(bb).lines{i2};
                set(hLn{bb}(i2), 'XData',vb(li,1), 'YData',vb(li,2), ...
                    'ZData',vb(li,3)+zoff);
            end
        end
        prevAct = Active(k);
        for jj = 1:n_cs
            if ~isempty(edgeIdx{jj}), edgeZ(jj,k) = mean(Vd(edgeIdx{jj},3)); end
        end

        if drawArrows
            pz = -pinS(:,3) + 2*zoff;               % base on the upper surface
            [sx, sy, sz, hv] = arrow_geometry(pinS(:,1), pinS(:,2), pz, ...
                arrL*Lr(:,k), wh, hh);
            set(hShaft, 'XData',sx, 'YData',sy, 'ZData',sz);
            set(hHead, 'Vertices', hv);
        end
        if labelOn
            if ~isempty(labels), str = labels{k};
            elseif Active(k),    str = char(opts.LabelOn);
            else,                str = char(opts.LabelOff);
            end
            if Active(k), col = S.ink; else, col = S.muted; end
            set(hLbl, 'String',str, 'Color',col);
        end
        if ~isempty(readout)
            if rdHi(k), col = S.coral; else, col = rdCol; end
            set(hRd, 'String',readout{k}, 'Color',col);
        end
    end

% =========================================================================
    function play_figure(nRep)
        % one frame every FrameDelay seconds, the timing of the GIF
        if playing, return; end
        playing = true;
        for rr = 1:nRep
            t0 = tic;
            for kk = 1:Nt
                if ~isvalid(fig), playing = false; return; end
                set_frame(kk);
                drawnow;
                pause(max(0, kk*opts.FrameDelay - toc(t0)));
            end
        end
        playing = false;
    end
end

% =========================================================================
function geom = wing_geometry(Structure)
%WING_GEOMETRY  Planform, control-surface hinge axes and default IMU positions.
Ps = Structure.Ps;
npan = numel(Ps);
xcn = zeros(4*npan,1);  ycn = zeros(4*npan,1);
for i = 1:npan
    for k = 1:4
        xcn(4*(i-1)+k) = Ps{i}{k}(1);
        ycn(4*(i-1)+k) = Ps{i}{k}(2);
    end
end
geom.xcn = xcn;  geom.ycn = ycn;  geom.npan = npan;
geom.span = max(ycn);
geom.xLE  = max(xcn);                  % structural x points to the leading edge
geom.xTE  = min(xcn);
geom.c    = geom.xLE - geom.xTE;
tol = 1e-6*max(1, geom.c);
ys = unique(round(ycn, 9));
for i = 1:numel(ys)
    sel = abs(ycn - ys(i)) < 1e-9*max(1, geom.span) + 1e-12;
    if abs(max(xcn(sel)) - geom.xLE) > tol || abs(min(xcn(sel)) - geom.xTE) > tol
        error('animate_wing:planform', ...
            ['animate_wing draws rectangular unswept wings only: at y = %g the ' ...
             'chord runs from x = %g to %g, not from %g to %g.'], ys(i), ...
            min(xcn(sel)), max(xcn(sel)), geom.xTE, geom.xLE);
    end
end

cs = Structure.cspanels;
n_cs = numel(cs);
geom.cspanels = cs;
geom.hinges = struct('xh',{}, 'y1',{}, 'y2',{}, 'isTE',{});
geom.imu = zeros(n_cs, 2);
for j = 1:n_cs
    blk = cs{j};
    ids = reshape(4*(blk(:)'-1) + (1:4)', [], 1);
    isTE = abs(min(xcn(ids)) - geom.xTE) < tol;
    isLE = abs(max(xcn(ids)) - geom.xLE) < tol;
    if isTE
        P1 = Ps{blk(1)}{1};      P2 = Ps{blk(end)}{4};     % axis root -> tip
    elseif isLE
        P1 = Ps{blk(end)}{3};    P2 = Ps{blk(1)}{2};       % axis tip -> root
    else
        error('animate_wing:surface', ...
            ['control surface %d touches neither the leading nor the trailing ' ...
             'edge; animate_wing draws edge surfaces only.'], j);
    end
    geom.hinges(j) = struct('xh',P1(1), 'y1',P1(2), 'y2',P2(2), 'isTE',isTE);
    HP = 0.5*(P1 + P2);
    geom.imu(j,:) = [HP(1), HP(2)];
end
end

% =========================================================================
function M = deflection_map(x, y, Structure, PHIgf, span)
%DEFLECTION_MAP  Modal displacement (positive down) at arbitrary (x,y) points.
y = min(max(y(:), 0), span - 1e-9);         % build_IMU needs points in an element
M = build_IMU(x(:), y, Structure.E, Structure.ele) * PHIgf;
end

% =========================================================================
function V = rotate_about_hinge(V, P0, aax, ang)
%ROTATE_ABOUT_HINGE  Rodrigues rotation of [n x 3] points about the hinge line.
if ang == 0, return, end
W = V - P0;
V = P0 + cos(ang)*W + sin(ang)*cross(repmat(aax, size(W,1), 1), W, 2) ...
    + (1-cos(ang))*(W*aax')*aax;
end

% =========================================================================
function [sx, sy, sz, hv] = arrow_geometry(px, py, pz, Lk, wh, hh)
%ARROW_GEOMETRY  NaN-separated arrow shafts and triangular heads (display coords).
n  = numel(px);
hk = min(hh, 0.6*Lk);
zt = pz + Lk;
sx = reshape([px'; px';           nan(1,n)], [], 1);
sy = reshape([py'; py';           nan(1,n)], [], 1);
sz = reshape([pz'; (zt-0.8*hk)';  nan(1,n)], [], 1);
hv = zeros(3*n, 3);
for j = 1:n
    hv(3*j-2,:) = [px(j),       py(j),       zt(j)];
    hv(3*j-1,:) = [px(j)+wh(1), py(j)+wh(2), zt(j)-hk(j)];
    hv(3*j  ,:) = [px(j)-wh(1), py(j)-wh(2), zt(j)-hk(j)];
end
end

% =========================================================================
function [zupf, zlof] = naca_funs(code)
%NACA_FUNS  Upper and lower surface of a NACA 4-digit foil (chord fractions).
if numel(code) ~= 4 || any(code < '0' | code > '9')
    error('animate_wing:NACA', 'NACA must be a 4-digit code, e.g. ''2408''.');
end
mC = str2double(code(1))/100;
pC = str2double(code(2))/10;
tk = str2double(code(3:4))/100;
ytf = @(xi) 5*tk*(0.2969*sqrt(xi) - 0.1260*xi - 0.3516*xi.^2 ...
                + 0.2843*xi.^3 - 0.1036*xi.^4);            % thickness, closed TE
if mC == 0
    ycf = @(xi) zeros(size(xi));                           % symmetric section
elseif pC <= 0 || pC >= 1
    error('animate_wing:NACA', ...
        'NACA %s: a cambered section needs the camber position 1..9.', code);
else
    ycf = @(xi) (mC/pC^2)*(2*pC*xi - xi.^2).*(xi <= pC) ...
              + (mC/(1-pC)^2)*((1-2*pC) + 2*pC*xi - xi.^2).*(xi > pC);
end
zupf = @(xi) ycf(xi) + ytf(xi);
zlof = @(xi) ycf(xi) - ytf(xi);
end

% =========================================================================
function [B, hingeY, hingeX] = build_bodies_flat(geom)
%BUILD_BODIES_FLAT  One patch of panel quads per body: main plus one per surface.
xcn = geom.xcn;  ycn = geom.ycn;
n_cs = numel(geom.hinges);
hingeY = zeros(n_cs, 2);
hingeX = zeros(1, n_cs);
for j = 1:n_cs
    hingeY(j,:) = [geom.hinges(j).y1, geom.hinges(j).y2];
    hingeX(j)   = geom.hinges(j).xh;       % the panel block edge is the hinge
end
blocks = cell(1, n_cs+1);
csIdx  = zeros(1, n_cs+1);
allcs  = [];
for j = 1:n_cs
    blocks{j+1} = geom.cspanels{j};
    csIdx(j+1)  = j;
    allcs = [allcs, geom.cspanels{j}(:)']; %#ok<AGROW>
end
blocks{1} = setdiff(1:geom.npan, allcs);
csIdx(1)  = 0;
keep = ~cellfun(@isempty, blocks);
blocks = blocks(keep);  csIdx = csIdx(keep);
B = repmat(empty_body(), 1, numel(blocks));
for b = 1:numel(blocks)
    pnl = blocks{b};
    ids = reshape(4*(pnl(:)'-1) + (1:4)', [], 1);
    B(b).x  = xcn(ids);
    B(b).y  = ycn(ids);
    B(b).a0 = zeros(numel(ids),1);
    B(b).faces = reshape(1:4*numel(pnl), 4, [])';
    B(b).faces = B(b).faces(:, end:-1:1);       % outward normals after the z flip
    B(b).lines = block_outline(pnl, xcn, ycn);
    B(b).cs    = csIdx(b);
    if csIdx(b) > 0, B(b).lw = 1.0; else, B(b).lw = 1.4; end
    B(b).hasCaps = false;
    B(b).Na = 0;  B(b).Ns = 0;  B(b).m = 0;
end
end

% =========================================================================
function loops = block_outline(pnl, xcn, ycn)
%BLOCK_OUTLINE  Closed boundary polylines (local vertex indices) of a panel block.
np = numel(pnl);
ids = reshape(4*(pnl(:)'-1) + (1:4)', [], 1);
key = round([xcn(ids), ycn(ids)], 6);
[~, ia, ic] = unique(key, 'rows');
ed = zeros(4*np, 2);
for i = 1:np
    n4 = ic(4*(i-1) + (1:4));  n4 = n4(:);
    ed(4*(i-1)+(1:4), :) = [n4([1 2]), n4([2 3]), n4([3 4]), n4([4 1])]';
end
[ue, ~, ie] = unique(sort(ed, 2), 'rows');
cnt = accumarray(ie, 1);
bnd = ue(cnt == 1, :);
loops = {};
used = false(size(bnd,1), 1);
while ~all(used)
    i0 = find(~used, 1);
    used(i0) = true;
    seq = bnd(i0,:);
    cur = seq(2);
    while cur ~= seq(1)
        nxt = find(~used & any(bnd == cur, 2), 1);
        if isempty(nxt), break, end
        used(nxt) = true;
        e = bnd(nxt,:);
        if e(1) == cur, cur = e(2); else, cur = e(1); end
        seq(end+1) = cur; %#ok<AGROW>
    end
    loops{end+1} = reshape(ia(seq), 1, []); %#ok<AGROW>
end
end

% =========================================================================
function [B, hingeY, hingeX] = build_bodies_foil(geom, isTEs, flapChord, slatChord, ...
    zupf, zlof, g_x, g_y)
%BUILD_BODIES_FOIL  Lofted skin: centre box, one body per surface, stubs.
span = geom.span;  xLE = geom.xLE;  cch = geom.c;
n_cs = numel(geom.hinges);
eps_tip = 1e-9;
hasTE = any(isTEs);  hasLE = any(~isTEs);
if hasTE, xiTE = 1 - flapChord; else, xiTE = 1; end
if hasLE, xiLE = slatChord;     else, xiLE = 0; end
if hasTE && hasLE && xiLE + 2*g_x >= xiTE
    error('animate_wing:chords', ...
        'FlapChord %.3g + SlatChord %.3g leave no centre body.', flapChord, slatChord);
end

% ---- chordwise stations --------------------------------------------------
xiF = linspace(xiTE + g_x, 1, 12);                         % TE surfaces
xiS = (xiLE - g_x)*(1 - cos(pi/2*(0:17)/17));              % LE surfaces
L_bn = 0.05;                                               % bullnose blend
xiCend = xiTE;
if hasTE, xiCend = xiTE - g_x; end
if hasLE
    xi0C = xiLE + g_x;
    xiC  = unique([xi0C + L_bn*(1 - cos(pi/2*(0:6)/6)), ...
                   linspace(xi0C + L_bn + 0.03, xiCend, 20)]);
    bnf  = @(xi) sqrt(1 - (1 - min(max((xi - xi0C)/L_bn, 0), 1)).^2);
    zupC = @(xi) 0.5*(zupf(xi)+zlof(xi)) + bnf(xi).*(0.5*(zupf(xi)-zlof(xi)));
    zloC = @(xi) 0.5*(zupf(xi)+zlof(xi)) - bnf(xi).*(0.5*(zupf(xi)-zlof(xi)));
    if hasTE, endC = 'none'; else, endC = 'last'; end
else
    xiC  = xiCend*(1 - cos(pi/2*(0:26)/26));               % LE-clustered, sharp LE
    zupC = zupf;  zloC = zlof;
    endC = 'first';
end

% ---- centre body ---------------------------------------------------------
ysC = linspace(0, span - eps_tip, 41);
B = repmat(empty_body(), 1, 1);
B(1) = fill_body(tube_body(xiC, endC, ysC, xLE, cch, zupC, zloC), 0, 0.9);
hingeY = zeros(n_cs, 2);
hingeX = zeros(1, n_cs);                    % the drawn cut is the render hinge

% ---- one body per control surface ---------------------------------------
for j = 1:n_cs
    ylo = min(geom.hinges(j).y1, geom.hinges(j).y2);
    yhi = max(geom.hinges(j).y1, geom.hinges(j).y2);
    ys  = linspace(max(ylo + g_y/2, 0), min(yhi - g_y/2, span - eps_tip), 9);
    if isTEs(j)
        bd = tube_body(xiF, 'last', ys, xLE, cch, zupf, zlof);
        hingeY(j,:) = [ys(1), ys(end)];                    % axis root -> tip
        hingeX(j)   = xLE - xiTE*cch;
    else
        bd = tube_body(xiS, 'first', ys, xLE, cch, zupf, zlof);
        hingeY(j,:) = [ys(end), ys(1)];                    % axis tip -> root
        hingeX(j)   = xLE - xiLE*cch;
    end
    B(end+1) = fill_body(bd, j, 0.9); %#ok<AGROW>
end

% ---- stubs: the parts of an edge strip that no surface covers ------------
for pass = 1:2
    if pass == 1
        sel = isTEs;  xiB = xiF;  endB = 'last';
    else
        sel = ~isTEs; xiB = xiS;  endB = 'first';
    end
    if ~any(sel), continue, end
    iv = zeros(0,2);
    for j = find(sel)
        iv(end+1,:) = sort([geom.hinges(j).y1, geom.hinges(j).y2]); %#ok<AGROW>
    end
    iv = sortrows(iv);
    y0 = 0;
    for i = 1:size(iv,1)
        if iv(i,1) - y0 > g_y
            B(end+1) = stub_body(xiB, endB, y0, iv(i,1), span, xLE, cch, ...
                zupf, zlof, g_y, eps_tip); %#ok<AGROW>
        end
        y0 = max(y0, iv(i,2));
    end
    if span - y0 > g_y
        B(end+1) = stub_body(xiB, endB, y0, span, span, xLE, cch, ...
            zupf, zlof, g_y, eps_tip); %#ok<AGROW>
    end
end
end

% =========================================================================
function b = stub_body(xiB, endB, ylo, yhi, span, xLE, cch, zupf, zlof, g_y, eps_tip)
%STUB_BODY  Fixed (non-rotating) part of an edge strip, drawn in the main colour.
ys = linspace(max(ylo + g_y/2, 0), min(yhi - g_y/2, span - eps_tip), 9);
b  = fill_body(tube_body(xiB, endB, ys, xLE, cch, zupf, zlof), 0, 0.9);
end

% =========================================================================
function b = fill_body(bd, csj, lw)
%FILL_BODY  Wrap a tube_body result as a render body (local faces plus caps).
b = empty_body();
b.x  = bd.x;   b.y = bd.y;   b.a0 = bd.a0;
nv   = bd.Na*bd.Ns;
rr   = (1:bd.Na)';  rrn = mod(rr, bd.Na) + 1;
Fcap1 = [repmat(nv+1, bd.Na, 1), rr, rrn];
Fcap2 = [repmat(nv+2, bd.Na, 1), (bd.Ns-1)*bd.Na + rrn, (bd.Ns-1)*bd.Na + rr];
F = [bd.faces; Fcap1; Fcap2];
b.faces = F(:, end:-1:1);                  % outward normals after the z flip
b.lines = { [1:bd.Na, 1], ...                                  % inboard rib
            (bd.Ns-1)*bd.Na + [1:bd.Na, 1], ...                % outboard rib
            (0:bd.Ns-1)*bd.Na + 1, ...                         % upper front
            (0:bd.Ns-1)*bd.Na + bd.m };                        % upper rear
b.cs = csj;  b.lw = lw;  b.hasCaps = true;
b.Na = bd.Na;  b.Ns = bd.Ns;  b.m = bd.m;
end

% =========================================================================
function b = empty_body()
b = struct('x',[], 'y',[], 'a0',[], 'faces',[], 'lines',{{}}, 'cs',0, ...
    'lw',1, 'hasCaps',false, 'Na',0, 'Ns',0, 'm',0, 'idx',[]);
end

% =========================================================================
function bd = tube_body(xiU, sharedEnd, ys, xLE, c_geo, zupf, zlof)
%TUBE_BODY  Closed ring-section tube: upper surface forward, lower back.
%   sharedEnd 'first': sharp point at xiU(1) (leading-edge body);
%   'last': sharp point at xiU(end) (trailing-edge body);
%   'none': open chordwise cut at both ends (centre body). The ring wraps
%   across the open cut(s), so the flat base faces fall out of the tube
%   triangulation automatically.
m = numel(xiU);
switch sharedEnd
    case 'first'
        xir = [xiU, fliplr(xiU(2:end))];
        up  = [true(1,m), false(1,m-1)];
    case 'last'
        xir = [xiU, fliplr(xiU(1:end-1))];
        up  = [true(1,m), false(1,m-1)];
    otherwise
        xir = [xiU, fliplr(xiU)];
        up  = [true(1,m), false(1,m)];
end
Na = numel(xir);
Ns = numel(ys);
a0r = zeros(Na,1);
a0r(up)  = zupf(xir(up))*c_geo;             % upward offset (physically up)
a0r(~up) = zlof(xir(~up))*c_geo;
x1 = xLE - xir(:)*c_geo;                    % structural x points to the LE
bd.x  = repmat(x1, Ns, 1);
bd.y  = reshape(repmat(ys(:)', Na, 1), [], 1);
bd.a0 = repmat(a0r, Ns, 1);
[Rg, Sg] = meshgrid(1:Na, 1:Ns-1);
Rg = Rg(:);  Sg = Sg(:);  Rn = mod(Rg,Na) + 1;
v00 = (Sg-1)*Na + Rg;  v01 = (Sg-1)*Na + Rn;
v10 =  Sg*Na    + Rg;  v11 =  Sg*Na    + Rn;
bd.faces = [v00 v01 v11; v00 v11 v10];
bd.Na = Na;  bd.Ns = Ns;  bd.m = m;
end

% =========================================================================
function img = grab_frame(fig, tmppng)
%GRAB_FRAME  Print at 2x and box-filter down to the nominal frame size.
print(fig, tmppng, '-dpng', '-r192');
img = imread(tmppng);
img = img(1:2*floor(end/2), 1:2*floor(end/2), :);
Fd  = double(img);
Fd  = (Fd(1:2:end,1:2:end,:) + Fd(2:2:end,1:2:end,:) + ...
       Fd(1:2:end,2:2:end,:) + Fd(2:2:end,2:2:end,:))/4;      % 2x2 box AA
img = uint8(Fd);
end

% =========================================================================
function box = union_crop(frames, bg, pad)
%UNION_CROP  Bounding box of every non-background pixel of all frames, plus pad.
msk = false(size(frames{1},1), size(frames{1},2));
bgv = reshape(255*bg(:)', 1, 1, 3);
for k = 1:numel(frames)
    msk = msk | (sum(abs(double(frames{k}) - bgv), 3) > 25);
end
rows = find(any(msk,2));  cols = find(any(msk,1));
if isempty(rows)
    box = [1, size(msk,1), 1, size(msk,2)];
    return
end
box = [max(1,rows(1)-pad), min(size(msk,1),rows(end)+pad), ...
       max(1,cols(1)-pad), min(size(msk,2),cols(end)+pad)];
end

% =========================================================================
function bytes = write_gif(frames, file, delay, loop, pf)
%WRITE_GIF  Animated GIF with ONE 256-colour palette (no palette flicker), no dithering:
% dither noise on the shaded skin doubles the file size and the 256 colours cover it.
[~, map] = rgb2ind(cat(1, frames{pf}), 256, 'nodither');
for k = 1:numel(frames)
    im = rgb2ind(frames{k}, map, 'nodither');
    if k == 1
        imwrite(im, map, file, 'gif', 'LoopCount',loop, 'DelayTime',delay);
    else
        imwrite(im, map, file, 'gif', 'WriteMode','append', 'DelayTime',delay);
    end
end
d = dir(file);
bytes = d.bytes;
end
