%MAKE_README_FIGURES  Render the README figures: docs/figures/quickstart_vg.png and afs_hero.gif.
%   Run from the repository root after startup.m; about a minute (54 s measured) on R2026a. Light theme only, no mp4.
%   Section 1 repeats the six quick-start lines of the README verbatim and exports the V-g figure
%   (fw_figure standard size, 200 dpi). Section 2 renders the hero scene of QUICKSTART.m section 6
%   with the constants of the simulate_afs_switch header example (RectWing at 108 m/s, 3.7 m/s past
%   the 104.3 m/s onset, 6 open-loop cycles then 8 closed-loop cycles at 12 frames per cycle, NACA
%   2408 skin, 7 % semispan tip amplitude, 33 deg flaps and 25 deg slats) closed by the shipped
%   static gain RectWing_Cont, writes the 168-frame GIF, asserts that it stays under 7 MB and prints
%   the per-surface activity, so the "mainly works ..." claim of the README caption is checked on
%   every run. The GIF adds the velocity readout and the "deflections exaggerated" footnote that
%   the figure-window version in QUICKSTART.m does not need.
%   See also QUICKSTART, simulate_afs_switch, animate_wing, Vg_plot, fw_export.

t_run = tic;
outdir = fileparts(mfilename('fullpath'));

%% 1. Quick start (the README lines, verbatim) -> quickstart_vg.png
V = linspace(20,160,32);  G = build_G_RectWing(V, 1:8, 1:8);   % wing model: 8 IMUs, 8 control surfaces, one ss per velocity
EV = getEigenvalueModeshape(G, 5);                             % aeroelastic modes tracked over velocity
load RectWing_Cont_imu18_ail18.mat                             % shipped AFS controller: static gain, 8 IMUs -> 8 surfaces
EV_cl = getEigenvalueModeshape(feedback(G, -RectWing_Cont), 5);   % u = K*y
[fig, cr] = Vg_plot({V,V}, {EV,EV_cl}, {'open loop','closed loop'}, 'Band',[2 10]);   % flutter at about 104.3 m/s; none with the controller
fw_export(fig, fullfile(outdir, 'quickstart_vg.png'));
close(fig)

%% 2. Hero scene (QUICKSTART.m section 6) -> afs_hero.gif
[Structure, ~] = define_RectWing_Structure_Aero(5, 6);
G108 = build_G_RectWing(108, 1:8, 1:8);
sim = simulate_afs_switch(G108, RectWing_Cont, 'GrowthCycles', 6, 'ControlCycles', 8, 'FramesPerCycle', 12);
gifFile = fullfile(outdir, 'afs_hero.gif');
Readout = repmat({'$V_\infty = 108$ m/s'}, 1, sim.Nt);
S = fw_style();
[~, ~, info] = animate_wing(Structure, sim.q_f, 'Surfaces', sim.delta, 'IMUAccel', sim.acc, 'Active', sim.active, ...
    'Skin', 'foil', 'NACA', '2408', 'FlapChord', 0.25, 'SlatChord', 0.20, ...
    'TipAmplitude', 0.07, 'SurfaceAmplitude', deg2rad([33 25]), ...
    'Readout', Readout, 'ReadoutColor', S.muted, 'LabelOff', 'AFS OFF', 'LabelOn', 'AFS ON', ...
    'Footnote', 'deflections exaggerated', ...
    'Output', 'gif', 'GifFile', gifFile, 'Width', 900, 'CropPad', 16);
d = dir(gifFile);
assert(d.bytes < 7e6, 'afs_hero.gif exceeds 7 MB: reduce FramesPerCycle to 10 before shrinking Width')

% per-surface activity in the hero scene: the README caption quotes it
a = rms(sim.delta, 2);  a = a/max(a);
fprintf('Hero activity (rms surface angle, normalised to the largest):\n')
for i = 1:numel(a), fprintf('  %-8s %5.2f\n', G108.InputName{i}, a(i)); end
fprintf('Controller off: +%.1f%% per cycle; on: %.1f%% per cycle, %.0f%% of the peak left after 8 cycles\n', ...
    100*(sim.growthPerCycle-1), 100*(sim.decayPerCycle-1), 100*sim.residual)
fprintf('afs_hero.gif: %d frames, %d x %d px, %.0f kB; quickstart_vg.png written; total %.0f s\n', ...
    sim.Nt, info.size_px(2), info.size_px(1), d.bytes/1024, toc(t_run))
