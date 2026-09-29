% make_virtual_flight_gif  Virtual flight of the Goland wing through its flutter
% boundary: writes docs/figures/virtual_flight.gif (light theme only, no mp4).
%
% A genuine time-varying simulation (build_LPV_P_Goland + lsim, the flight of
% docs/TUTORIAL.md chapter 11 on a compressed ramp): the freestream
% velocity ramps from 163 to 184 m/s in 1.8 s while a multisine disturbance
% (light turbulence) shakes the wing. Below the 165 m/s flutter boundary the wing
% only buffets; beyond it the flutter mode grows; when the tip deflection exceeds
% 8 cm the AFS controller from data/Goland_Cont_nO4_lpf.mat switches on (the
% label flips from AFS OFF to AFS ON, the flap turns blue) and damps the motion back to buffet level
% while the velocity keeps rising. Every frame comes from one continuous
% integration; only the ramp is compressed and the playback is slow motion
% (on-screen tag).
%
% Rendering: util/animate_wing with a NACA 0012 envelope, a drawing aid over the
% flat-plate DLM; tip deflection and flap throw are exaggerated (13 % semispan,
% 26 deg). 198 frames, 880 x 403 px, about 4.4 MB. Runtime about 40 s on R2026a.
% Run from the repo root after startup.m.

t_run = tic;
outdir = fileparts(mfilename('fullpath'));

num_modes = 5;
num_poles = 6;
[Structure, ~] = define_Goland_Structure_Aero(num_modes, num_poles);
G_act = define_PT2Actuator(1);
load('Goland_Cont_nO4_lpf.mat')                % Goland_Cont_nO (AFS controller)
LPV_P = build_LPV_P_Goland(100, 1, 1, 1:num_modes);

% --- geometry from the model definition ------------------------------------
Ps  = Structure.Ps;
xLE = max(cellfun(@(p) p{1}(1), Ps));          % structural x points upstream
xTE = min(cellfun(@(p) p{2}(1), Ps));
span  = max(cellfun(@(p) p{4}(2), Ps));
Kff = Structure.Kff;
q_f_scale = 1./(sqrt(diag(Kff))/sqrt(Kff(1,1)));   % LPV_P modal outputs are scaled

% flutter frequency and boundary (for timing and the readout color)
V_bound = 165;                                 % flutter boundary (interpolated)
e175 = eig(psample(LPV_P, [], 175));
cand = e175(real(e175) > 0 & imag(e175) > 20);
[~, imax] = max(real(cand));
f_flut = imag(cand(imax))/(2*pi);

% --- condensed virtual flight (genuine LTV sim, compressed ramp) -----------
fpc  = 10;                                     % frames per flutter cycle
sub  = 8;                                      % sim steps per frame
dt   = 1/(f_flut*fpc*sub);
T    = 1.8;                                    % flight duration (s), ~20 cycles
t_sim = (0:dt:T)';  Nt = numel(t_sim);
V0 = 163;  V1 = 184;
V_t = V0 + (V1-V0)*t_sim/T;

rng(7)                                         % reproducible turbulence
n_sin = 12;  sigma_d = 0.15;
f_sin = 3 + 17*rand(n_sin,2);  ph_sin = 2*pi*rand(n_sin,2);
uOL = zeros(Nt, num_modes+2);                  % [dist_q_f1..5, noise, flap1_d]
uOL(:,1) = sigma_d/sqrt(n_sin)*sum(sin(2*pi*t_sim*f_sin(:,1)' + ph_sin(:,1)'),2);
uOL(:,2) = sigma_d/sqrt(n_sin)*sum(sin(2*pi*t_sim*f_sin(:,2)' + ph_sin(:,2)'),2);

% tip deflection map (true LE/TE at the tip)
PHIz_tip = build_IMU([xLE; xTE], [span; span]-1e-9, Structure.E, Structure.ele);
Ctip = PHIz_tip*Structure.PHIgf;

% open loop until the tip deflection crosses the threshold
[yOL,~,xOL] = lsim(LPV_P, uOL, t_sim, [], V_t);
z_OL = (yOL(:,1:num_modes).*q_f_scale') * Ctip(1,:)';
z_thresh = 0.08;                               % trigger threshold (m)
i_on = find(abs(z_OL) > z_thresh & V_t > V_bound, 1);
assert(~isempty(i_on) && t_sim(i_on) < 0.8*T, 'trigger too late - retune sigma_d/z_thresh')

% closed loop from the switch-on state (u = K*y via lft; plant states first)
CL = lft(LPV_P, Goland_Cont_nO);
x_on = [xOL(i_on,:)'; zeros(order(Goland_Cont_nO),1)];
[yCL,~] = lsim(CL, uOL(i_on:end,1:num_modes+1), t_sim(i_on:end), x_on, V_t(i_on:end));

% stitch histories: modal coords (descaled), flap demand, measured accel
qrec = [yOL(1:i_on-1,1:num_modes); yCL(:,1:num_modes)].*q_f_scale';
urec = [zeros(i_on-1,1); yCL(:,2*num_modes+1)];          % flap1_d demand (rad)
% reconstruct the measured IMU acceleration by replaying the plant open loop
% with the recorded flap demand (same trajectory up to the FOH interpolation
% of the demand samples, but with the u_z1_ddot output exposed)
uRe = [uOL(:,1:num_modes+1), urec];
[yRe,~] = lsim(LPV_P, uRe, t_sim, [], V_t);
acc_rec = yRe(:,2*num_modes+2);                          % scaled u_z1_ddot
q_err = max(abs(yRe(:,1)*q_f_scale(1) - qrec(:,1))) / max(abs(qrec(:,1)));
fprintf('replay consistency (FOH interpolation error): %.2e relative\n', q_err);
assert(q_err < 0.05, 'closed-loop replay mismatch')

% actual flap deflection: demand through the PT2 actuator (output 1: angle)
flap_act = lsim(G_act{1}(1,1), urec, t_sim);

i_b = find(V_t > V_bound, 1);
fprintf(['Virtual flight: %.1f s = %.0f cycles at %.2f Hz | buffet %.0f%%, ' ...
    'growth %.0f%%, recovery %.0f%% of flight\n'], T, T*f_flut, f_flut, ...
    100*i_b/Nt, 100*(i_on-i_b)/Nt, 100*(Nt-i_on)/Nt);
fprintf('trigger at t = %.2f s, V = %.1f m/s | peak |z_tip| %.3f m, peak flap %.1f deg (physical)\n', ...
    t_sim(i_on), V_t(i_on), max(abs(z_OL(1:i_on))), max(abs(flap_act))*180/pi);

% --- render -----------------------------------------------------------------
frameIdx = 1:sub:Nt;
Nf = numel(frameIdx);
DelayTime = 0.05;
slowmo = DelayTime/(dt*sub);

q    = qrec(frameIdx,:)';
flap = flap_act(frameIdx)';
acc  = acc_rec(frameIdx)';
Active = frameIdx >= i_on;
Readout = arrayfun(@(k) sprintf('$V_\\infty = %5.1f$ m/s', V_t(k)), frameIdx, ...
    'UniformOutput', false);
Footnote = sprintf('%.0fx slow motion', slowmo);
S = fw_style();

[~, ~, info] = animate_wing(Structure, q, 'Surfaces', flap, 'IMUAccel', acc, ...
    'Active', Active, 'Readout', Readout, 'ReadoutColor', S.muted, ...
    'LabelOff', 'AFS OFF', 'LabelOn', 'AFS ON', ...
    'Footnote', Footnote, 'Skin', 'foil', 'NACA', '0012', ...
    'TipAmplitude', 0.13, 'SurfaceAmplitude', deg2rad(26), ...
    'Output', 'gif', 'GifFile', fullfile(outdir, 'virtual_flight.gif'), ...
    'FrameDelay', DelayTime);
assert(info.gifBytes < 7e6, 'virtual_flight.gif exceeds 7 MB')

fprintf('virtual_flight.gif: %d frames, %d x %d px, %.0f kB | total runtime %.0f s\n', ...
    Nf, info.size_px(2), info.size_px(1), info.gifBytes/1024, toc(t_run));
