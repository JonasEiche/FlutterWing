%% R05: Quadrature Mismatch Delta Sweep on the Full FlutterWing Model
%
% Validates the quadrature mismatch mechanism (derived analytically in R04
% on a 2x2 toy plant) on the full-order FlutterWing aeroservoelastic model
% with 8 accelerometers, 4 trailing-edge flaps, and 4 leading-edge slats.
%
%
% EXPERIMENTAL DESIGN
% -------------------
% The script sweeps the quadrature mismatch angle Delta over a uniform
% grid 0:5:170 degrees plus the exact Delta values of two H2-optimal
% SISO methods (injected for fair comparison). At each Delta:
%
%   1. Blending vectors k_y and k_u are constructed by rotating within
%      the modal coupling subspaces of the critical flutter mode.
%      The output direction theta_y is fixed at the H2-opt separate
%      angle; only the input direction theta_u = theta_y + Delta varies.
%
%   2. A 4th-order SISO controller is synthesised via systune with:
%      - Hard constraint: gain margin >= 6 dB, phase margin >= 45 deg
%      - Soft goals: modal velocity attenuation (q_f1, q_f2), control
%        effort under disturbance, control effort under sensor noise
%      - Weighting: Theis-style frequency shaping for actuator limits
%
%   3. The maximum soft goal (worst-case among the three objectives)
%      and the hard goal satisfaction are recorded.
%
% The H2-optimal methods go through the identical synthesis pipeline,
% ensuring any performance difference is due to the blending direction
% alone, not controller tuning.
%
% MODAL COUPLING SUBSPACES
% ------------------------
% The flutter mode's coupling to sensors and actuators is characterised
% by the modal input/output matrices:
%
%       B_modal = W^H * B     (left eigenvector projection of input matrix)
%       C_modal = C * V       (right eigenvector projection of output matrix)
%
% For a complex eigenvalue lambda = sigma + j*omega with eigenvector v,
% the real and imaginary parts of C_modal(:, idx) span a 2D subspace in
% measurement space. Similarly for B_modal(idx, :)^T in actuator space.
% These are normalised to form orthonormal bases:
%
%       [cvr, cvi] = normalised [Re(C*v), Im(C*v)]   (output coupling)
%       [bwr, bwi] = normalised [Re(B'*w), Im(B'*w)] (input coupling)
%
% Any blending vector in this subspace can be parameterised by a single
% angle (theta_y for output, theta_u for input). The mismatch
% Delta = theta_u - theta_y governs the zero placement exactly as in R04,
% up to perturbations from the non-targeted modes.
%
%
% PLANT CONSTRUCTION
% ------------------
% Two plant models are used:
%
%   G = build_G_RectWing(V_inf_ref, imuIDX, ailIDX)
%       Output-only plant at V_inf_ref = 90 m/s (stable regime).
%       Used for eigenanalysis, modal coupling extraction, and H2
%       blending vector computation. The stable reference speed ensures
%       well-defined eigenvectors and avoids numerical issues with
%       near-zero damping at the flutter speed.
%
%   P = build_P_RectWing(V_inf, imuIDX, ailIDX, modesIDX)
%       Full performance plant at V_inf = 130 m/s (above flutter speed).
%       Includes disturbance and noise input channels, modal velocity
%       performance outputs, and the generalized plant structure needed
%       for systune. modesIDX = 1:2 retains the first two structural
%       modes (first bending and first torsion).
%
% CONTROLLER SYNTHESIS DETAILS
% ----------------------------
%   - Controller order: nO = 4 (4th order, as in R06_Locked_Comparison)
%   - Stability margins: TuningGoal.Margins('ud', 6, 45) at plant input
%   - Performance goals (soft):
%       ReqAttenQF:      disturbance -> modal velocities, gain < Vp/Vd
%       ReqContEffNoise: sensor noise -> actuator, gain < invWu*Vu/Vn
%       ReqContEffDist:  disturbance -> actuator, gain < invWu*Vu/Vd
%   - Actuator weighting: Theis band-pass shape invW_theis with
%       corner frequencies w1 = 12, w2 = 64 rad/s
%
% SCRIPT STRUCTURE
% ----------------
%   S1: Configuration & Model Setup
%       Defines sensor/actuator indices (8 IMUs, 8 control surfaces),
%       builds G at V_inf_ref = 90 m/s and P at V_inf = 130 m/s.
%
%   S2: Flutter Mode Extraction, Modal Coupling & H2 Blending Vectors
%       Eigendecomposition of G.A to identify the flutter mode (sigma > -1,
%       5 < omega < 40). Extracts modal input/output coupling matrices,
%       normalises to form orthonormal bases. Computes H2 Blending and
%       H2-opt separate blending vectors, projects into modal basis to
%       determine Delta.
%
%   S3: Systune Delta Sweep
%       Sweeps Delta over the grid, synthesises controllers, records
%       soft/hard goal values. H2 methods are injected at their exact
%       display-Delta values for direct comparison.
%
%   S4: Summary
%       Reports best sweep point and H2 method performance. Identifies
%       the optimal Delta range for the FlutterWing.
%
%   S5: Save Results
%       Saves all sweep data to Data/Quadrature_Mismatch_DeltaSweep.mat
%       for use by fig06_Quadrature_Mismatch_Sweep.m.
%
% OUTPUTS
% -------
%   Data/Quadrature_Mismatch_DeltaSweep.mat containing:
%       - V_inf, V_inf_ref: operating and reference velocities
%       - Delta_synth_deg: sweep grid in degrees
%       - synth_softGoal: [nSynth x 3] soft goals (attenuation, noise CE, dist CE)
%       - synth_hardGoal: [nSynth x 1] hard margin constraint
%       - synth_maxSoft: [nSynth x 1] max soft goal per Delta (primary metric)
%       - h2_names: {'H2 Blending', 'H2-opt sep.'}
%       - h2_Delta_deg: raw Delta of H2 methods (may be negative)
%       - h2_Delta_display: Delta mapped to [0, 180) for display
%       - h2_sweep_idx: indices into sweep grid for H2 methods
%
%   Figures generated by fig06_Quadrature_Mismatch_Sweep.m (separate script):
%       - fig06_flutterwing_delta_sweep: MaxSoft vs Delta with H2 marker
%
%
% DEPENDENCIES
% ------------
%   Toolboxes: Control System Toolbox (systune, TuningGoal, tunableSS, AnalysisPoint)
%   Functions: build_G_RectWing, build_P_RectWing, 
%              h2_opt_output_siso, h2_opt_input_siso
%
% SEE ALSO
%   R04_SISO_RHP_Zeros.m       - Analytical derivation on toy plant
%   fig06_Quadrature_Mismatch_Sweep.m - Paper-ready figure from R05 results
%   R06_Locked_Comparison.m - locked re-synthesis of the best member (Table 2)
%
% -------------------------------------------------------------------------

clearvars

%% ========================================================================
%  S1: Configuration & Model Setup
%  ========================================================================

%  -5-6-7-8-
% |         |
%  -1-2-3-4-
imuIDX = 1:8;
ailIDX = 1:8;
modesIDX = 1:2;

V_inf = 130;
V_inf_ref = 90; % stable regime for blending (as in R06_Locked_Comparison)
G = build_G_RectWing(V_inf_ref, imuIDX, ailIDX);
P = build_P_RectWing(V_inf, imuIDX, ailIDX, modesIDX);
nym = length(imuIDX);
nud = length(ailIDX);


%% ========================================================================
%  S2: Flutter Mode Extraction, Modal Coupling & H2 Blending Vectors
%  ========================================================================

[V_eig, DD] = eig(G.A);
flutIDX = find((real(diag(DD)) > -1) & (abs(imag(diag(DD))) < 40) & (abs(imag(diag(DD))) > 5));

Bm_flut = (V_eig \ G.B); Bm_flut = Bm_flut(flutIDX, :);
Cm_flut = G.C * V_eig;    Cm_flut = Cm_flut(:, flutIDX);

% Modal coupling vectors (normalized real/imag parts)
cvr = normalize(real(Cm_flut(:,1)));
cvi = normalize(imag(Cm_flut(:,1)));
bwr = normalize(real(Bm_flut(1,:))');
bwi = normalize(imag(Bm_flut(1,:))');

lambda_flut = DD(flutIDX(1), flutIDX(1));
fprintf('Flutter mode at V=%d: lambda = %.4f +/- %.4fj\n\n', ...
    V_inf_ref, real(lambda_flut), imag(lambda_flut));

% --- H2-Optimal SISO Blending Vectors ---

% Method 1: H2 Optimal Blending
% The code for the calculation of the H_2 optimal blending vectors is proprietary, hence hardcoded here:
ky_h2 = [-0.0995;-0.3201;-0.5322;-0.6997;0.0764;0.180;0.2131;0.1764];
ku_h2 = [0.0450;0.3120;0.4477;0.5603;0.1156;0.2680;0.3858;0.3902];

% Method 2: H2-optimal separate input/output
W_eig = inv(V_eig)';
w_flut = W_eig(:, flutIDX(1));
[~, h] = h2_opt_output_siso(G.C, V_eig(:, flutIDX(1)));
[~, d] = h2_opt_input_siso(G.B, w_flut);
ky_h2sep = h;
ku_h2sep = d;

% Project into modal basis to find Delta
ay = [cvr, cvi] \ ky_h2;
theta_y_h2 = atan2(ay(2), ay(1));
au = [bwr, bwi] \ ku_h2;
theta_u_h2 = atan2(au(2), au(1));
Delta_h2 = wrapToPiLocal(theta_u_h2 - theta_y_h2);

ay2 = [cvr, cvi] \ ky_h2sep;
theta_y_h2sep = atan2(ay2(2), ay2(1));
au2 = [bwr, bwi] \ ku_h2sep;
theta_u_h2sep = atan2(au2(2), au2(1));
Delta_h2sep = wrapToPiLocal(theta_u_h2sep - theta_y_h2sep);

h2_Delta_deg = [rad2deg(Delta_h2); rad2deg(Delta_h2sep)];
h2_names = {'$H_2$ Blending', '$H_2$-opt sep.'};

% Map to display-Delta in [0, 180)
h2_Delta_display = h2_Delta_deg;
h2_Delta_display(h2_Delta_display < 0) = h2_Delta_display(h2_Delta_display < 0) + 180;

fprintf('=== H2 Blending Vectors ===\n');
fprintf('%-18s %12s %12s\n', 'Method', 'Delta [deg]', 'Display [deg]');
fprintf('%s\n', repmat('-', 1, 44));
for i = 1:2
    fprintf('%-18s %12.1f %12.1f\n', h2_names{i}, h2_Delta_deg(i), h2_Delta_display(i));
end
fprintf('\n');

%% ========================================================================
%  S3: Systune Delta Sweep
%  ========================================================================

fprintf('=== Systune Delta Sweep ===\n');

% Build sweep grid: regular steps plus H2 display-Deltas
Delta_synth_deg = unique(sort([0:5:170, ...
                               h2_Delta_display(:)']));
nSynth = numel(Delta_synth_deg);

% Identify which sweep indices correspond to H2 methods
[~, h2_sweep_idx] = ismember(h2_Delta_display, Delta_synth_deg);

% Weighting functions
w1 = 12;
w2 = 64;
s = tf('s');
invW_theis = ((s+0.01*w1)*(0.01*s+w2)) / ((s+w1)*(s+w2));

Vp = 0.2;
Vn = 0.1;
Vd = 0.5;
Vu = 0.5;
invWu = invW_theis;

ReqMarg = TuningGoal.Margins('ud', 6, 45);

noise_u_z_ddot = {'noise_u_z1_ddot','noise_u_z2_ddot','noise_u_z3_ddot','noise_u_z4_ddot', ...
                  'noise_u_z5_ddot','noise_u_z6_ddot','noise_u_z7_ddot','noise_u_z8_ddot'};
flap_d_slat_d = {'flap1_d','flap2_d','flap3_d','flap4_d', ...
                 'slat1_d','slat2_d','slat3_d','slat4_d'};

ReqContEffDist  = TuningGoal.Gain({'dist_q_f1_ddot','dist_q_f2_ddot'}, flap_d_slat_d, invWu*Vu*(1/Vd));
ReqAttenQF      = TuningGoal.Gain({'dist_q_f1_ddot','dist_q_f2_ddot'}, {'q_f1','q_f2'}, Vp*(1/Vd));
ReqContEffNoise = TuningGoal.Gain(noise_u_z_ddot, flap_d_slat_d, invWu*Vu*(1/Vn));

nO = 4;
% theta_y_fixed = 0;
theta_y_fixed = theta_y_h2sep;
rngSeed = 20260305;
rng(rngSeed, 'twister');
opt = systuneOptions('RandomStart', 3, 'UseParallel', true, 'Display', 'off');

% Storage
synth_softGoal = nan(nSynth, 3);
synth_hardGoal = nan(nSynth, 1);
synth_maxSoft  = nan(nSynth, 1);

fprintf('%-12s %12s %12s %8s\n', 'Delta [deg]', 'MaxSoft', 'HardGoal', 'Status');
disp('------------------------------------------------')

for k = 1:nSynth
    Delta_k = Delta_synth_deg(k) * pi/180;
    theta_u_k = theta_y_fixed + Delta_k;

    ky_k = normalize(cos(theta_y_fixed)*cvr + sin(theta_y_fixed)*cvi);
    ku_k = normalize(cos(theta_u_k)*bwr + sin(theta_u_k)*bwi);

    % Systune
    tuneCont_k = ku_k * tunableSS('cont', nO, 1, 1) * ky_k';
    tuneCL_k = lft(P, AnalysisPoint('ud', nud) * tuneCont_k * AnalysisPoint('ym', nym));
    rng(rngSeed, 'twister');
    [~, fSoft_k, gHard_k] = systune(tuneCL_k, ...
        [ReqAttenQF, ReqContEffNoise, ReqContEffDist], ReqMarg, opt);

    synth_softGoal(k,:) = fSoft_k;
    synth_hardGoal(k) = gHard_k;
    synth_maxSoft(k) = max(fSoft_k);

    status_k = "OK";
    if gHard_k > 1, status_k = "FAIL"; end

    label = '';
    if k == h2_sweep_idx(1), label = '  <-- H2 Blending'; end
    if k == h2_sweep_idx(2), label = '  <-- H2-opt sep.'; end
    fprintf('%12.1f %12.4f %12.4f %8s%s\n', ...
        Delta_synth_deg(k), synth_maxSoft(k), gHard_k, status_k, label);
end
fprintf('\n');

%% ========================================================================
%  S4: Summary
%  ========================================================================

% Extract H2 method results from the sweep
h2_maxSoft  = synth_maxSoft(h2_sweep_idx);
h2_hardGoal = synth_hardGoal(h2_sweep_idx);

[bestSoft, bestIdx] = min(synth_maxSoft);

fprintf('=== Summary ===\n');
fprintf('%-12s %12s %12s %8s\n', 'Delta [deg]', 'MaxSoft', 'HardGoal', 'Status');
disp('------------------------------------------------')
for k = 1:nSynth
    status_k = "OK";
    if synth_hardGoal(k) > 1, status_k = "FAIL"; end

    label = '';
    if k == h2_sweep_idx(1), label = '  <-- H2 Blending'; end
    if k == h2_sweep_idx(2), label = '  <-- H2-opt sep.'; end
    fprintf('%-12.1f %12.4f %12.4f %8s%s\n', ...
        Delta_synth_deg(k), synth_maxSoft(k), synth_hardGoal(k), status_k, label);
end

fprintf('\nBest SISO at Delta=%.1f deg (MaxSoft=%.4f)\n', ...
    Delta_synth_deg(bestIdx), bestSoft);
for i = 1:2
    fprintf('%-18s display-Delta = %+.1f deg, MaxSoft = %.4f, HardGoal = %.4f\n', ...
        h2_names{i}, h2_Delta_display(i), h2_maxSoft(i), h2_hardGoal(i));
end
fprintf('\n');

% See fig06_Quadrature_Mismatch_Sweep.m for plots

scriptDir = fileparts(mfilename('fullpath'));

%% ========================================================================
%  S5: Save Results
%  ========================================================================

dataDir = fullfile(scriptDir, 'Data');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

results.V_inf = V_inf;
results.V_inf_ref = V_inf_ref;
results.Delta_synth_deg = Delta_synth_deg;
results.synth_softGoal = synth_softGoal;
results.synth_hardGoal = synth_hardGoal;
results.synth_maxSoft = synth_maxSoft;
results.h2_names = h2_names;
results.h2_Delta_deg = h2_Delta_deg;
results.h2_Delta_display = h2_Delta_display;
results.h2_sweep_idx = h2_sweep_idx;

save(fullfile(dataDir, 'Quadrature_Mismatch_DeltaSweep.mat'), 'results');
fprintf('Results saved to Data/Quadrature_Mismatch_DeltaSweep.mat\n');

%% ========================================================================
%  Helper Functions
%  ========================================================================

function x = normalize(x)
    x = x / norm(x);
end

function a = wrapToPiLocal(a)
    a = mod(a + pi, 2*pi) - pi;
end
