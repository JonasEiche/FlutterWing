%% Sensor Fault Tolerance: Progressive Failures with Structural Blending
%
%
%  Structural Blending Reconfigurability
%  --------------------------------------
%  The SB controller has the factored structure
%
%      $K = K_U \, K_\mathrm{int}(s) \, K_Y^\top$
%
%  where $K_U \in \mathbb{R}^{n_u \times n_c}$ and
%  $K_Y \in \mathbb{R}^{n_y \times n_c}$ are the static input and output
%  blending matrices, and $K_\mathrm{int}(s) \in \mathbb{R}^{n_c \times n_c}$
%  is the inner MIMO controller designed for the blended
%  $n_c$-input/$n_c$-output loop. The output blending matrix $K_Y$ is
%  derived from the output mode shape matrix $\Phi_{zf}$ via pseudoinverse:
%
%      $K_Y = [\Phi_{zf}(:,1{:}n_c)^+]^\top$
%
%  where $\Phi_{zf} \in \mathbb{R}^{n_y \times n_f}$ maps structural modes
%  to IMU sensor outputs, and $n_c = 2$ selects the flutter-critical
%  bending and torsion modes.
%
%
%  Requires: Data/controller_imu18_ail18_structural_blending.mat
%            (from R0_afs_structural_blending_imu18_ail18_synthesis.m)
%
%  Output: Data/FaultTolerance_SensorFailure.mat
%  ------  (read by Fig_Sensor_Fault_Tolerance_Region.m, which draws
%           Figures/Fig26_Sensor_Fault_Tolerance_Region.pdf)
%  V_SB_min, V_SB_max     — SB worst/best $V_\mathrm{flutter}$ per $k$  [max_fail x 1]
%  V_H2_min, V_H2_max     — $\mathcal{H}_2$ worst/best $V_\mathrm{flutter}$ per $k$  [max_fail x 1]
%  V_flutter_OL            — Open-loop $V_\mathrm{flutter}$              [scalar]
%  V_flutter_SB_nom        — Nominal SB $V_\mathrm{flutter}$             [scalar]
%  V_flutter_H2_nom        — Nominal $\mathcal{H}_2$ $V_\mathrm{flutter}$ [scalar]
%  ny ($n_y$)              — Number of IMU sensors                       [scalar]
%  n_modes_blend ($n_c$)   — Virtual channels (blending dimension)       [scalar]
%  max_fail                — Maximum simultaneous failures tested        [scalar]
%
%  NOTATION
%  =========================================================================
%  Symbol                    Code variable      Description
%  -------------------------------------------------------------------------
%  $\Phi_{zf}$               PHIzf              Output mode shape matrix
%  $\Phi_{zf}'$              PHIzf_red          Reduced mode shape matrix (failed rows removed)
%  $K_Y$                     ky_SB              Output blending (nominal)
%  $K_Y'$                    ky_full            Output blending (reconfigured)
%  $K_U$                     ku_SB              Input blending (unchanged under faults)
%  $K_\mathrm{int}(s)$       nOcont_SB          Inner controller (unchanged under faults)
%  $n_y = 8$                 ny                 Number of IMU sensors
%  $n_c = 2$                 n_modes_blend      Virtual channels (flutter modes)
%  $k$                       k                  Number of simultaneous sensor failures
%  $V_\mathrm{flutter}$      V_flutter_*        Flutter velocity (from eigenvalue analysis)

clearvars

script_path = mfilename('fullpath');
script_dir  = fileparts(script_path);
data_dir    = fullfile(script_dir, 'Data');

% ---- Load controllers ---------------------------------------------------
ContPath = fullfile(data_dir, 'controller_imu18_ail18_structural_blending.mat');
load(ContPath, 'cont_H2', 'cont_SB')

cont_SB_nom = cont_SB;   % 8x8 nominal SB controller
cont_H2_nom = cont_H2;   % 8x8 nominal H2 controller

% ---- Structural/Aero parameters ----------------------------------------
num_modes = 5;
num_poles = 6;
imuIDX = 1:8;
ailIDX = 1:8;
ny = length(imuIDX);
n_modes_blend = 2;  % number of flutter-critical modes targeted by blending

[Structure, Aero] = define_RectWing_Structure_Aero(num_modes, num_poles);
Sfj    = Structure.Sfj;
PHIgf  = Structure.PHIgf;
DRe_jx = Structure.DRe_jx;
DIm_jx = Structure.DIm_jx;
PHIzg  = Structure.PHIzg;
c_ref  = Aero.c_ref;
poles  = Aero.poles;
Q0jj   = Aero.Q0jj;
QLpjj  = Aero.QLpjj;
num_panels = size(Q0jj,1);
num_AIL = size(DRe_jx,2);

% ---- Nominal blending vectors -------------------------------------------
PHIzf = PHIzg * PHIgf;                       % 8 x 5
ky_SB = pinv(PHIzf(:,1:n_modes_blend))';     % 8 x 2

sumQLpjjBjx = zeros(num_panels, num_AIL);
for i = 1:num_poles
    sumQLpjjBjx = sumQLpjjBjx + QLpjj(:,:,i)*(DRe_jx - DIm_jx*poles(i)*2/c_ref);
end
Bgx_struct = Sfj*(Q0jj*DRe_jx + sumQLpjjBjx);
ku_SB = pinv(Bgx_struct(1:n_modes_blend,:)); % 8 x 2

% ---- Extract inner controller -------------------------------------------
nOcont_SB = pinv(ku_SB) * cont_SB_nom * pinv(ky_SB');   % 2x2 ss

% ---- Velocity sweep setup -----------------------------------------------
V_inf = linspace(20,160,32);
G = build_G_RectWing(V_inf,imuIDX,ailIDX);

find_flutter_V = @(EV, vel) find_flutter_velocity(EV, vel, num_modes);

% ---- Nominal flutter velocities ----------------------------------------
[EV_OL,~] = getEigenvalueModeshape(G, num_modes);
CL_SB_nom = feedback(G, -cont_SB_nom);
CL_H2_nom = feedback(G, -cont_H2_nom);
[EV_SB_nom,~] = getEigenvalueModeshape(CL_SB_nom, num_modes);
[EV_H2_nom,~] = getEigenvalueModeshape(CL_H2_nom, num_modes);

V_flutter_OL     = find_flutter_V(EV_OL,     V_inf);
V_flutter_SB_nom = find_flutter_V(EV_SB_nom, V_inf);
V_flutter_H2_nom = find_flutter_V(EV_H2_nom, V_inf);

fprintf('=== Nominal Flutter Velocities ===\n');
fprintf('Open Loop:  %.1f m/s\n', V_flutter_OL);
fprintf('SB Nominal: %.1f m/s\n', V_flutter_SB_nom);
fprintf('H2 Nominal: %.1f m/s\n', V_flutter_H2_nom);

% =========================================================================
%% Progressive Sensor Failure Analysis
% =========================================================================
max_fail = ny - n_modes_blend;  % 6 (need >= 2 sensors for blending)

V_SB_min = zeros(max_fail,1); V_SB_max = zeros(max_fail,1);
V_H2_min = zeros(max_fail,1); V_H2_max = zeros(max_fail,1);

for k = 1:max_fail
    combos = nchoosek(1:ny, k);
    nc = size(combos, 1);
    V_SB_k = zeros(nc, 1);
    V_H2_k = zeros(nc, 1);

    for c = 1:nc
        failed = combos(c, :);
        active = setdiff(1:ny, failed);

        % SB: recompute output blending vector for reduced sensor set
        PHIzf_red = PHIzf(active, :);
        ky_red = pinv(PHIzf_red(:,1:n_modes_blend))';
        ky_full = zeros(ny, n_modes_blend);
        ky_full(active, :) = ky_red;
        cont_SB_red = ku_SB * nOcont_SB * ky_full';

        CL = feedback(G, -cont_SB_red);
        [EV,~] = getEigenvalueModeshape(CL, num_modes);
        V_SB_k(c) = find_flutter_V(EV, V_inf);

        % H2: zero failed sensor channels (no reconfiguration)
        cont_H2_red = cont_H2_nom;
        cont_H2_red.B(:, failed) = 0;
        cont_H2_red.D(:, failed) = 0;

        CL = feedback(G, -cont_H2_red);
        [EV,~] = getEigenvalueModeshape(CL, num_modes);
        V_H2_k(c) = find_flutter_V(EV, V_inf);
    end

    V_SB_min(k) = min(V_SB_k); V_SB_max(k) = max(V_SB_k);
    V_H2_min(k) = min(V_H2_k); V_H2_max(k) = max(V_H2_k);

    fprintf('k=%d (%3d combos): SB [%.1f, %.1f], H2 [%.1f, %.1f] m/s\n', ...
        k, nc, V_SB_min(k), V_SB_max(k), V_H2_min(k), V_H2_max(k));
end

% =========================================================================
%% Save results
% =========================================================================
save(fullfile(data_dir, 'FaultTolerance_SensorFailure.mat'), ...
    'V_SB_min', 'V_SB_max', 'V_H2_min', 'V_H2_max', ...
    'V_flutter_OL', 'V_flutter_SB_nom', 'V_flutter_H2_nom', ...
    'ny', 'n_modes_blend', 'max_fail');

fprintf('\nResults saved to Data/FaultTolerance_SensorFailure.mat\n');


% =========================================================================
%% Local function
% =========================================================================
function V_flutter = find_flutter_velocity(EV, V_inf, num_modes)
    num_modes_disp = 2;
    EV_disp = EV(1:2*num_modes_disp,:);
    tol = 0.1; divtol = 1.0;
    flut_ind = (cumsum((real(EV_disp) > tol),2) == 1 & imag(EV_disp) > divtol);
    VEL = repmat(V_inf, size(EV_disp,1), 1);
    flut_v = VEL(flut_ind);
    if isempty(flut_v)
        V_flutter = V_inf(end);
    else
        V_flutter = min(flut_v);
    end
end
