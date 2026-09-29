function G = build_G_noAct_EnergyNormalized(V_inf, imuIDX, ailIDX)
% build_G_noAct_EnergyNormalized  Energy-normalized RectWing plant without actuator dynamics.
%
% Paper-local helper of the AS2026 folder, not part of build/. Used by
% R02_Modal_Observability_NoActuators.m.
%
% Simplified variant of build_G_RectWing_EnergyNormalized:
%   - No actuator dynamics: inputs are direct surface deflections (no _dot, _ddot)
%   - No output/input scaling: physical units preserved
%   - Energy normalization via similarity transform T on states:
%       q_f     states scaled by sqrt(Kff_ii)        (Strain Energy)
%       q_f_dot states scaled by sqrt(Mff_ii)        (Kinetic Energy)
%       x_L     states scaled by 1/sqrt(Kff_ii)      (Complementary Strain Energy)
%     This makes eigenvector components directly comparable across modes,
%     so visibility/actuatability metrics from modal decomposition are
%     absolute and decoupled.
%
% INPUT
%   V_inf       :       Vector of free stream Velocities
%   imuIDX      :       Indices of selected IMUs e.g. [4,8]
%   ailIDX      :       Indices of selected Ailerons  [4,8]
%
%   Flap [1234] / Slat [5678] Aileron & Acceleration Sensor Location:
%
%                        -5-6-7-8-
%                       |         |
%                        -1-2-3-4-
%
% OUTPUT
%   G(s)        :       [ny, nu, nv] Array of plant state space models

% --- Model Order Properties ----------------------------------------------
num_modes = 5;
num_poles = 6; % must be 6 to hit the RectWing RFA cache (define_RectWing_Structure_Aero); other values trigger a multi-minute DLM solve and warning FlutterWing:RFAcacheMiss
num_x_L   = num_modes * num_poles;

% --- HARD CODED ----------------------------------------------------------
num_AIL = 8;
num_IMU = 8;
% --- HARD CODED ----------------------------------------------------------

% --- Structure & Aero (no actuator) -------------------------------------
[Structure, Aero] = define_RectWing_Structure_Aero(num_modes, num_poles);

% --- Name Definition -----------------------------------------------------
assert(num_AIL == 8, 'Hard Coded 4 Slats & 4 Flaps');
InputName = cell(1, num_AIL);
for i = 1:4
    InputName{i}   = ['flap', num2str(i)];
    InputName{4+i} = ['slat', num2str(i)];
end
OutputName = cell(1, num_IMU);
for i = 1:num_IMU
    OutputName{i} = ['u_z', num2str(i), '_ddot'];
end
StateName = cell(1, 2*num_modes + num_x_L);
for i = 1:num_modes
    StateName{i}           = ['q_f', num2str(i)];
    StateName{num_modes+i} = ['q_f', num2str(i), '_dot'];
end
for i = 1:num_x_L
    StateName{2*num_modes+i} = ['aero_lag', num2str(i)];
end

% --- Energy Transformation -----------------------------------------------
% Diagonal similarity transform T such that ||x_tilde||_2 corresponds to
% the total energy of the system. Since Kff and Mff are strictly diagonal
% (from define_RectWing_Structure_Aero), the energy for each mode i is
% isolated and the transform is diagonal.
%
% 1. Structural States (q_f and q_f_dot)
%    Strain Energy:   U = 1/2 * q_f' * Kff * q_f
%    Kinetic Energy:  T = 1/2 * q_f_dot' * Mff * q_f_dot
%
%    Scaling:  q_f_tilde     = sqrt(Kff) * q_f
%              q_f_dot_tilde = sqrt(Mff) * q_f_dot
%
%    Then 1/2 * (q_f_tilde_i^2 + q_f_dot_tilde_i^2) is exactly the total
%    mechanical energy [J] in the i-th structural mode.
%
% 2. Aerodynamic Lag States (x_L)
%    The Roger RFA lag states have no intrinsic Hamiltonian energy — they
%    are mathematical poles representing the wake's memory. However, from
%    build_ABCD_G the generalized aero force from the lag states is
%    F_lag = D_til * x_L, where D_til = repmat(eye(num_modes),1,num_poles).
%    Thus x_L has the physical units of generalized force.
%
%    To give a force an energy-equivalent scaling, we use the complementary
%    strain energy — the elastic energy stored in the structure if this
%    aero force were applied statically:
%
%       U* = 1/2 * x_L' * Kff^{-1} * x_L
%
%    Scaling:  x_L_tilde = (1/sqrt(Kff)) * x_L
%
%    So x_L_tilde_i^2 reflects the "potential energy impact" of the i-th
%    lag state on the structure.
Kff_diag = diag(Structure.Kff);
Mff_diag = diag(Structure.Mff);

T_diag = [sqrt(Kff_diag);                             % q_f     - Strain Energy
          sqrt(Mff_diag);                              % q_f_dot - Kinetic Energy
          repmat(1 ./ sqrt(Kff_diag), num_poles, 1)];  % x_L     - Complementary Strain Energy

T     = diag(T_diag);
T_inv = diag(1 ./ T_diag);


% OutScale_diag = repmat(Kff_diag(1)/Mff_diag(1), num_IMU, 1);
OutScale_diag = ones(num_IMU, 1);
So_inv = diag(1 ./ OutScale_diag);

V_inf_ref = 50;
rho = Aero.rho;
q_bar_ref = 0.5*rho*V_inf_ref^2;
InScale_diag = repmat(1/(q_bar_ref), num_AIL,1);
Si = diag(InScale_diag);

% --- Build State Space Model Array ---------------------------------------
ny = length(imuIDX);
nu = length(ailIDX);
nv = length(V_inf);
G  = ss(zeros(ny, nu, nv));

for i_v = 1:nv
    [A_noS, B_noS, C_noS, D_noS] = build_ABCD_G(V_inf(i_v), Structure, Aero);

    % Energy-normalized state-space (only deflection columns from B, D)
    A_eng = T * A_noS * T_inv;
    B_eng = T * B_noS(:, 1:num_AIL) * Si;
    C_eng = So_inv * C_noS * T_inv;
    D_eng = So_inv * D_noS(:, 1:num_AIL) * Si;

    G_iv = ss(A_eng, B_eng(:, ailIDX), C_eng(imuIDX, :), D_eng(imuIDX, ailIDX));
    G_iv.InputName  = InputName(ailIDX);
    G_iv.OutputName = OutputName(imuIDX);
    G_iv.StateName  = StateName;

    G(:,:,i_v) = G_iv;
end
