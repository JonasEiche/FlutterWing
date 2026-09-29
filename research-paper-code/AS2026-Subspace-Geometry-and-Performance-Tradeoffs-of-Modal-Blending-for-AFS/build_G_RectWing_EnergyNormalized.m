function G = build_G_RectWing_EnergyNormalized(V_inf,imuIDX,ailIDX)
% build_G_RectWing_EnergyNormalized  Energy-normalized RectWing plant (with actuators) for modal pole-vector analysis.
%
% Paper-local helper of the AS2026 folder, not part of build/. Used by
% R01_Modal_Observability.m.
%
% Drop-in replacement for build_G_RectWing. Applies an additional similarity
% transform so that states are energy-normalized (kinetic, strain, and
% complementary strain energy). This makes eigenvector analysis physically
% meaningful without post-processing.
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
num_x_L         = num_modes*num_poles;

% --- HARD CODED ----------------------------------------------------------
V_inf_ref = 100;
num_AIL = 8;
num_IMU = 8;
% --- HARD CODED ----------------------------------------------------------

% --- Structure, Aero, Actuator PT2 Model ---------------------------------
[Structure, Aero] = define_RectWing_Structure_Aero(num_modes,num_poles);
[G_act] = define_RectWing_PT2Actuator(num_AIL);

% --- Extract fields for energy normalization -----------------------------
Sfj = Structure.Sfj;
DRe_jx = Structure.DRe_jx;
DIm_jx = Structure.DIm_jx;
Q0jj = Aero.Q0jj;
Q1jj = Aero.Q1jj;
QLpjj = Aero.QLpjj;
poles_aero = Aero.poles;

% --- Name Definition -----------------------------------------------------
assert(num_AIL==8,'Hard Coded 4 Slats & 4 Flaps');
InputNameNoAct = cell(1,3*num_AIL);
InputNameDesired = cell(1,num_AIL);
for i = 1:4
    InputNameDesired{i} = ['flap',num2str(i),'_d'];
    InputNameDesired{4+i} = ['slat',num2str(i),'_d'];

    InputNameNoAct{i} = ['flap',num2str(i)];
    InputNameNoAct{num_AIL+i} = ['flap',num2str(i),'_dot'];
    InputNameNoAct{2*num_AIL+i} = ['flap',num2str(i),'_ddot'];

    InputNameNoAct{4+i} = ['slat',num2str(i)];
    InputNameNoAct{num_AIL+4+i} = ['slat',num2str(i),'_dot'];
    InputNameNoAct{2*num_AIL+4+i} = ['slat',num2str(i),'_ddot'];
end
OutputNameNoAct = cell(1,num_IMU);
for i = 1:num_IMU
    OutputNameNoAct{i} = ['u_z',num2str(i),'_ddot'];
end
OutputNameDesired = cell(1,num_IMU);
for i = 1:num_IMU
    OutputNameDesired{i} = ['u_z',num2str(i),'_ddot'];
end
StateNameNoAct = cell(1,2*num_modes+num_x_L);
for i = 1:num_modes
    StateNameNoAct{i} = ['q_f',num2str(i)];
    StateNameNoAct{num_modes+i} = ['q_f',num2str(i),'_dot'];
end
for i = 1:num_x_L
    StateNameNoAct{2*num_modes+i} = ['aero_lag',num2str(i)];
end

% --- Scale & Energy Normalization ----------------------------------------

rho = Aero.rho;
c_ref = Aero.c_ref;
OMEGA = Structure.OMEGA;
q_bar_ref           = 0.5*rho*V_inf_ref^2;
Kff = Structure.Kff;
Mff = Structure.Mff;
Kff_diag = diag(Kff);

% Combined state scale: original StateScale * sqrt(energy_metric / Kff(1,1))
%   q_f:       1 * sqrt(Kff_ii / Kff(1,1))
%   q_f_dot:   OMEGA_i * sqrt(Kff_ii / Kff(1,1))
%   aero_lag:  (q_bar_ref*c_ref^2)^2 / sqrt(Kff_ii * Kff(1,1))
EnergyStateScale = [sqrt(Kff_diag / Kff(1,1));
                    OMEGA .* sqrt(Kff_diag / Kff(1,1));
                    repmat((q_bar_ref*c_ref^2)^2 ./ sqrt(Kff_diag * Kff(1,1)), num_poles, 1)];

% OutScale = repmat(Kff(1,1)/Mff(1,1), num_IMU, 1);
OutScale = ones(num_IMU, 1);    % No output scaling to keep sigmin contr and obsrv in @Modal_Obsrv_Contr at similar orders of magnitude

Kff_inv = inv(Kff);
w0_act = 32*2*pi;  % from define_RectWing_PT2Actuator
n_aero = 2*num_modes + num_x_L;

% --- Build State Space Model Array ---------------------------------------
aildddIDX = [ailIDX,num_AIL+ailIDX,2*num_AIL+ailIDX];
InputNameDesired = InputNameDesired(ailIDX);
OutputNameDesired = OutputNameDesired(imuIDX);
G_act = G_act(ailIDX);

ny = length(imuIDX);
nu = length(ailIDX);
nv = length(V_inf);
G = ss(zeros(ny,nu,nv));
for i_v = 1:length(V_inf)
    V_inf_iv = V_inf(i_v);
    [A_noS,B_noS,C_noS,D_noS] = build_ABCD_G(V_inf_iv,Structure,Aero);
    G_noAct = ss(diag(1./EnergyStateScale)*A_noS*diag(EnergyStateScale), ...
                 diag(1./EnergyStateScale)*B_noS, ...
                 diag(1./OutScale)*C_noS*diag(EnergyStateScale), ...
                 diag(1./OutScale)*D_noS);

    G_noAct.InputName = InputNameNoAct;
    G_noAct.OutputName = OutputNameNoAct;
    G_noAct.StateName = StateNameNoAct;

    G_noAct = G_noAct(imuIDX,aildddIDX);

    G_iv = connect(G_noAct,G_act{1:end},InputNameDesired,OutputNameDesired);

    % --- Energy normalization: actuator states (V_inf-dependent) ----------
    % Generalized aero forces from actuator deflection/rate
    num_panels = size(Q0jj, 1);
    sumQLpjjBjx = zeros(num_panels, num_AIL);
    for i = 1:length(poles_aero)
        sumQLpjjBjx = sumQLpjjBjx + QLpjj(:,:,i)*(DRe_jx - DIm_jx*poles_aero(i)*2/c_ref);
    end
    Bgx = 0.5*rho*V_inf_iv^2*Sfj*(Q0jj*DRe_jx + sumQLpjjBjx);

    sumQLpjjBjx_dot = zeros(num_panels, num_AIL);
    for i = 1:length(poles_aero)
        sumQLpjjBjx_dot = sumQLpjjBjx_dot + QLpjj(:,:,i)*DIm_jx;
    end
    Bgx_dot = 0.5*rho*Sfj*V_inf_iv*(Q0jj*DIm_jx + Q1jj*0.5*c_ref*DRe_jx + sumQLpjjBjx_dot);
    % Complementary strain energy for actuator states (selected ailerons)
    % Interleaved [delta_1, delta_dot_1/w0, delta_2, ...] after connect()
    q_diag_act_pos = zeros(length(ailIDX), 1);
    q_diag_act_vel = zeros(length(ailIDX), 1);
    for j = 1:length(ailIDX)
        F_delta     = Bgx(:, ailIDX(j));
        F_delta_dot = Bgx_dot(:, ailIDX(j));
        q_diag_act_pos(j) = F_delta' * Kff_inv * F_delta;
        q_diag_act_vel(j) = w0_act^2 * (F_delta_dot' * Kff_inv * F_delta_dot);
    end
    q_diag_act = reshape([q_diag_act_pos'; q_diag_act_vel'], [], 1) / Kff(1,1);

    % Similarity transform only on actuator states (aero already scaled)
    T_vec = [ones(n_aero, 1); sqrt(q_diag_act)];
    T_inv_vec = 1 ./ T_vec;
    G_iv.A = diag(T_vec) * G_iv.A * diag(T_inv_vec);
    G_iv.B = diag(T_vec) * G_iv.B;
    G_iv.C = G_iv.C * diag(T_inv_vec);

    G(:,:,i_v) = G_iv;
end
