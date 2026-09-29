function P = build_P_Goland(V_inf,imuIDX,ailIDX,modesIDX)
% build_P_Goland        Build the generalized plant state space model array of
%                       the Goland benchmark wing (1 trailing edge flap, 1 IMU)
%                       for controller synthesis: build_G_Goland plus
%                       disturbance on the structural modes, noise on the IMU
%                       (acceleration sensor), and the performance outputs
%                       modal deflection, modal rate and actuator demand.
%
% INPUT
%   V_inf       :       Vector of free stream velocities (m/s); one model per entry
%   imuIDX      :       Index of selected IMU (the Goland model has 1: use 1)
%   ailIDX      :       Index of selected aileron (1 trailing edge flap: use 1)
%   modesIDX    :       Indices of selected structural modes, e.g. [1,2]
%
% OUTPUT
%   P(s)        :       [ny, nu, nv] array of generalized plant state space
%                       models, nv = numel(V_inf), PT2 actuator included,
%                       states, inputs and outputs diagonally scaled
%                       (StateScale/InScale/OutScale). For modesIDX = [1,2]:
%
%                       q_f1                 dist_q_f1_ddot
%                       q_f2                 dist_q_f2_ddot
%                       q_f1_dot    <----    noise_u_z1_ddot
%                       q_f2_dot    <----    flap1_d
%                       flap1_d
%                       u_z1_ddot
%
%                       i.e. inputs [w; u] = [dist_q_f_ddot; noise_u_z1_ddot; flap1_d]
%                       and outputs [z; y] = [q_f; q_f_dot; flap1_d; u_z1_ddot],
%                       the partition lft(P,K) expects.
%
% Runs the full pipeline (FEM, DLM solve, RFA fit) on every call - there is
% no cached RFA fit for the Goland wing; expect about a second.
%
% Example:  P = build_P_Goland([90,110,130,150,170,190],1,1,[1,2]);
%           (the synthesis grid of R00_Synthesize_Goland_Controller)

% --- Model Order Properties ----------------------------------------------
num_modes = 5;
num_poles = 6;
num_x_L         = num_modes*num_poles;

% --- HARD CODED ----------------------------------------------------------
V_inf_ref = 100;    % (m/s) reference velocity used ONLY for the aero-lag state
                    % scaling q_bar_ref*c_ref^2 in StateScale below. A state
                    % scaling is a similarity transform: it changes neither the
                    % physics nor the I/O behaviour, an output-feedback
                    % controller is unaffected, and it need not match the
                    % V_inf grid. Keep 100 so that state trajectories match
                    % the shipped results (the LPV builders take V_inf_ref as
                    % an argument instead).
num_AIL = 1;        % 1 trailing edge flap, fixed by define_Goland_Structure_Aero
num_IMU = 1;        % 1 IMU at the flap hinge midpoint, fixed there as well
% --- HARD CODED ----------------------------------------------------------

[Structure, Aero] = define_Goland_Structure_Aero(num_modes,num_poles);
G_act = define_PT2Actuator(num_AIL);


% --- Name Definition -----------------------------------------------------
InputNameDesired = cell(1,num_modes + num_IMU + num_AIL);
for i = 1:num_modes
    InputNameDesired{i} = ['dist_q_f', num2str(i), '_ddot'];
end
for i = 1:num_IMU
    InputNameDesired{num_modes + i} = ['noise_u_z', num2str(i), '_ddot'];
end
for i = 1:num_AIL
    idx = num_modes + num_IMU + i;
    InputNameDesired{idx} = ['flap',num2str(i),'_d'];
end

OutputNameDesired = cell(1, 2 * num_modes + num_IMU + num_AIL);
for i = 1:num_modes
    OutputNameDesired{i} = ['q_f', num2str(i)];
    OutputNameDesired{num_modes + i} = ['q_f', num2str(i), '_dot'];
end
for i = 1:num_AIL
    idx = 2 * num_modes + i;
    OutputNameDesired{idx} = ['flap',num2str(i),'_d'];
end
for i = 1:num_IMU
    idx = 2 * num_modes + num_AIL + i;
    OutputNameDesired{idx} = ['u_z',num2str(i),'_ddot'];
end

InputNameNoAct = cell(1,num_modes + num_IMU + 3*num_AIL);
for i = 1:num_modes
    InputNameNoAct{i} = ['dist_q_f', num2str(i), '_ddot'];
end
for i = 1:num_IMU
    InputNameNoAct{num_modes + i} = ['noise_u_z', num2str(i), '_ddot'];
end
for i = 1:num_AIL
    idx = num_modes + num_IMU + i;
    InputNameNoAct{idx} = ['flap',num2str(i)];
    InputNameNoAct{num_AIL+idx} = ['flap',num2str(i),'_dot'];
    InputNameNoAct{2*num_AIL+idx} = ['flap',num2str(i),'_ddot'];
end

OutputNameNoAct = cell(1, 2 * num_modes + num_IMU);
for i = 1:num_modes
    OutputNameNoAct{i} = ['q_f', num2str(i)];
    OutputNameNoAct{num_modes + i} = ['q_f', num2str(i), '_dot'];
end
for i = 1:num_IMU
    idx = 2 * num_modes + i;
    OutputNameNoAct{idx} = ['u_z',num2str(i),'_ddot'];
end

StateNameNoAct = cell(1,2*num_modes+num_x_L);
for i = 1:num_modes
    StateNameNoAct{i} = ['q_f',num2str(i)];
    StateNameNoAct{num_modes+i} = ['q_f',num2str(i),'_dot'];
end
for i = 1:num_x_L
    StateNameNoAct{2*num_modes+i} = ['aero_lag',num2str(i)];
end

% --- Scale Definition ----------------------------------------------------
rho = Aero.rho;
c_ref = Aero.c_ref;
OMEGA = Structure.OMEGA;
q_bar_ref           = 0.5*rho*V_inf_ref^2;
Kff = Structure.Kff;
Mff = Structure.Mff;

StateScale = [ones(num_modes,1);
              OMEGA;
              repmat(q_bar_ref*c_ref^2,num_x_L,1)];

dist_q_f_ddot_scale = 1./(sqrt(diag(Mff))/sqrt(Kff(1,1)));
noise_u_z_ddot_scale = repmat(Kff(1,1)/Mff(1,1), num_IMU,1);
phi_scale = ones(num_AIL,1);
phi_dot_scale = ones(num_AIL,1);
phi_ddot_scale = ones(num_AIL,1);

InScale = [ dist_q_f_ddot_scale;
            noise_u_z_ddot_scale;
            phi_scale;
            phi_dot_scale;
            phi_ddot_scale ];


q_f_scale = 1./(sqrt(diag(Kff))/sqrt(Kff(1,1)));
q_f_dot_scale = 1./(sqrt(diag(Mff))/sqrt(Kff(1,1)));
u_z_ddot_scale   = repmat(Kff(1,1)/Mff(1,1), num_IMU,1); 

OutScale = [ q_f_scale;
             q_f_dot_scale;
             u_z_ddot_scale ];

% --- Build State Space Model Array ---------------------------------------
y_q_fIDX            = 1:num_modes;
y_q_f_dotIDX        = num_modes + (1:num_modes);
y_phi_dIDX          = 2*num_modes + (1:num_AIL);
y_u_z_ddotIDX       = 2*num_modes + num_AIL + (1:num_IMU);

u_dist_q_f_ddotIDX  = 1:num_modes;
u_noise_u_z_ddotIDX = num_modes + (1:num_IMU);
u_phi_dIDX          = num_modes + num_IMU + (1:num_AIL);

u_phiIDX            = num_modes + num_IMU + (1:3*num_AIL);
w_u_z_ddotIDX       = 2*num_modes + (1:num_IMU);
aildddIDX           = [ailIDX,num_AIL+ailIDX,2*num_AIL+ailIDX];

G_act = G_act(ailIDX);

InputNameDesired = InputNameDesired([u_dist_q_f_ddotIDX(modesIDX), u_noise_u_z_ddotIDX(imuIDX), u_phi_dIDX(ailIDX)]);
OutputNameDesired = OutputNameDesired([y_q_fIDX(modesIDX),y_q_f_dotIDX(modesIDX),y_phi_dIDX(ailIDX),y_u_z_ddotIDX(imuIDX)]);

uIDX = [u_dist_q_f_ddotIDX(modesIDX), u_noise_u_z_ddotIDX(imuIDX), u_phiIDX(aildddIDX)];
yIDX = [y_q_fIDX(modesIDX),y_q_f_dotIDX(modesIDX),w_u_z_ddotIDX(imuIDX)];

n_mode = length(modesIDX);
ny_imu = length(imuIDX);       % u_z_ddot                 (Measured output y)
nu_ail = length(ailIDX);       % flap, slat               (Control input u)
nz = n_mode+n_mode+nu_ail;     % q_f + q_f_dot + phi_d    (Performance output z)
nw = n_mode+ny_imu;            % dist + noise             (Disturbance input w)
nu = nw+nu_ail;
ny = nz+ny_imu;
nv = length(V_inf);
P = ss(zeros(ny,nu,nv));
for i_v = 1:length(V_inf)
    V_inf_iv = V_inf(i_v);
    [A_noS,B_noS,C_noS,D_noS] = build_ABCD_P(V_inf_iv,Structure,Aero);
    P_noAct = ss(diag(1./StateScale)*A_noS*diag(StateScale),diag(1./StateScale)*B_noS*diag(InScale),diag(1./OutScale)*C_noS*diag(StateScale),diag(1./OutScale)*D_noS*diag(InScale));
    P_noAct.InputName = InputNameNoAct;
    P_noAct.OutputName = OutputNameNoAct;
    P_noAct.StateName = StateNameNoAct;

    P_noAct = P_noAct(yIDX,uIDX);

    P_iv = connect(P_noAct,G_act{1:end},InputNameDesired,OutputNameDesired);
    P(:,:,i_v) = P_iv;
end