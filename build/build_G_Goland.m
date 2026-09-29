function G = build_G_Goland(V_inf,imuIDX,ailIDX)
% build_G_Goland        Build the state space model array of the Goland
%                       benchmark wing (1 trailing edge flap, 1 IMU at the
%                       midpoint of the flap hinge line), gridded over V_inf.
%
% INPUT
%   V_inf       :       Vector of free stream velocities (m/s); one model per entry
%   imuIDX      :       Index of selected IMU (the Goland model has 1: use 1)
%   ailIDX      :       Index of selected aileron (1 trailing edge flap: use 1)
%
% OUTPUT
%   G(s)        :       [ny, nu, nv] = [1, 1, numel(V_inf)] array of plant
%                       state space models, PT2 actuator included, states and
%                       outputs diagonally scaled (StateScale/OutScale).
%                       Order 2*num_modes + num_modes*num_poles + 2 = 42.
%
%                       u_z1_ddot    <----      flap1_d
%
% Runs the full pipeline (FEM, DLM solve, RFA fit) on every call - there is
% no cached RFA fit for the Goland wing; expect about a second.
%
% Example:  G = build_G_Goland(linspace(20,190,42),1,1);

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

InputNameNoAct = cell(1,3*num_AIL);
InputNameDesired = cell(1,num_AIL);
for i = 1:num_AIL
    InputNameDesired{i} = ['flap',num2str(i),'_d'];

    InputNameNoAct{i} = ['flap',num2str(i)];
    InputNameNoAct{num_AIL+i} = ['flap',num2str(i),'_dot'];
    InputNameNoAct{2*num_AIL+i} = ['flap',num2str(i),'_ddot'];
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

u_z_ddot_scale   = repmat(Kff(1,1)/Mff(1,1), num_IMU,1); 
OutScale = u_z_ddot_scale;

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
    G_noAct = ss(diag(1./StateScale)*A_noS*diag(StateScale), ...
                 diag(1./StateScale)*B_noS, ...
                 diag(1./OutScale)*C_noS*diag(StateScale), ...
                 diag(1./OutScale)*D_noS);

    G_noAct.InputName = InputNameNoAct;
    G_noAct.OutputName = OutputNameNoAct;
    G_noAct.StateName = StateNameNoAct;

    G_noAct = G_noAct(imuIDX,aildddIDX);

    G_iv = connect(G_noAct,G_act{1:end},InputNameDesired,OutputNameDesired);
    G(:,:,i_v) = G_iv;
end