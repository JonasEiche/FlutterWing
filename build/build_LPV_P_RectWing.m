function LPV_P = build_LPV_P_RectWing(V_inf_ref,imuIDX,ailIDX,modesIDX)
% build_LPV_P_RectWing  Build the linear parameter varying state space model 
%                       (lpvss)of the rectengular wing generalized reference 
%                       model with selected IMUs and Ailerons.
%                       Disturbance is applied on the structural modes,
%                       noise on the IMUs (acceleration sensors). The 
%                       performance output are structural modes deflection,
%                       rate and actuator demand. LPV Parameter is free
%                       stream velocity V_inf: P90 = psample(LPV_P, [], 90)
% 
% INPUT
%   V_inf_ref   :       Reference free stream velocity for Scaling
%   imuIDX      :       Indices of selected IMUs e.g. [4,8]
%   ailIDX      :       Indices of selected Ailerons  [4,8]
%   modesIDX    :       Indices of selected structural modes  [1,2]
% 
%   Flap [1234] / Slat [5678] Aileron & Acceleration Sensor Location:
%
%                        -5-6-7-8-
%                       |         |
%                        -1-2-3-4-
%
% OUTPUT
%   G(s)        :       [ny, nu, nv] Array of plant state space models
% 
%                       q_f1                dist_q_f1_dd
%                       q_f2                dist_q_f2_dd
%                       q_f1_dot            noise_u_z4_d
%                       q_f2_dot    <----   noise_u_z8_d
%                       flap4_d     <----   flap4_d
%                       slat4_d             slat4_d
%                       u_z4_ddot
%                       u_z8_ddot
%                       
% V_inf_ref = 100;
% imuIDX = [4,8];
% ailIDX = [4,8];
% modesIDX = [1,2];
% 
% V_inf=90;
% P90 = psample(LPV_P, [], V_inf)
% 
% --- Model Order Properties ----------------------------------------------
num_modes = 5;
num_poles = 6;
% num_x_L         = num_modes*num_poles;

% --- HARD CODED ----------------------------------------------------------
num_AIL = 8;
num_IMU = 8;
% --- HARD CODED ----------------------------------------------------------

% --- Structure, Aero, Actuator PT2 Model ---------------------------------
[Structure, Aero] = define_RectWing_Structure_Aero(num_modes,num_poles);
[G_act] = define_RectWing_PT2Actuator(num_AIL);

% --- Name Definition -----------------------------------------------------
assert(num_AIL==8,'Hard Coded 4 Slats & 4 Flaps');
InputNameDesired = cell(1,num_modes + num_IMU + num_AIL);
for i = 1:num_modes
    InputNameDesired{i} = ['dist_q_f', num2str(i), '_ddot'];
end
for i = 1:num_IMU
    InputNameDesired{num_modes + i} = ['noise_u_z', num2str(i), '_ddot'];
end
for i = 1:4
    idx = num_modes + num_IMU + i;
    InputNameDesired{idx} = ['flap',num2str(i),'_d'];
    InputNameDesired{4+idx} = ['slat',num2str(i),'_d'];
end

OutputNameDesired = cell(1, 2 * num_modes + num_IMU + num_AIL);
for i = 1:num_modes
    OutputNameDesired{i} = ['q_f', num2str(i)];
    OutputNameDesired{num_modes + i} = ['q_f', num2str(i), '_dot'];
end
for i = 1:4
    idx = 2 * num_modes + i;
    OutputNameDesired{idx} = ['flap',num2str(i),'_d'];
    OutputNameDesired{4+idx} = ['slat',num2str(i),'_d'];
end
for i = 1:num_IMU
    idx = 2 * num_modes + num_AIL + i;
    OutputNameDesired{idx} = ['u_z',num2str(i),'_ddot'];
end
y_q_fIDX            = 1:num_modes;
y_q_f_dotIDX        = num_modes + (1:num_modes);
y_phi_dIDX          = 2*num_modes + (1:num_AIL);
y_u_z_ddotIDX       = 2*num_modes + num_AIL + (1:num_IMU);

u_dist_q_f_ddotIDX  = 1:num_modes;
u_noise_u_z_ddotIDX = num_modes + (1:num_IMU);
u_phi_dIDX          = num_modes + num_IMU + (1:num_AIL);

InputNameDesired = InputNameDesired([u_dist_q_f_ddotIDX(modesIDX), u_noise_u_z_ddotIDX(imuIDX), u_phi_dIDX(ailIDX)]);
OutputNameDesired = OutputNameDesired([y_q_fIDX(modesIDX),y_q_f_dotIDX(modesIDX),y_phi_dIDX(ailIDX),y_u_z_ddotIDX(imuIDX)]);

DF = @(t,p) dataFcnABCD_P(t,p,imuIDX,ailIDX,modesIDX,Structure,Aero, G_act, V_inf_ref);
LPV_P   = lpvss("V_inf",DF,0,0,V_inf_ref,'InputName',InputNameDesired, ...
                            'OutputName',OutputNameDesired);
end

function [A,B,C,D,E,dx0,x0,u0,y0,Delay] = dataFcnABCD_P(~,V_inf,imuIDX,ailIDX,modesIDX,Structure,Aero, G_act, V_inf_ref)
    % --- HARD CODED ------------------------------------------------------
    num_AIL = 8;
    num_IMU = 8;
    % --- HARD CODED ------------------------------------------------------    
    num_modes = size(Structure.Kff,1);
    num_poles = length(Aero.poles);
    num_x_L = num_modes*num_poles;
  
    % --- Name Definition -----------------------------------------------------
    assert(num_AIL==8,'Hard Coded 4 Slats & 4 Flaps');
    InputNameDesired = cell(1,num_modes + num_IMU + num_AIL);
    for i = 1:num_modes
        InputNameDesired{i} = ['dist_q_f', num2str(i), '_ddot'];
    end
    for i = 1:num_IMU
        InputNameDesired{num_modes + i} = ['noise_u_z', num2str(i), '_ddot'];
    end
    for i = 1:4
        idx = num_modes + num_IMU + i;
        InputNameDesired{idx} = ['flap',num2str(i),'_d'];
        InputNameDesired{4+idx} = ['slat',num2str(i),'_d'];
    end
    OutputNameDesired = cell(1, 2 * num_modes + num_IMU + num_AIL);
    for i = 1:num_modes
        OutputNameDesired{i} = ['q_f', num2str(i)];
        OutputNameDesired{num_modes + i} = ['q_f', num2str(i), '_dot'];
    end
    for i = 1:4
        idx = 2 * num_modes + i;
        OutputNameDesired{idx} = ['flap',num2str(i),'_d'];
        OutputNameDesired{4+idx} = ['slat',num2str(i),'_d'];
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
    for i = 1:4
        idx = num_modes + num_IMU + i;
    
        InputNameNoAct{idx} = ['flap',num2str(i)];
        InputNameNoAct{num_AIL+idx} = ['flap',num2str(i),'_dot'];
        InputNameNoAct{2*num_AIL+idx} = ['flap',num2str(i),'_ddot'];
        
        InputNameNoAct{4+idx} = ['slat',num2str(i)];
        InputNameNoAct{num_AIL+4+idx} = ['slat',num2str(i),'_dot'];
        InputNameNoAct{2*num_AIL+4+idx} = ['slat',num2str(i),'_ddot'];
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
    
    % StateNameNoAct = cell(1,2*num_modes+num_x_L);
    % for i = 1:num_modes
    %     StateNameNoAct{i} = ['q_f',num2str(i)];
    %     StateNameNoAct{num_modes+i} = ['q_f',num2str(i),'_dot'];
    % end
    % for i = 1:num_x_L
    %     StateNameNoAct{2*num_modes+i} = ['aero_lag',num2str(i)];
    % end
    
    % --- Scale Definition ----------------------------------------------------
    rho = Aero.rho;
    c_ref = Aero.c_ref;
    OMEGA = Structure.OMEGA;
    q_bar_ref           = 0.5*rho*V_inf_ref^2;
    Kff = Structure.Kff;
    Mff = Structure.Mff;
    
    sqrt_gamma_sc = diag( sqrt(Kff)/sqrt(Kff(1,1)) ); % k
    sqrt_mue_sc = diag( sqrt(Mff)/sqrt(Kff(1,1)) );   % m
    
    StateScale = [ones(num_modes,1);
                  OMEGA;
                  repmat(q_bar_ref*c_ref^2,num_x_L,1)];
    
    
    dist_q_f_ddot_scale = 1./sqrt_mue_sc;
    noise_u_z_ddot_scale = repmat((sqrt_gamma_sc(1)./sqrt_mue_sc(1)) .* 1./sqrt_mue_sc(1), num_IMU,1); 
    phi_scale = ones(num_AIL,1);
    phi_dot_scale = ones(num_AIL,1);
    phi_ddot_scale = ones(num_AIL,1);
    
    InScale = [ dist_q_f_ddot_scale;
                noise_u_z_ddot_scale;
                phi_scale;
                phi_dot_scale;
                phi_ddot_scale ];
    
    q_f_scale = 1./sqrt_gamma_sc;
    q_f_dot_scale = 1./sqrt_mue_sc;
    u_z_ddot_scale   = repmat((sqrt_gamma_sc(1)./sqrt_mue_sc(1)) .* 1./sqrt_mue_sc(1), num_IMU,1); 
    
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

    [A_noS,B_noS,C_noS,D_noS] = build_ABCD_P(V_inf,Structure,Aero);
    P_noAct = ss(diag(1./StateScale)*A_noS*diag(StateScale),diag(1./StateScale)*B_noS*diag(InScale),diag(1./OutScale)*C_noS*diag(StateScale),diag(1./OutScale)*D_noS*diag(InScale));
    P_noAct.InputName = InputNameNoAct;
    P_noAct.OutputName = OutputNameNoAct;
    % P_noAct.StateName = StateNameNoAct;

    P_noAct = P_noAct(yIDX,uIDX);

    P = connect(P_noAct,G_act{1:end},InputNameDesired,OutputNameDesired);
    [A,B,C,D] = ssdata(P);
    E = [];             % Optional descriptor matrix
    % Offsets and delays
    dx0   = [];         % derivative offset
    x0    = [];         % state offset
    u0    = [];         % input offset
    y0    = [];         % output offset
    Delay = [];         % structure with .Input and .Output fields
end