function LPV_G = build_LPV_G_RectWing(V_inf_ref,imuIDX,ailIDX)
% build_G_RectWing      Build the linear parameter varying state space model 
%                       (lpvss) of the rectengular wing reference model 
%                       with selected IMUs and Ailerons. LPV Parameter is free
%                       stream velocity V_inf: G90 = psample(LPV_G, [], 90)
% 
% INPUT
%   V_inf_ref   :       Reference free stream velocity for Scaling
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
% 
%                       u_z1_ddot               flap1_d
%                       u_z2_ddot               flap2_d
%                       u_z3_ddot               flap3_d
%                       u_z4_ddot    <----      flap4_d
%                       u_z5_ddot    <----      slat1_d
%                       u_z6_ddot               slat2_d
%                       u_z7_ddot               slat3_d
%                       u_z8_ddot               slat4_d
%                       
% V_inf_ref = 100;
% imuIDX = [4,8];
% ailIDX = [4,8];
% 
% V_inf=90;
% G90 = psample(LPV_G, [], V_inf)
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

% --- Name Definition Desired ---------------------------------------------
InputNameDesired = cell(1,num_AIL);
for i = 1:4
    InputNameDesired{i} = ['flap',num2str(i),'_d'];
    InputNameDesired{4+i} = ['slat',num2str(i),'_d'];
end
OutputNameDesired = cell(1,num_IMU);
for i = 1:num_IMU
    OutputNameDesired{i} = ['u_z',num2str(i),'_ddot'];
end
InputNameDesired = InputNameDesired(ailIDX);
OutputNameDesired = OutputNameDesired(imuIDX);

% --- LPVSS ---------------------------------------------------------------
DF = @(t,p) dataFcnABCD_G(t,p,imuIDX,ailIDX,Structure,Aero, G_act, V_inf_ref);
LPV_G   = lpvss("V_inf",DF,0,0,V_inf_ref,'InputName',InputNameDesired, ...
                            'OutputName',OutputNameDesired);

end

% --- dataFcn for LPVSS ---------------------------------------------------
function [A,B,C,D,E,dx0,x0,u0,y0,Delay] = dataFcnABCD_G(~,V_inf,imuIDX,ailIDX,Structure,Aero, G_act, V_inf_ref)
    % --- HARD CODED ------------------------------------------------------
    num_AIL = 8;
    num_IMU = 8;
    % --- HARD CODED ------------------------------------------------------

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
    % StateNameNoAct = cell(1,2*num_modes+num_x_L);
    % for i = 1:num_modes
    %     StateNameNoAct{i} = ['q_f',num2str(i)];
    %     StateNameNoAct{num_modes+i} = ['q_f',num2str(i),'_dot'];
    % end
    % for i = 1:num_x_L
    %     StateNameNoAct{2*num_modes+i} = ['aero_lag',num2str(i)];
    % end
    
    % --- Scale Definition ------------------------------------------------
    
    rho = Aero.rho;
    c_ref = Aero.c_ref;
    OMEGA = Structure.OMEGA;
    q_bar_ref           = 0.5*rho*V_inf_ref^2;
    Kff = Structure.Kff;
    Mff = Structure.Mff;
    sqrt_gamma_sc = diag( sqrt(Kff)/sqrt(Kff(1,1)) ); % k
    sqrt_mue_sc = diag( sqrt(Mff)/sqrt(Kff(1,1)) );   % m
    num_modes = size(Kff,1);
    num_poles = length(Aero.poles);
    num_x_L         = num_modes*num_poles;
    StateScale = [ones(num_modes,1);
                  OMEGA;
                  repmat(q_bar_ref*c_ref^2,num_x_L,1)];
    
    u_z_ddot_scale   = repmat((sqrt_gamma_sc(1)./sqrt_mue_sc(1)) .* 1./sqrt_mue_sc(1), num_IMU,1); 
    OutScale = u_z_ddot_scale;
    
    % --- Build State Space Model Array -----------------------------------
    aildddIDX = [ailIDX,num_AIL+ailIDX,2*num_AIL+ailIDX];
    InputNameDesired = InputNameDesired(ailIDX);
    OutputNameDesired = OutputNameDesired(imuIDX);
    G_act = G_act(ailIDX);
    [A_noS,B_noS,C_noS,D_noS] = build_ABCD_G(V_inf,Structure,Aero);    
    G_noAct = ss(diag(1./StateScale)*A_noS*diag(StateScale), ...
                 diag(1./StateScale)*B_noS, ...
                 diag(1./OutScale)*C_noS*diag(StateScale), ...
                 diag(1./OutScale)*D_noS);
    % G_noAct = ss(A_noS,B_noS,C_noS,D_noS);
    
    G_noAct.InputName = InputNameNoAct;
    G_noAct.OutputName = OutputNameNoAct;
    % G_noAct.StateName = StateNameNoAct;
    
    G_noAct = G_noAct(imuIDX,aildddIDX);
    
    G = connect(G_noAct,G_act{1:end},InputNameDesired,OutputNameDesired);
    
    [A,B,C,D] = ssdata(G);
    
    E = [];             % Optional descriptor matrix
    % Offsets and delays
    dx0   = [];         % derivative offset
    x0    = [];         % state offset
    u0    = [];         % input offset
    y0    = [];         % output offset
    Delay = [];         % structure with .Input and .Output fields
end