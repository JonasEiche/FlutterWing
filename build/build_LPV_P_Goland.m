function LPV_P = build_LPV_P_Goland(V_inf_ref,imuIDX,ailIDX,modesIDX)
% build_LPV_P_Goland    Build the linear parameter varying state space model
%                       (lpvss) of the Goland generalized plant. Disturbance
%                       is applied on the structural modes, noise on the IMU
%                       (acceleration sensor). The performance outputs are
%                       structural mode deflection, rate and actuator demand.
%                       LPV parameter is the free stream velocity V_inf.
%
% INPUT
%   V_inf_ref   :       Reference free stream velocity for scaling.
%                       NOTE: the controller shipped in data/Goland_Cont_nO4_lpf.mat
%                       was synthesized against V_inf_ref = 100 - use 100
%                       whenever that controller closes the loop.
%   imuIDX      :       Index of selected IMU (the Goland model has 1: use 1)
%   ailIDX      :       Index of selected aileron (1 trailing edge flap: use 1)
%   modesIDX    :       Indices of selected structural modes, e.g. [1,2]
%
% OUTPUT
%   LPV_P       :       lpvss
%
%                       q_f1                 dist_q_f1_ddot
%                       q_f2                 dist_q_f2_ddot
%                       q_f1_dot    <----    noise_u_z1_ddot
%                       q_f2_dot    <----    flap1_d
%                       flap1_d
%                       u_z1_ddot
%
% Freeze at a velocity:         P90 = psample(LPV_P, [], 90)
% Simulate a velocity ramp:     t = (0:0.002:10)';  V_t = 100 + 8*t;
%                               [y,~,x] = lsim(LPV_P, u, t, [], V_t)
% Close the loop (u = K*y):     CL = lft(LPV_P, K)
%                               (states: plant first, then controller)
%
% The data function assembles the plant + PT2 actuator with explicit block
% matrices instead of connect() so that a single evaluation costs well under
% a millisecond - lsim on an lpvss evaluates the data function roughly twice
% per time step (TR-BDF2), i.e. tens of thousands of times per simulation.
%
% --- Model Order Properties ----------------------------------------------
num_modes = 5;
num_poles = 6;

% --- HARD CODED ----------------------------------------------------------
num_AIL = 1;
num_IMU = 1;
% --- HARD CODED ----------------------------------------------------------
assert(isequal(imuIDX,1) && isequal(ailIDX,1), ...
    'Goland model has exactly 1 IMU and 1 flap: imuIDX = ailIDX = 1')

% --- Structure, Aero, Actuator PT2 Model (built ONCE - DLM solve here) ----
[Structure, Aero] = define_Goland_Structure_Aero(num_modes,num_poles);
G_act = define_PT2Actuator(num_AIL);

% --- everything V-independent, precomputed for the data function ----------
pre.Structure = Structure;
pre.Aero = Aero;
[pre.Aact,pre.Bact,pre.Cact,pre.Dact] = ssdata(G_act{1});
pre.modesIDX = modesIDX(:)';

num_x_L = num_modes*num_poles;
q_bar_ref = 0.5*Aero.rho*V_inf_ref^2;
Kff = Structure.Kff;
Mff = Structure.Mff;
pre.StateScale = [ones(num_modes,1);
                  Structure.OMEGA;
                  repmat(q_bar_ref*Aero.c_ref^2,num_x_L,1)];

dist_q_f_ddot_scale  = 1./(sqrt(diag(Mff))/sqrt(Kff(1,1)));
noise_u_z_ddot_scale = repmat(Kff(1,1)/Mff(1,1), num_IMU,1);
pre.InScale = [ dist_q_f_ddot_scale;
                noise_u_z_ddot_scale;
                ones(3*num_AIL,1) ];        % phi, phi_dot, phi_ddot

q_f_scale     = 1./(sqrt(diag(Kff))/sqrt(Kff(1,1)));
q_f_dot_scale = 1./(sqrt(diag(Mff))/sqrt(Kff(1,1)));
u_z_ddot_scale = repmat(Kff(1,1)/Mff(1,1), num_IMU,1);
pre.OutScale = [ q_f_scale;
                 q_f_dot_scale;
                 u_z_ddot_scale ];

% --- Name Definition ------------------------------------------------------
nm = numel(modesIDX);
InputNameDesired = cell(1,nm+2);
OutputNameDesired = cell(1,2*nm+2);
for i = 1:nm
    InputNameDesired{i} = ['dist_q_f',num2str(modesIDX(i)),'_ddot'];
    OutputNameDesired{i} = ['q_f',num2str(modesIDX(i))];
    OutputNameDesired{nm+i} = ['q_f',num2str(modesIDX(i)),'_dot'];
end
InputNameDesired{nm+1} = 'noise_u_z1_ddot';
InputNameDesired{nm+2} = 'flap1_d';
OutputNameDesired{2*nm+1} = 'flap1_d';
OutputNameDesired{2*nm+2} = 'u_z1_ddot';

% --- LPVSS ---------------------------------------------------------------
DF = @(t,p) dataFcnABCD_P(t,p,pre);
LPV_P = lpvss("V_inf",DF,0,0,V_inf_ref, ...
              'InputName',InputNameDesired,'OutputName',OutputNameDesired);
end

% --- dataFcn for LPVSS ----------------------------------------------------
function [A,B,C,D,E,dx0,x0,u0,y0,Delay] = dataFcnABCD_P(~,V_inf,pre)
    [A_noS,B_noS,C_noS,D_noS] = build_ABCD_P(V_inf,pre.Structure,pre.Aero);

    sx = pre.StateScale;
    si = pre.InScale;
    so = pre.OutScale;
    As = A_noS./sx.*(sx');              % diag(1/sx)*A*diag(sx)
    Bs = B_noS./sx.*(si');
    Cs = C_noS./so.*(sx');
    Ds = D_noS./so.*(si');

    % unscaled build_ABCD_P channel layout (num_modes = 5, 1 IMU, 1 flap):
    % u = [dist_q_f_ddot(1:5), noise(6), flap(7), flap_dot(8), flap_ddot(9)]
    % y = [q_f(1:5), q_f_dot(6:10), u_z_ddot(11)]
    num_modes = 5;
    mIDX = pre.modesIDX;
    nm = numel(mIDX);
    uD = mIDX;  uN = num_modes+1;  uF = num_modes+1+(1:3);
    yQ = mIDX;  yQd = num_modes+mIDX;  yU = 2*num_modes+1;

    % wire in the PT2 actuator: x = [x_plant; x_act]
    % u = [dist(modesIDX), noise, flap1_d]
    % y = [q_f(modesIDX), q_f_dot(modesIDX), flap1_d, u_z1_ddot]
    Aa = pre.Aact; Ba = pre.Bact; Ca = pre.Cact; Da = pre.Dact;
    n_p = size(As,1);
    n_a = size(Aa,1);
    A = [As,               Bs(:,uF)*Ca;
         zeros(n_a,n_p),   Aa];
    B = [Bs(:,uD),         Bs(:,uN),        Bs(:,uF)*Da;
         zeros(n_a,nm),    zeros(n_a,1),    Ba];
    C = [Cs(yQ,:),         Ds(yQ,uF)*Ca;
         Cs(yQd,:),        Ds(yQd,uF)*Ca;
         zeros(1,n_p),     zeros(1,n_a);    % flap1_d output = demand passthrough
         Cs(yU,:),         Ds(yU,uF)*Ca];
    D = [Ds(yQ,uD),        Ds(yQ,uN),       Ds(yQ,uF)*Da;
         Ds(yQd,uD),       Ds(yQd,uN),      Ds(yQd,uF)*Da;
         zeros(1,nm),      0,               1;
         Ds(yU,uD),        Ds(yU,uN),       Ds(yU,uF)*Da];

    E = [];             % Optional descriptor matrix
    % Offsets and delays
    dx0   = [];         % derivative offset
    x0    = [];         % state offset
    u0    = [];         % input offset
    y0    = [];         % output offset
    Delay = [];         % structure with .Input and .Output fields
end
