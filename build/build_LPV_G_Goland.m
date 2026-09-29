function LPV_G = build_LPV_G_Goland(V_inf_ref,imuIDX,ailIDX)
% build_LPV_G_Goland    Build the linear parameter varying state space model
%                       (lpvss) of the Goland benchmark wing. LPV parameter
%                       is the free stream velocity V_inf.
%
% INPUT
%   V_inf_ref   :       Reference free stream velocity for scaling.
%                       NOTE: the controller shipped in data/Goland_Cont_nO4_lpf.mat
%                       was synthesized against V_inf_ref = 100 - use 100
%                       whenever that controller closes the loop.
%   imuIDX      :       Index of selected IMU (the Goland model has 1: use 1)
%   ailIDX      :       Index of selected aileron (1 trailing edge flap: use 1)
%
% OUTPUT
%   LPV_G       :       lpvss, u_z1_ddot  <----  flap1_d
%
% Freeze at a velocity:         G90 = psample(LPV_G, [], 90)
% Simulate a velocity ramp:     t = (0:0.002:10)';  V_t = 100 + 8*t;
%                               y = lsim(LPV_G, u, t, [], V_t)
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

num_x_L = num_modes*num_poles;
q_bar_ref = 0.5*Aero.rho*V_inf_ref^2;
pre.StateScale = [ones(num_modes,1);
                  Structure.OMEGA;
                  repmat(q_bar_ref*Aero.c_ref^2,num_x_L,1)];
pre.OutScale = repmat(Structure.Kff(1,1)/Structure.Mff(1,1), num_IMU,1);

% --- LPVSS ---------------------------------------------------------------
DF = @(t,p) dataFcnABCD_G(t,p,pre);
LPV_G = lpvss("V_inf",DF,0,0,V_inf_ref, ...
              'InputName',{'flap1_d'},'OutputName',{'u_z1_ddot'});
end

% --- dataFcn for LPVSS ----------------------------------------------------
function [A,B,C,D,E,dx0,x0,u0,y0,Delay] = dataFcnABCD_G(~,V_inf,pre)
    [A_noS,B_noS,C_noS,D_noS] = build_ABCD_G(V_inf,pre.Structure,pre.Aero);

    sx = pre.StateScale;
    so = pre.OutScale;
    As = A_noS./sx.*(sx');              % diag(1/sx)*A*diag(sx)
    Bs = B_noS./sx;                     % inputs [flap, flap_dot, flap_ddot]
    Cs = C_noS./so.*(sx');
    Ds = D_noS./so;

    % wire in the PT2 actuator: x = [x_plant; x_act], u = flap1_d
    Aa = pre.Aact; Ba = pre.Bact; Ca = pre.Cact; Da = pre.Dact;
    n_p = size(As,1);
    n_a = size(Aa,1);
    A = [As,               Bs*Ca;
         zeros(n_a,n_p),   Aa];
    B = [Bs*Da;
         Ba];
    C = [Cs, Ds*Ca];
    D = Ds*Da;

    E = [];             % Optional descriptor matrix
    % Offsets and delays
    dx0   = [];         % derivative offset
    x0    = [];         % state offset
    u0    = [];         % input offset
    y0    = [];         % output offset
    Delay = [];         % structure with .Input and .Output fields
end
