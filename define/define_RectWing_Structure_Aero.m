function [Structure, Aero] = define_RectWing_Structure_Aero(num_modes, num_poles)

% Kff = Structure.Kff
% Mff = Structure.Mff
% Dff = Structure.Dff
% Sfj = Structure.Sfj
% DRe_jf = Structure.DRe_jf
% DIm_jf = Structure.DIm_jf
% DRe_jx = Structure.DRe_jx
% DIm_jx = Structure.DIm_jx
% HM_xj = Structure.HM_xj
% PHIzg = Structure.PHIzg
% OMEGA = Structure.OMEGA
% E = Structure.E
% ele = Structure.ele
% Pa = Structure.Pa
% Ps = Structure.Ps
% Sgj = Structure.Sgj
% DRe_jg = Structure.DRe_jg
% DIm_jg = Structure.DIm_jg
% Mgg = Structure.Mgg;
% Kgg = Structure.Kgg;
% DRe_jx = Structure.DRe_jx;
% DIm_jx = Structure.DIm_jx;
% cspanels = Structure.cspanels;
% 
% c_ref = Aero.c_ref
% rho = Aero.rho
% poles = Aero.poles
% Q0jj = Aero.Q0jj
% Q1jj = Aero.Q1jj
% QLpjj = Aero.QLpjj

% ---- Structure ----------------------------------------------------------

% num_modes = 5;          % Structural modes to keep
% num_poles = 6;          % Number of poles in the Rogers RFA Fit 

% num_IMU = 8;            % HARD CODED
% num_AIL = 8;            % HARD CODED

s = 7.5;                % semi span
c = 2;                  % root chord
c_ref = c;              % reference chord
np_s = 20;              % number panels spanwise
np_c = 10;              % number panels chordwise
sw = 0;                 % sweep angle at leading edge
dh = 0;                 % dihedral
tr = 1;                 % taper ratio
xf = 0.48;              % flex axis position relative to chord
%xf=0.5
dampRatio = 0.01;       % Dff = Mff*diag(2*dampRatio*OMEGA);
num_ele = 16;
rho = 1.225;                % air density  (kg/m^3)

% rho_bar = 200;          % Mass per unit length = int_A rho dA
% I_Tym = 66.9867;       % Torsionsmoment = int_A rho r^2 dA
% I_zm = -8.0;            % Deviationsmoment = int_A rho * x dA   (dynamische Kopplung von Biegung und Torsion)
% EI_xxa = 3.1228e+07;    % Flexural rigidity  : E * int_A z^2 dA
% GI_Tya = 1.7851e+06;    % Torsional rigidity : G * int_A r^2 dA
xfc=xf*c;
rho_wing = 100;                             % mass per unit wing area
rho_bar = rho_wing*c;                       % mass per unit length = int_A rho dA           
I_Tym = rho_wing*1/3*(xfc^3-(xfc-c)^3);    % Torsionsmoment = int_A rho r^2 dA
I_zm = rho_wing*c*(xfc-0.5*c);            % Deviationsmoment = int_A rho * x dA   (dynamische Kopplung von Biegung und Torsion)
bending_freq = 5;   % 5 bending freq in Hz - approximate - ignores coupling term; build_PHIgf_wrightcooper
torsion_freq = 6;   % 10 torsion freq in Hz - approximate - ignores coupling term; real frequency will be approx. 27% lower due to stiffening assumption of qudratic shape function see. Wright, Cooper P.55
EI_xxa = (bending_freq * pi * 2)^2 * (rho_wing * s * c / 5) / 4 * s^3;   % flexural rigidity  : E * int_A z^2 dA
GI_Tya = (torsion_freq * pi * 2)^2 * (rho_wing * s / 3 * (c^3 / 3 - c^2 * xfc + xfc^2 * c)) * s;         % torsional rigidity : G * int_A r^2 dA

% ---- Structural Dynamics & Coupling -------------------------------------
[E, ele] = build_E_y(s, 0, num_ele);
[Pa,Ps] = build_PaPs(s,c,np_s,np_c,sw,dh,tr,xf);
Sgj = build_Sgj(E,ele,Ps);
[DRe_jg, DIm_jg] = build_DReDIm_jg(E,ele,Ps);
Mgg = build_Mgg(E, ele, rho_bar, I_Tym, I_zm);
Kgg = build_Kgg(E, ele, EI_xxa, GI_Tya);

% ---- Project FEM to Modeshapes q_f = PHIgf'*q_g -------------------------
[PHIgf,OMEGA] = build_PHIgf(Kgg,Mgg,num_modes);
Mff = PHIgf'*Mgg*PHIgf;
max_offdiag = max(max(abs(Mff - diag(diag(Mff))))); % test diagonality
assert(max_offdiag < 1e-8,"Generalized Mass Matrix Mff is not diagonal")
% Mff = spdiags(diag(Mff),0,num_modes,num_modes);
Mff = diag(diag(Mff));
% Kff = PHIgf'*Kgg*PHIgf;
Kff = Mff*diag(OMEGA.^2);
Dff = Mff*diag(2*dampRatio*OMEGA);
Sfj = PHIgf'*Sgj;
DRe_jf = DRe_jg*PHIgf;
DIm_jf = DIm_jg*PHIgf;
% Keep FEM DoFs and apply Dirichlet BC on left beam root (fixed)
%     Mff = Mgg(4:end,4:end);
%     Kff = Kgg(4:end,4:end);
%     Dgg = zeros(size(Kgg));
%     Dff = Dgg(4:end,4:end);
%     Sfj = Sgj(4:end,:);
%     DRe_jf = DRe_jg(:,4:end);
%     DIm_jf = DIm_jg(:,4:end);
%     num_modes = size(Kff,1);
%     PHIgf = [zeros(3,num_modes);
%              eye(num_modes)];

% ---- Safe to Structure --------------------------------------------------

Structure.Kff = Kff;
Structure.Mff = Mff;
Structure.Dff = Dff;
Structure.Sfj = Sfj;
Structure.DRe_jf = DRe_jf;
Structure.DIm_jf = DIm_jf;
Structure.PHIgf = PHIgf;
Structure.OMEGA = OMEGA;
Structure.E = E;
Structure.ele = ele;
Structure.Pa = Pa;
Structure.Ps = Ps;
Structure.Sgj = Sgj;
Structure.DRe_jg = DRe_jg;
Structure.DIm_jg = DIm_jg;
Structure.Mgg = Mgg;
Structure.Kgg = Kgg;


% ---- Unsteady Aerodynamics ----------------------------------------------
if num_poles == 6
    script_path = mfilename('fullpath');
    script_dir = fileparts(script_path);
    load(fullfile(script_dir,'Q0jjQ1jjQLpjj_k_red11_p6.mat'), 'poles', 'Q0jj', 'Q1jj', 'QLpjj')
else
    Ma = 0.0;
    k_red = [0, 0.001, 0.01, 0.02, 0.05, 0.07, 0.1, 0.2, 0.5, 0.7, 0.9, 1.1];
    SYM =1;
    Qjj = build_Qjj(Ma,k_red,c_ref,Pa,SYM);
    [poles,Q0jj,Q1jj,~,QLpjj,~,~,~] = rogersRFA_magW(k_red, Qjj, num_poles);
end
% ---- Safe to Aero -------------------------------------------------------

Aero.c_ref = c_ref;
Aero.rho = rho;
Aero.poles = poles;
Aero.Q0jj = Q0jj;
Aero.Q1jj = Q1jj;
Aero.QLpjj = QLpjj;

% ---- Control Surface Definition -----------------------------------------

%    Slats
%  -5-6-7-8-
% |         |
%  -1-2-3-4-
%    Flaps

% Hard coded Panel numbers assume np_s = 20 (spanwise) & np_c = 10 (chordwise) panels

% ---- Define Flaps -------------------------------------------------------
cspanels{1} = [1,2,3,4,5];
rot_axis_a{1}{1} = Pa{cspanels{1}(1)}{1};
rot_axis_a{1}{2} = Pa{cspanels{1}(end)}{4};

cspanels{2} = [6,7,8,9,10];
rot_axis_a{2}{1} = Pa{cspanels{2}(1)}{1};
rot_axis_a{2}{2} = Pa{cspanels{2}(end)}{4};

cspanels{3} = [11,12,13,14,15];
rot_axis_a{3}{1} = Pa{cspanels{3}(1)}{1};
rot_axis_a{3}{2} = Pa{cspanels{3}(end)}{4};

cspanels{4} = [16,17,18,19,20];
rot_axis_a{4}{1} = Pa{cspanels{4}(1)}{1};
rot_axis_a{4}{2} = Pa{cspanels{4}(end)}{4};

% ---- Define Slats -------------------------------------------------------
cspanels{5} = 180+[1,2,3,4,5];
rot_axis_a{5}{1} = Pa{cspanels{5}(end)}{3};
rot_axis_a{5}{2} = Pa{cspanels{5}(1)}{2};

cspanels{6} = 180+[6,7,8,9,10];
rot_axis_a{6}{1} = Pa{cspanels{6}(end)}{3};
rot_axis_a{6}{2} = Pa{cspanels{6}(1)}{2};

cspanels{7} = 180+[11,12,13,14,15];
rot_axis_a{7}{1} = Pa{cspanels{7}(end)}{3};
rot_axis_a{7}{2} = Pa{cspanels{7}(1)}{2};

cspanels{8} = 180+[16,17,18,19,20];
rot_axis_a{8}{1} = Pa{cspanels{8}(end)}{3};
rot_axis_a{8}{2} = Pa{cspanels{8}(1)}{2};

% ---- DRe_jx,DIm_jx,HM_xj ------------------------------------------------
[DRe_jx,DIm_jx] = build_DReDIm_jx(cspanels,rot_axis_a,Pa);

Structure.DRe_jx = DRe_jx;
Structure.DIm_jx = DIm_jx;
Structure.cspanels = cspanels;

% ---- IMU Acceleration Collocated Sensor Definition ----------------------

%    Slats
%  -5-6-7-8-
% |         |
%  -1-2-3-4-
%    Flaps

% Measure z-dir. acceleration in the center of each control surface hinge axis

num_AIL=length(cspanels);
num_IMU=num_AIL;
x_pos_IMU=zeros(num_IMU,1);
y_pos_IMU=zeros(num_IMU,1);
for i_cs = 1:4
    % Flaps
    rot_axis_flap_1_s = Ps{cspanels{i_cs}(1)}{1};
    rot_axis_flap_2_s = Ps{cspanels{i_cs}(end)}{4};
    
    HP_flap = 0.5*(rot_axis_flap_1_s+rot_axis_flap_2_s);
    x_pos_IMU(i_cs)=HP_flap(1);
    y_pos_IMU(i_cs)=HP_flap(2);
    % Slats
    rot_axis_slat_1_s = Ps{cspanels{4+i_cs}(end)}{3};
    rot_axis_slat_2_s = Ps{cspanels{4+i_cs}(1)}{2};

    HP_slat = 0.5*(rot_axis_slat_1_s+rot_axis_slat_2_s);
    x_pos_IMU(4+i_cs)=HP_slat(1);
    y_pos_IMU(4+i_cs)=HP_slat(2);
end
[PHIzg, ~] = build_IMU(x_pos_IMU, y_pos_IMU, E, ele);

Structure.PHIzg = PHIzg;
