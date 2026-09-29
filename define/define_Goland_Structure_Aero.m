function [Structure, Aero] = define_Goland_Structure_Aero(num_modes,num_poles)
% define_Goland_Structure_Aero  Goland benchmark wing: 1 trailing-edge flap, 1 IMU.
%   [Structure, Aero] = define_Goland_Structure_Aero(num_modes, num_poles)
%
%   Classic validation case (semi span 6.096 m, chord 1.8288 m, flexural axis at
%   33 % chord, mass axis at 43 % chord); flutter near 175.6 m/s in the literature.
%   10 beam elements, 5 x 10 DLM panels (chordwise x spanwise). Runs the pipeline
%   FEM -> modal truncation -> DLM + Roger RFA -> coupling and returns the two
%   dicts consumed by build_ABCD_G / build_ABCD_P (called through build_G_Goland,
%   build_P_Goland, build_LPV_G_Goland, build_LPV_P_Goland). No RFA cache: the
%   DLM solve runs every call (about 1 s).
%
%   num_modes : number of structural modes kept
%   num_poles : number of Roger RFA lag poles
%
%   Structure  (f = num_modes modal DOFs, g = 3*(num_ele+1) FEM DOFs,
%               j = 50 panels, x = 1 flap, z = 1 IMU)
%     Kff, Mff, Dff   generalized stiffness / mass / damping, diagonal,
%                     Dff = Mff*diag(2*dampRatio*OMEGA)
%     OMEGA           eigenfrequencies (rad/s), ascending
%     PHIgf           modeshapes, largest entry of each column = 1 (see build_PHIgf)
%     Sfj, Sgj        panel pressure coefficients -> modal / nodal forces
%     DRe_jf, DIm_jf  modal displacement / velocity -> panel downwash
%     DRe_jg, DIm_jg  same on FEM DOFs
%     DRe_jx, DIm_jx  flap deflection / rate -> panel downwash
%     PHIzg           FEM DOF acceleration -> IMU vertical acceleration
%     E, ele          beam node coordinates and element-to-DOF table
%     Pa, Ps          panel corner points in aero / structural coordinates
%     Mgg, Kgg        physical FEM mass / stiffness matrices
%     cspanels        panel numbers of the flap
%   Aero
%     c_ref, rho      reference chord (m), air density (kg/m^3)
%     poles           Roger RFA poles in k_red, 1 x num_poles
%     Q0jj, Q1jj      quasi-steady and first-order AIC terms, j x j
%     QLpjj           lag terms, j x j x num_poles
%
% Jonas * July 2025
% _________________________________________________________________________

% ---- Dimensions ---------------------------------------------------------
rho = 1.020;                % air density  (kg/m^3)
s = 6.096;                  % semi span (m)
c = 1.8288;                 % root chord (m)
c_ref = c;                  % reference chord (m)
sw = 0;                     % sweep angle at leading edge (deg)
dh = 0;                     % dihedral (deg)
tr = 1;                     % taper ratio

xf = 0.33;                  % flexural axis position relative to chord: xf=0.5 is mid chord flexural axis
dampRatio = 0.00; %0.01;    % structural damping ratio (0 for the undamped benchmark): Dff = Mff*diag(2*dampRatio*OMEGA)
xm = 0.43;                  % mass axis position relative to chord

rho_bar = 35.71;            % mass per unit length (kg/m) : int_A rho dA
I_Tym = 8.64;               % polar mass moment of inertia per unit length (kg*m) : int_A rho r^2 dA
I_zm = rho_bar*-(xm-xf)*c;  % static mass moment per unit length (kg) : int_A rho x dA  (bending-torsion coupling; negative because the mass axis lies behind the flexural axis)
EI_xxa = 9.77e6;            % flexural rigidity (N*m^2)   : E(N/m^2) * int_A z^2 dA
GI_Tya = 0.99e6;            % torsional rigidity (N*m^2)  : G(N/m^2) * int_A r^2 dA


np_s = 10;                  % number of DLM panels spanwise
np_c = 5;                   % number of DLM panels chordwise

num_ele = 10;               % number of FEM elements
% num_modes and num_poles are passed in by the entry points (build_G_Goland etc.)

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
% alt (unused): sparse Mff = spdiags(diag(Mff),0,num_modes,num_modes);
Mff = diag(diag(Mff));
% equivalent: Kff = PHIgf'*Kgg*PHIgf;
Kff = Mff*diag(OMEGA.^2);
Dff = Mff*diag(2*dampRatio*OMEGA);
Sfj = PHIgf'*Sgj;
DRe_jf = DRe_jg*PHIgf;
DIm_jf = DIm_jg*PHIgf;

% ---- Unsteady Aerodynamics ----------------------------------------------
% No cache for the Goland wing: 50 panels (SYM = 1 mirrors them) x 14 reduced
% frequencies solve in about 1 s.
Ma = 0.0;
k_red = [0, 0.001, 0.01, 0.02, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 1.1];
SYM =1;
Qjj = build_Qjj(Ma,k_red,c_ref,Pa,SYM);
[poles,Q0jj,Q1jj,~,QLpjj,~,~,~] = rogersRFA_magW(k_red, Qjj, num_poles);

% ---- Save to Structure --------------------------------------------------
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
% ---- Save to Aero -------------------------------------------------------
Aero.c_ref = c_ref;
Aero.rho = rho;
Aero.poles = poles;
Aero.Q0jj = Q0jj;
Aero.Q1jj = Q1jj;
Aero.QLpjj = QLpjj;

% ---- Define Flaps -------------------------------------------------------
% One trailing-edge flap on the outer 3 of the 10 spanwise panels of the TE row
% (build_PaPs numbers panels rowwise from the trailing edge, root -> tip).
cspanels{1} = [8,9,10];
rot_axis_a{1}{1} = Pa{cspanels{1}(1)}{1};
rot_axis_a{1}{2} = Pa{cspanels{1}(end)}{4};
[DRe_jx,DIm_jx] = build_DReDIm_jx(cspanels,rot_axis_a,Pa);

Structure.DRe_jx = DRe_jx;
Structure.DIm_jx = DIm_jx;
Structure.cspanels = cspanels;

% ---- Define IMUs --------------------------------------------------------
% One IMU measuring z-acceleration at the midpoint of the flap hinge axis.
rot_axis_flap_1_s = Ps{cspanels{1}(1)}{1};
rot_axis_flap_2_s = Ps{cspanels{1}(end)}{4};
HP_flap = 0.5*(rot_axis_flap_1_s+rot_axis_flap_2_s);
x_pos_IMU=HP_flap(1);
y_pos_IMU=HP_flap(2);
[PHIzg, ~] = build_IMU(x_pos_IMU, y_pos_IMU, E, ele);

Structure.PHIzg = PHIzg;