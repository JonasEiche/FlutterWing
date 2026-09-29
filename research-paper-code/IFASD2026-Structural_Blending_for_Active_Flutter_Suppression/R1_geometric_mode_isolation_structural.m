%% Geometric Mode Isolation via Structural Blending
%
%  Research script for IFASD 2026 paper:
%    "Structural Blending for Active Flutter Suppression"
%
%  PAPER SECTION: 3.5 Geometric Mode Isolation
%
%  FIGURES: Fig_Coupling_Pinv_vs_GeomIso.m
%    - Figures/Fig16_Output_Coupling_Pinv_vs_GeomIso.pdf
%    - Figures/Fig17_Input_Coupling_Pinv_vs_GeomIso.pdf
%
%  OUTPUTS (saved to Data/):
%    - R1_geometric_mode_isolation_results.mat   (blending vectors, coupling
%      matrices; read by Fig_Coupling_Pinv_vs_GeomIso.m)
%    - controller_imu18_ail18_geometric_isolation.mat  (synthesized
%      controllers: cont_SB, cont_ISO, inner controllers, blending vectors;
%      kept for inspection, read by no other script)
%
% =========================================================================
%  MOTIVATION
% =========================================================================
%  The structural blending approach (R0) computes output and input blending
%  matrices from the pseudoinverse of the flutter-mode submatrix alone:
%
%    $K_Y = [\Phi_{zf}(:,1{:}n_c)^+]^\top
%           \in \mathbb{R}^{n_y \times n_c}$     (output, $8 \times 2$)
%    $K_U = \Phi_{fx}(1{:}n_c,:)^+
%           \in \mathbb{R}^{n_u \times n_c}$      (input, $8 \times 2$)
%
%  These guarantee flutter mode recovery
%  ($K_Y^\top C_\mathrm{flut} = I_{n_c}$ and
%  $\Phi_{fx}(1{:}n_c,:) K_U = I_{n_c}$)
%  but do NOT guarantee rejection of the residual modes 3--5.
%  The pseudoinverse of the submatrix is agnostic to the existence of
%  other modes, so spillover occurs:
%
%    $K_Y^{\mathrm{SB},\top} \Phi_{zf} = [I_{n_c} \mid \neq 0]$
%     (residual columns nonzero)
%    $\Phi_{fx} K_U^\mathrm{SB} = [I_{n_c};\, \neq 0]$
%     (residual rows nonzero)
%
%  This script investigates geometric mode isolation, which additionally
%  enforces exact residual rejection:
%
%    $K_Y^{\mathrm{iso},\top} \Phi_{zf} = [I_{n_c} \mid 0_{n_c \times (n_f - n_c)}]$
%     (zero coupling to modes 3--5)
%    $\Phi_{fx} K_U^\mathrm{iso} = [I_{n_c};\, 0_{(n_f - n_c) \times n_c}]$
%     (zero excitation of modes 3--5)
%
% =========================================================================
%  MATHEMATICAL FRAMEWORK
% =========================================================================
%  Two equivalent formulations are implemented and verified:
%
%  --- Method B: Full pseudoinverse (computationally simplest) ---
%
%  When $\Phi_{zf}$ ($n_y \times n_f = 8 \times 5$) has full column rank,
%  $\Phi_{zf}^+$ is an $n_f \times n_y$ left inverse satisfying
%  $\Phi_{zf}^+ \Phi_{zf} = I_{n_f}$. Taking the first $n_c$ rows yields
%  an $n_c \times n_y$ matrix that recovers modes 1 and 2 with identically
%  zero coupling to modes 3--5:
%
%    $K_Y^{\mathrm{iso},\top} = \Phi_{zf}^+(1{:}n_c,\,:)$
%     ($n_c \times n_y = 2 \times 8$)
%    $K_U^\mathrm{iso} = (\Phi_{fx}^+)(:,\,1{:}n_c)$
%     ($n_u \times n_c = 8 \times 2$; dual, columns of right pseudoinverse)
%
%  This is the minimum-norm solution to the combined constraint
%  $K_Y^\top \Phi_{zf} = [I_{n_c},\, 0_{n_c \times (n_f - n_c)}]$.
%
%  --- Method C: Orthogonal projection (theoretical derivation) ---
%
%  Partition $\Phi_{zf} = [C_\mathrm{flut},\, C_\mathrm{res}]$ where
%  $C_\mathrm{flut} = \Phi_{zf}(:,1{:}n_c) \in \mathbb{R}^{n_y \times n_c}$
%  and $C_\mathrm{res} = \Phi_{zf}(:,n_c{+}1{:}n_f) \in \mathbb{R}^{n_y \times (n_f - n_c)}$.
%  The isolation operator $P_\mathrm{iso}$ is constructed by:
%
%    1. Projector onto residual subspace:
%       $P_\mathrm{res} = C_\mathrm{res}
%        (C_\mathrm{res}^\top C_\mathrm{res})^{-1}
%        C_\mathrm{res}^\top$
%
%    2. Orthogonal complement (annihilates $C_\mathrm{res}$ by construction):
%       $P^\perp = I_{n_y} - P_\mathrm{res}$
%
%    3. Isolation operator (left pseudoinverse of projected flutter map):
%       $P_\mathrm{iso} = (C_\mathrm{flut}^\top P^\perp C_\mathrm{flut})^{-1}
%        C_\mathrm{flut}^\top P^\perp$
%
%  The simplification from the standard left pseudoinverse
%  $(P^\perp C_\mathrm{flut})^+$ uses the idempotence and symmetry of
%  $P^\perp$: $(P^\perp)^\top P^\perp = P^\perp$.
%
% =========================================================================
%  WELL-POSEDNESS
% =========================================================================
%  The construction requires:
%    - $\Phi_{zf}$ has full column rank ($\mathrm{rank} = n_f$), ensured
%      when $n_y \geq n_f$ ($8 \geq 5$)
%    - $C_\mathrm{res}$ has full column rank ($\mathrm{rank} = n_f - n_c$),
%      so $P_\mathrm{res}$ is well-defined
%    - $C_\mathrm{flut}$ and $C_\mathrm{res}$ span distinct subspaces in
%      measurement space, so $P^\perp C_\mathrm{flut}$ retains rank $n_c$
%
%  A single sufficient condition for all three:
%    $\sigma_\mathrm{min}(\Phi_{zf}) \gg 0$
%
%  The condition number analysis (Fig_Condition_Number_Blending.m)
%  confirms $\kappa(\Phi_{zf}) = 1.24$ (all $n_f = 5$ modes), indicating
%  excellent numerical conditioning. For $\Phi_{fx}$:
%  $\kappa(\Phi_{fx}) = 6.85$ (all $n_f = 5$ modes).
%
%
%  CONTROLLER SYNTHESIS (\texttt{systune}, identical framework to R0):
%
%    Standard SB:           fSoft = [1.05, 1.05, 1.05],  gHard = 0.999
%    Geometric isolation:   fSoft = [1.09, 1.13, 1.13],  gHard = 1.000
%
%
% =========================================================================
%  SCRIPT STRUCTURE
% =========================================================================
%  Section 1: Setup, load structural/aero parameters, compute
%             $\Phi_{zf}$, $\Phi_{fx}$
%  Section 2: Output blending comparison (Methods A/B/C, coupling verify)
%  Section 3: Input blending comparison (Methods A/B, coupling verify)
%  Section 4: Save analysis results ->
%             Data/R1_geometric_mode_isolation_results.mat
%  Section 5: Controller synthesis (\texttt{systune}, identical to R0
%             framework) ->
%             Data/controller_imu18_ail18_geometric_isolation.mat
%
% =========================================================================
%  DEPENDENCIES
% =========================================================================
%  define_RectWing_Structure_Aero()   structural/aero parameters
%  build_P_RectWing()                 generalized plant construction
%
%  RELATED SCRIPTS:
%  R0_afs_structural_blending_imu18_ail18_synthesis.m  — baseline synthesis
%  Fig_Condition_Number_Blending.m  — well-posedness validation
%  Fig_Coupling_Pinv_vs_GeomIso.m  — coupling comparison figures
%
% =========================================================================
%  NOTATION
% =========================================================================
%  Symbol                     Code variable   Dim        Description
%  -----------------------------------------------------------------------
%  $\Phi_{zf}$                PHIzf           8 x 5      Sensor-to-mode map
%  $\Phi_{fx}$                PHIfx           5 x 8      Mode-to-actuator force map
%  $C_\mathrm{flut}$          C_flut          8 x 2      Flutter columns of $\Phi_{zf}$
%  $C_\mathrm{res}$           C_res           8 x 3      Residual columns of $\Phi_{zf}$
%  $P_\mathrm{res}$           P_res           8 x 8      Projector onto residual subspace
%  $P^\perp$                  P_perp          8 x 8      Orthogonal complement projector
%  $P_\mathrm{iso}$           P_iso           2 x 8      Isolation operator
%  $K_Y^\mathrm{SB}$          ky_SB           8 x 2      Output blending (standard SB)
%  $K_U^\mathrm{SB}$          ku_SB           8 x 2      Input blending (standard SB)
%  $K_Y^\mathrm{iso}$         ky_ISO          8 x 2      Output blending (isolation)
%  $K_U^\mathrm{iso}$         ku_ISO          8 x 2      Input blending (isolation)
%  $K_\mathrm{int}$           nOcont_SB       2x2 ss     Inner controller (standard SB)
%  $K_\mathrm{int}$           nOcont_ISO      2x2 ss     Inner controller (isolation)
%  $n_y = 8$                  nym                        Number of IMU sensors
%  $n_u = 8$                  nud                        Number of actuators
%  $n_f = 5$                  num_modes                  Retained structural modes
%  $n_c = 2$                  nm                         Virtual channels (flutter modes)
%  $\kappa(\cdot)$            cond()                     Condition number
%  $(\cdot)^+$                pinv()                     Moore--Penrose pseudoinverse
%  $\sigma_\mathrm{min}$      —                          Smallest singular value

clearvars

script_path = mfilename('fullpath');
script_dir  = fileparts(script_path);
data_dir = fullfile(script_dir, 'Data');
if ~exist(data_dir, 'dir'), mkdir(data_dir); end

% ---- Model parameters ----------------------------------------------------
num_modes = 5;
num_poles = 6;
imuIDX   = 1:8;
ailIDX   = 1:8;
modesIDX = [1,2];
nm = length(modesIDX);

% ---- Load structural/aero parameters ------------------------------------
[Structure, Aero] = define_RectWing_Structure_Aero(num_modes, num_poles);
Sfj    = Structure.Sfj;
PHIgf  = Structure.PHIgf;
DRe_jx = Structure.DRe_jx;
DIm_jx = Structure.DIm_jx;
PHIzg  = Structure.PHIzg;

c_ref  = Aero.c_ref;
poles  = Aero.poles;
Q0jj   = Aero.Q0jj;
QLpjj  = Aero.QLpjj;

num_panels = size(Q0jj,1);
num_AIL    = size(DRe_jx,2);

% ---- Compute blending base matrices -------------------------------------
PHIzf = PHIzg * PHIgf;           % 8 x 5 (IMU -> modes)

sumQLpjjBjx = zeros(num_panels, num_AIL);
for i = 1:num_poles
    sumQLpjjBjx = sumQLpjjBjx + QLpjj(:,:,i)*(DRe_jx - DIm_jx*poles(i)*2/c_ref);
end
PHIfx = Sfj*(Q0jj*DRe_jx + sumQLpjjBjx);  % 5 x 8 (modes <- actuators)

%% ======================================================================
%  OUTPUT BLENDING: Simple pinv vs Geometric Isolation
%  ======================================================================
C_flut = PHIzf(:,1:2);   % 8 x 2
C_res  = PHIzf(:,3:5);   % 8 x 3

% Method A: Simple pseudo-inverse (current SB approach from R0)
ky_pinv = pinv(C_flut);  % 2 x 8

% Method B: Geometric isolation via full pseudo-inverse
ky_iso_full = pinv(PHIzf);     % 5 x 8
ky_iso = ky_iso_full(1:2,:);   % 2 x 8

% Method C: Orthogonal projection (explicit derivation per notes)
P_res  = C_res * ((C_res' * C_res) \ C_res');
P_perp = eye(size(PHIzf,1)) - P_res;
P_iso  = (C_flut' * P_perp * C_flut) \ (C_flut' * P_perp);  % 2 x 8

% ---- Verify coupling matrices -------------------------------------------
coupling_pinv = ky_pinv * PHIzf;   % 2 x 5
coupling_iso  = ky_iso  * PHIzf;   % 2 x 5
coupling_proj = P_iso   * PHIzf;   % 2 x 5

fprintf('=== Output Blending: ky * PHIzf (2x5) ===\n');
fprintf('\n  Simple pinv:\n');
disp(coupling_pinv)
fprintf('  Geometric isolation (full pinv):\n');
disp(coupling_iso)
fprintf('  Orthogonal projection:\n');
disp(coupling_proj)

% ---- Check equivalence of Methods B and C --------------------------------
diff_iso_proj = max(abs(ky_iso - P_iso), [], 'all');
fprintf('Max |ky_iso - P_iso| = %.2e\n\n', diff_iso_proj);

%% ======================================================================
%  INPUT BLENDING: Simple pinv vs Geometric Isolation
%  ======================================================================
B_flut = PHIfx(1:2,:);   % 2 x 8
B_res  = PHIfx(3:5,:);   % 3 x 8

% Method A: Simple pseudo-inverse
ku_pinv = pinv(B_flut);          % 8 x 2

% Method B: Geometric isolation via full pseudo-inverse
ku_iso = pinv(PHIfx);            % 8 x 5
ku_iso = ku_iso(:,1:2);          % 8 x 2

% ---- Verify coupling matrices -------------------------------------------
coupling_ku_pinv = PHIfx * ku_pinv;  % 5 x 2
coupling_ku_iso  = PHIfx * ku_iso;   % 5 x 2

fprintf('=== Input Blending: PHIfx * ku (5x2) ===\n');
fprintf('\n  Simple pinv:\n');
disp(coupling_ku_pinv)
fprintf('  Geometric isolation (full pinv):\n');
disp(coupling_ku_iso)

%% ======================================================================
%  SAVE ANALYSIS RESULTS
%  ======================================================================
ResultsPath = fullfile(data_dir, 'R1_geometric_mode_isolation_results.mat');
save(ResultsPath, ...
    'PHIzf', 'PHIfx', ...
    'ky_pinv', 'ky_iso', 'P_iso', 'ku_pinv', 'ku_iso', ...
    'coupling_pinv', 'coupling_iso', 'coupling_proj', ...
    'coupling_ku_pinv', 'coupling_ku_iso', ...
    'diff_iso_proj');
fprintf('Analysis results saved to %s\n', ResultsPath);

%% ======================================================================
%  CONTROLLER SYNTHESIS: Simple SB vs Geometric Isolation
%  ======================================================================
fprintf('\n=== Controller Synthesis ===\n');

V_inf = [90,100,110,120,130,140,150,160];
P = build_P_RectWing(V_inf, imuIDX, ailIDX, modesIDX);
nym = length(imuIDX);
nud = length(ailIDX);

% Weighting functions (identical to R0)
s = tf('s');
w1 = 12; w2 = 64;
invW_theis = ((s+0.01*w1)*(0.01*s+w2))/((s+w1)*(s+w2));

Vp = 0.2; Vn = 0.1; Vd = 0.5; Vu = 0.5;
invWu = invW_theis;
nO = 4;

% Tuning goals (identical to R0)
ReqMarg = TuningGoal.Margins('ud', 6, 45);
noise_u_z_ddot = {'noise_u_z1_ddot','noise_u_z2_ddot','noise_u_z3_ddot','noise_u_z4_ddot', ...
                  'noise_u_z5_ddot','noise_u_z6_ddot','noise_u_z7_ddot','noise_u_z8_ddot'};
flap_d_slat_d  = {'flap1_d','flap2_d','flap3_d','flap4_d', ...
                  'slat1_d','slat2_d','slat3_d','slat4_d'};
ReqContEffNoise = TuningGoal.Gain(noise_u_z_ddot, flap_d_slat_d, invWu*Vu*(1/Vn));
ReqContEffDist  = TuningGoal.Gain({'dist_q_f1_ddot','dist_q_f2_ddot'}, flap_d_slat_d, invWu*Vu*(1/Vd));
ReqAttenQF      = TuningGoal.Gain({'dist_q_f1_ddot','dist_q_f2_ddot'}, {'q_f1','q_f2'}, Vp*(1/Vd));

rng(20260331,'twister');
% Open a local pool if the Parallel Computing Toolbox is installed; without it
% systune ignores 'UseParallel' and evaluates the random starts serially (slower).
if ~isempty(ver('parallel')) && isempty(gcp('nocreate')), parpool('Processes'); end
opt = systuneOptions('RandomStart', 3, 'UseParallel', true);

% ---- Blending vectors (R0 convention: ky is 8x2) ------------------------
ky_SB  = ky_pinv';   % 8 x 2
ku_SB  = ku_pinv;     % 8 x 2
ky_ISO = ky_iso';     % 8 x 2
ku_ISO = ku_iso;      % 8 x 2

% ---- Synthesis: Simple SB ------------------------------------------------
fprintf('Synthesizing Simple SB controller...\n');
tuneCont = ku_SB * tunableSS('nOcont', nO, nm, nm) * ky_SB';
tuneCL = lft(P, AnalysisPoint('ud', nud) * tuneCont * AnalysisPoint('ym', nym));
[CL_SB, fSoft_SB, gHard_SB] = systune(tuneCL, ...
    [ReqAttenQF, ReqContEffNoise, ReqContEffDist], ReqMarg, opt);
nOcont_SB = getBlockValue(CL_SB, 'nOcont');
cont_SB = ku_SB * nOcont_SB * ky_SB';

% ---- Synthesis: Geometric Isolation --------------------------------------
fprintf('Synthesizing Geometric Isolation controller...\n');
tuneCont = ku_ISO * tunableSS('nOcont', nO, nm, nm) * ky_ISO';
tuneCL = lft(P, AnalysisPoint('ud', nud) * tuneCont * AnalysisPoint('ym', nym));
[CL_ISO, fSoft_ISO, gHard_ISO] = systune(tuneCL, ...
    [ReqAttenQF, ReqContEffNoise, ReqContEffDist], ReqMarg, opt);
nOcont_ISO = getBlockValue(CL_ISO, 'nOcont');
cont_ISO = ku_ISO * nOcont_ISO * ky_ISO';

% ---- Name controllers ----------------------------------------------------
ActuatorName = P.InputName(nm+nym+(1:nud));
SensorName   = P.OutputName(2*nm+nud+(1:nym));
cont_SB.InputName   = SensorName;
cont_SB.OutputName  = ActuatorName;
cont_ISO.InputName  = SensorName;
cont_ISO.OutputName = ActuatorName;

% ---- Print comparison ----------------------------------------------------
fprintf('\n=== Synthesis Results ===\n');
fprintf('Simple SB:            fSoft = [%.4f, %.4f, %.4f],  gHard = %.4f\n', fSoft_SB, gHard_SB);
fprintf('Geometric Isolation:  fSoft = [%.4f, %.4f, %.4f],  gHard = %.4f\n', fSoft_ISO, gHard_ISO);

% ---- Save controllers ----------------------------------------------------
ContPath = fullfile(data_dir, 'controller_imu18_ail18_geometric_isolation.mat');
save(ContPath, "cont_SB", "cont_ISO", "nOcont_SB", "nOcont_ISO", ...
               "ky_SB", "ku_SB", "ky_ISO", "ku_ISO");
fprintf('Controllers saved to %s\n', ContPath);
