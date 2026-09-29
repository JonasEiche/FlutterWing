%% Structural Blending for Active Flutter Suppression — Controller Synthesis
%
%  Principal research script for IFASD 2026 paper:
%    "Structural Blending for Active Flutter Suppression"
%    Paper ID: IFASD-2026-0068
%
%  PAPER SECTIONS: 2 (Benchmark Model), 3 (Structural Blending),
%                  4 (Controller Synthesis), 5.1--5.4 (Results data)
%
%  OUTPUTS (saved to Data/):
%    controller_imu18_ail18_structural_blending.mat
%      cont_SB       8x8 ss   Structural blending controller
%      cont_H2       8x8 ss   $\mathcal{H}_2$-optimal blending controller (baseline)
%      nOcont_SB     2x2 ss   Inner MIMO controller (structural blending)
%      nOcont_H2     1x1 ss   Inner SISO controller ($\mathcal{H}_2$-optimal)
%      fSoft_SB      1x3      Soft goal values (structural blending)
%      gHard_SB      scalar   Hard constraint value (structural blending)
%      fSoft_H2      1x3      Soft goal values ($\mathcal{H}_2$-optimal)
%      gHard_H2      scalar   Hard constraint value ($\mathcal{H}_2$-optimal)
%
%    Read by R2_Sensor_Fault_Tolerance.m, Fig_Vg_Vpz_OLvsH2vsSB.m,
%    Fig_SigmaPlot_H2vsSB.m and Fig_FreqSwepTimesim_H2vsSB.m. The shipped
%    file holds 8x8 controllers (imuIDX = ailIDX = 1:8).
%
% =========================================================================
%  STRUCTURAL BLENDING METHODOLOGY
% =========================================================================
%
%  --- Output Blending $K_Y$ (from structural mode shapes) ---
%
%  The output mode shape matrix
%    $\Phi_{zf} = \Phi_{zg} \, \Phi_{gf} \in \mathbb{R}^{n_y \times n_f}$
%  maps structural eigenmodes to IMU sensor outputs. Column $i$ contains
%  the $i$-th wind-off eigenmode evaluated at the $n_y = 8$ sensor
%  locations. The measurement equation is:
%    $a_z = \Phi_{zf} \, \ddot{q}_f + \text{noise}$
%
%  Output blending extracts the two flutter modal coordinates from the
%  8-sensor array via least-squares inversion of the flutter submatrix:
%    $K_Y = [\Phi_{zf}(:,1{:}n_c)^+]^\top
%           \in \mathbb{R}^{n_y \times n_c}$
%
%  The blended output $y_\mathrm{blend} = K_Y^\top a_z$ is a least-squares
%  estimate of the flutter modal accelerations
%  $[\ddot{q}_{f1};\, \ddot{q}_{f2}]$. With 8 sensors and 2 unknowns,
%  the overdetermined system provides noise rejection through spatial
%  averaging. The pseudoinverse yields the minimum-norm blending weights.
%
%  --- Input Blending $K_U$ (from generalized aerodynamic forces) ---
%
%  The generalized aerodynamic force matrix
%    $\Phi_{fx} = S_{fj}\!\left(Q_{0,jj}\,D^{Re}_{jx}
%      + \sum_{i=1}^{n_p} Q^{(L_i)}_{jj}\,\hat{D}^{(i)}_{jx}\right)
%      \in \mathbb{R}^{n_f \times n_u}$
%
%  maps actuator deflections to generalized modal forces, where:
%    $S_{fj}$:          modal spline matrix (virtual work projection)
%    $Q_{0,jj}$:        steady-state AIC matrix (quasi-steady aero)
%    $Q^{(L_i)}_{jj}$:  RFA lag coefficients for pole $p_i$
%    $D^{Re}_{jx}$:     downwash from actuator deflection
%    $D^{Im}_{jx}$:     downwash from actuator deflection rate
%    $\hat{D}^{(i)}_{jx} = D^{Re}_{jx} - D^{Im}_{jx}\,p_i\,2/c_\mathrm{ref}$
%
%  Entry $(i,j)$ of $\Phi_{fx}$ quantifies how effectively a unit
%  deflection of actuator $j$ excites structural mode $i$.
%
%  Input blending distributes virtual modal commands across the 8
%  actuators via right pseudoinverse of the flutter-mode rows:
%    $K_U = \Phi_{fx}(1{:}n_c,:)^+
%           \in \mathbb{R}^{n_u \times n_c}$
%
%  The right pseudoinverse yields the minimum-norm actuator deflection
%  pattern that produces a prescribed generalized force on the two
%  flutter modes.
%
%  --- Inner Controller Dimension ---
%
%  The structural blending preserves the two-mode structure: the inner
%  controller $K_\mathrm{int}(s)$ is $2 \times 2$ MIMO ($n_c = 2$),
%  allowing independent targeting of the bending and torsion flutter
%  modes. The $\mathcal{H}_2$-optimal baseline collapses to a single
%  scalar loop: $K_\mathrm{int}(s)$ is $1 \times 1$ SISO ($n_c = 1$).
%
%
% =========================================================================
%  SYNTHESIS FRAMEWORK
% =========================================================================
%  Both controllers are synthesized via structured $\mathcal{H}_\infty$
%  optimization (\texttt{systune}) within identical weighting and constraint
%  frameworks. The multi-model design spans 8 velocities simultaneously.
%
%  --- Weighting Functions ---
%
%  Theis bandpass weighting [Theis2020] confines actuator activity to
%  the flutter frequency band:
%    $W_u(s) = \frac{(s+\omega_1)(s+\omega_2)}{(s+0.01\omega_1)(0.01s+\omega_2)}$
%  with $\omega_1 = 12$ rad/s, $\omega_2 = 64$ rad/s. This passband
%  covers the flutter frequency ($\approx 28$ rad/s) with margin.
%  Performance scalings: $V_p = 0.2$, $V_n = 0.1$, $V_d = 0.5$,
%  $V_u = 0.5$.
%
%  --- Tuning Goals ---
%
%  Hard constraint (must be satisfied):
%    Stability margins: gain $\geq 6$ dB, phase $\geq 45°$ at the
%    plant input (analysis point 'ud').
%
%  Soft constraints (minimized):
%    1. Modal displacement attenuation from disturbance:
%       $\|T_{w \to q_f}\|_\infty \leq V_p / V_d$
%    2. Control effort bound for sensor noise:
%       $\|T_{n \to u}\|_\infty \leq W_u^{-1} \cdot V_u / V_n$
%    3. Control effort bound for disturbances:
%       $\|T_{w \to u}\|_\infty \leq W_u^{-1} \cdot V_u / V_d$
%
%  --- Controller Structure ---
%
%  Structural blending:
%    $K = K_U \, K_\mathrm{int}(s) \, K_Y^\top$ where
%    $K_\mathrm{int}(s) \in \mathbb{R}^{n_c \times n_c}$ ($2 \times 2$)
%    is a tunable state-space controller of order $n_O = 4$.
%  $\mathcal{H}_2$ baseline:
%    $K = k_u \, K_\mathrm{int}(s) \, k_y^\top$ where
%    $K_\mathrm{int}(s) \in \mathbb{R}^{1 \times 1}$ (SISO)
%    is a tunable state-space controller of order $n_O = 4$.
%
%  Inner controller order $n_O = 4$ for both.
%  \texttt{systune} optimization: 3 random starts, parallel evaluation.
%  Random seed: 20260331 for reproducibility.
%
% =========================================================================
%  SCRIPT STRUCTURE
% =========================================================================
%  Section 1: Setup — model dimensions, velocity sweep, sensor/actuator
%             indices, number of flutter modes to target.
%  Section 2: Build generalized plant $P(s)$ for synthesis.
%  Section 3: Define weighting functions ($W_u$) and tuning goals
%             (margins, attenuation, control effort bounds).
%  Section 4: Compute structural blending matrices $K_Y$, $K_U$ from
%             wind-off mode shapes and generalized aerodynamic forces.
%  Section 5: Compute $\mathcal{H}_2$-optimal blending vectors via
%             eigendecomposition and real parametric decomposition.
%  Section 6: Controller synthesis — $\mathcal{H}_2$ blending (\texttt{systune},
%             $1 \times 1$ inner controller).
%  Section 7: Controller synthesis — structural blending (\texttt{systune},
%             $2 \times 2$ inner controller).
%  Section 8: Assign I/O names and save controllers to Data/.
%
% =========================================================================
%  DEPENDENCIES
% =========================================================================
%  define_RectWing_Structure_Aero()   structural/aero parameters
%  build_P_RectWing()                 generalized plant construction
%  build_G_RectWing()                 bare input-output plant
%
%  RELATED SCRIPTS:
%  R1_geometric_mode_isolation_structural.m  — extends SB with exact
%      residual mode rejection (Sec. 3.5)
%  R2_Sensor_Fault_Tolerance.m  — evaluates fault tolerance under
%      progressive sensor failures (Sec. 5.5)
%  Fig_Vg_Vpz_OLvsH2vsSB.m     — V-g flutter stability plots (Sec. 5.1)
%  Fig_SigmaPlot_H2vsSB.m       — singular value analysis (Sec. 5.2--5.3)
%  Fig_FreqSwepTimesim_H2vsSB.m — time-domain validation (Sec. 5.4)
%
% =========================================================================
%  NOTATION
% =========================================================================
%  Symbol            Code variable   Dim        Description
%  -----------------------------------------------------------------------
%  $\Phi_{zf}$       PHIzf           8 x 5      Sensor-to-mode map
%  $\Phi_{fx}$       PHIfx           5 x 8      Mode-to-actuator force map
%  $K_Y$             ky_SB           8 x 2      Output blending (SB)
%  $K_U$             ku_SB           8 x 2      Input blending (SB)
%  $k_y$             ky_H2           8 x 1      Output blending ($\mathcal{H}_2$)
%  $k_u$             ku_H2           8 x 1      Input blending ($\mathcal{H}_2$)
%  $K_\mathrm{int}$  nOcont_SB       2x2 ss     Inner controller (SB)
%  $K_\mathrm{int}$  nOcont_H2       1x1 ss     Inner controller ($\mathcal{H}_2$)
%  $n_y = 8$         nym                        Number of IMU sensors
%  $n_u = 8$         nud                        Number of actuators
%  $n_f = 5$         num_modes                  Retained structural modes
%  $n_c = 2$         nm                         Virtual channels (flutter modes)
%  $n_O = 4$         nO                         Inner controller order
%  $V_\mathrm{ref}$  V_inf_ref       90 m/s     Reference velocity
%  $S_{fj}$          Sfj             5 x 200    Modal spline matrix
%  $Q_{0,jj}$        Q0jj            200 x 200  Steady-state AIC
%  $Q^{(L_i)}_{jj}$  QLpjj           200x200x6  RFA lag coefficients
%  $D^{Re}_{jx}$     DRe_jx          200 x 8    actuator deflection downwash
%  $D^{Im}_{jx}$     DIm_jx          200 x 8    actuator deflection rate downwash
%
clearvars
num_modes = 5;
num_poles = 6;
V_inf_ref = 90;%130;          % design (synthesis) velocity passed as V_inf to build_G_RectWing - not the StateScale reference V_inf_ref = 100 hard-coded inside the builders
V_inf = [90,100,110,120,130,140,150,160];
ailIDX = 1:8;                 % AIL used
imuIDX = 1:8;               % IMU (acc) used
modesIDX = [1,2];
nm = length(modesIDX);

P = build_P_RectWing(V_inf,imuIDX,ailIDX,modesIDX);
nym = length(imuIDX);
nud = length(ailIDX);

w1 = 12;
w2 = 64;
s = tf('s');
W_theis=((s+w1)*(s+w2))/((s+0.01*w1)*(0.01*s+w2));
invW_theis = ((s+0.01*w1)*(0.01*s+w2))/((s+w1)*(s+w2));

Vp=0.2; % 1/8;
Vn=0.1;
Vd=0.5;
Vu=0.5; % 1
invWu = invW_theis;

nO=4;
ReqMarg = TuningGoal.Margins('ud',6,45);
% ReqMarg = TuningGoal.Margins('ud',6,30);
noise_u_z_ddot = {'noise_u_z1_ddot','noise_u_z2_ddot','noise_u_z3_ddot','noise_u_z4_ddot','noise_u_z5_ddot','noise_u_z6_ddot','noise_u_z7_ddot','noise_u_z8_ddot'};
flap_d_slat_d  = {'flap1_d','flap2_d','flap3_d','flap4_d','slat1_d','slat2_d','slat3_d','slat4_d'};
ReqContEffNoise = TuningGoal.Gain(noise_u_z_ddot,flap_d_slat_d,invWu*Vu*(1/Vn));
ReqContEffDist = TuningGoal.Gain({'dist_q_f1_ddot','dist_q_f2_ddot'},flap_d_slat_d,invWu*Vu*(1/Vd));
ReqAttenQF = TuningGoal.Gain({'dist_q_f1_ddot','dist_q_f2_ddot'},{'q_f1','q_f2'},Vp*(1/Vd));

rng(20260331,'twister'); 
% Open a local pool if the Parallel Computing Toolbox is installed; without it
% systune ignores 'UseParallel' and evaluates the random starts serially (slower).
if ~isempty(ver('parallel')) && isempty(gcp('nocreate')), parpool('Processes'); end
opt = systuneOptions('RandomStart',3,'UseParallel',true);

% ------------------------------------------------------------------------
% Structural Blending Matrices
% ------------------------------------------------------------------------

[Structure, Aero] = define_RectWing_Structure_Aero(num_modes,num_poles);
Sfj = Structure.Sfj;
PHIgf = Structure.PHIgf;
DRe_jx = Structure.DRe_jx;
DIm_jx = Structure.DIm_jx;
PHIzg = Structure.PHIzg;

c_ref = Aero.c_ref;
rho = Aero.rho;
poles = Aero.poles;
Q0jj = Aero.Q0jj;
% Q1jj = Aero.Q1jj;
QLpjj = Aero.QLpjj;

num_panels = size(Q0jj,1);
num_AIL = size(DRe_jx,2);

sumQLpjjBjx = zeros(num_panels,num_AIL);
for i = 1:num_poles
    sumQLpjjBjx = sumQLpjjBjx + QLpjj(:,:,i)*(DRe_jx-DIm_jx*poles(i)*2/c_ref);
end
% Bgx         = 0.5*rho*V_inf_ref^2*Sfj*(Q0jj*DRe_jx+sumQLpjjBjx);
PHIfx         = Sfj*(Q0jj*DRe_jx+sumQLpjjBjx);
PHIxf12 = pinv( PHIfx(1:2,:));

PHIzf = PHIzg * PHIgf;
PHIf12z = pinv( PHIzf(:,1:2) );

ku_SB = PHIxf12;
ky_SB = PHIf12z';

% ------------------------------------------------------------------------
% H_2 Blending Vectors
% ------------------------------------------------------------------------

G = build_G_RectWing(V_inf_ref,imuIDX,ailIDX);
[V,DD] = eig(G.A);
flutIDX = find( (real(diag(DD)) > -1)  & (abs(imag(diag(DD))) < 40) & (abs(imag(diag(DD))) > 5) );
critIDX = find( (real(diag(DD)) > -10)  & (abs(imag(diag(DD))) < 40) & (abs(imag(diag(DD))) > 5) );
resIDX = find( (real(diag(DD)) > -10)  & (abs(imag(diag(DD))) < 120) & (abs(imag(diag(DD))) > 40) );
Am = DD;
Bm = V\G.B;
Cm = G.C*V;
Dm = G.D;

Am_flut = Am(flutIDX,flutIDX);
Bm_flut = Bm(flutIDX,:);
Cm_flut = Cm(:,flutIDX);

Am_crit = Am(critIDX,critIDX);
Bm_crit = Bm(critIDX,:);
Cm_crit = Cm(:,critIDX);

cmG_flut = ss(Am_flut,Bm_flut,Cm_flut,0);

% The code for the calculation of the H_2 optimal blending vectors is proprietary, hence hardcoded here: 
ky_H2 = [-0.0995;-0.3201;-0.5322;-0.6997;0.0764;0.1800;0.2131;0.1764];
ku_H2 = [0.0450;0.3120;0.4477;0.5603;0.1156;0.2680;0.3858;0.3902];

H2_NORM_cmG_flut = norm(cmG_flut,2)
H2_NORM_ky_cmG_flut_ku = norm(ky_H2'*cmG_flut*ku_H2,2)

% ------------------------------------------------------------------------
% Controller Synthesis H2 Blending
% ------------------------------------------------------------------------

tuneCont = ku_H2*tunableSS('nOcont',nO,1,1)*ky_H2';
tuneCL = lft(P,AnalysisPoint('ud',nud)*tuneCont*AnalysisPoint('ym',nym));
% wu wn & wu wd & wp wd
[CL_H2,fSoft_H2,gHard_H2] = systune(tuneCL,[ReqAttenQF,ReqContEffNoise,ReqContEffDist],ReqMarg,opt);
nOcont_H2 = getBlockValue(CL_H2,'nOcont');
cont_H2 = ku_H2*nOcont_H2*ky_H2';

% ------------------------------------------------------------------------
% Controller Synthesis Structural Blending
% ------------------------------------------------------------------------

tuneCont = ku_SB*tunableSS('nOcont',nO,2,2)*ky_SB';
tuneCL = lft(P,AnalysisPoint('ud',nud)*tuneCont*AnalysisPoint('ym',nym));
% wu wn & wu wd & wp wd
[CL_SB,fSoft_SB,gHard_SB] = systune(tuneCL,[ReqAttenQF,ReqContEffNoise,ReqContEffDist],ReqMarg,opt);
nOcont_SB = getBlockValue(CL_SB,'nOcont');
cont_SB = ku_SB*nOcont_SB*ky_SB';


% ------------------------------------------------------------------------
% Output (Actuator) Names & Input (Sensor) Names
% ------------------------------------------------------------------------

ActuatorName = P.InputName(nm+nym+(1:nud));
SensorName = P.OutputName(2*nm+nud+(1:nym));
cont_H2.InputName = SensorName;
cont_H2.OutputName = ActuatorName;
cont_SB.InputName = SensorName;
cont_SB.OutputName = ActuatorName;

% ------------------------------------------------------------------------
% Save to Data/controller_*.mat
% ------------------------------------------------------------------------

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
data_dir = fullfile(script_dir, 'Data');
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end
ContPath = fullfile(data_dir, 'controller_imu18_ail18_structural_blending.mat');
save(ContPath,"cont_H2","cont_SB","nOcont_H2","nOcont_SB", ...
    "fSoft_H2","gHard_H2","fSoft_SB","gHard_SB");




