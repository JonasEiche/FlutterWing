%% R01: Synthesize the RectWing quick-start controller (static 8x8 gain)
%
% Regenerates data/RectWing_Cont_imu18_ail18.mat: RectWing_Cont, the active
% flutter suppression controller that QUICKSTART.m section 5 closes the loop
% with and that the README hero animation (docs/figures/make_readme_figures.m)
% shows. All 8 IMUs feed all 8 control surfaces (4 trailing-edge flaps, 4
% leading-edge slats) through one static output-feedback gain, the simplest
% controller structure there is: u = RectWing_Cont * y, closed with
% feedback(G, -RectWing_Cont) (the lft(P,K) sign convention).
%
% Fixed-structure synthesis with systune on the generalized plant
% build_P_RectWing(90:10:160, 1:8, 1:8, [1,2]), modeled on
% R00_Synthesize_Goland_Controller.m 
% 
% Tuning goals (soft, minimized):
%   1. ReqAttenQF      : gain dist_q_f{1,2}_ddot -> q_f{1,2}
%                        (suppress the flutter-mode response to modal forces)
%   2. ReqContEffNoise : weighted gain noise_u_z{1..8}_ddot -> the 8 demands
%                        (limit actuator wear from sensor noise)
%   3. ReqContEffDist  : weighted gain dist_q_f{1,2}_ddot -> the 8 demands
%                        (limit actuator activity during disturbance rejection)
% Hard constraint:
%   4. ReqMarg         : disk margins 6 dB / 45 deg at the plant input 'ud'
% The Theis weight concentrates the allowed activity in the flutter band:
% w1 = 12, w2 = 64 rad/s around the 28 rad/s (4.5 Hz) flutter mode.
%
% Result: the margin goal is NOT met,
% gHard = 1.47. This tuning result does not establish infeasibility for all
% static gains. The returned gain's achieved multiloop
% disk margin at 'ud', worst case over the synthesis grid, is 0.58: gain
% margin 5.1 dB, phase margin 32 deg.

clearvars

% --- synthesis knobs ------------------------------------------------------
RS   = 0;                           % <<< systune random starts (0 = default start only)
SEED = 20260903;                    % <<< rng seed for reproducibility
rng(SEED,'twister');

num_modes = 5;
num_poles = 6;
imuIDX = 1:8;                       % all 8 IMUs
ailIDX = 1:8;                       % all 8 surfaces: flaps 1-4, slats 1-4
modesIDX = [1,2];                   % target the first two structural modes
V_inf = 90:10:160;                  % LPV synthesis grid (m/s), flutter onset 104.3 m/s

P = build_P_RectWing(V_inf,imuIDX,ailIDX,modesIDX);
nym = length(imuIDX);
nud = length(ailIDX);

% actuator effort weight (Theis bandpass-like filter) centered on the flutter band
w1 = 12;
w2 = 64;
s = tf('s');
invW_theis = ((s+0.01*w1)*(0.01*s+w2))/((s+w1)*(s+w2));

Vp=0.2;
Vn=0.1;
Vd=0.5;
Vu=0.5;
invWu = invW_theis;

noise_u_z_ddot = {'noise_u_z1_ddot','noise_u_z2_ddot','noise_u_z3_ddot','noise_u_z4_ddot', ...
                  'noise_u_z5_ddot','noise_u_z6_ddot','noise_u_z7_ddot','noise_u_z8_ddot'};
flap_d_slat_d  = {'flap1_d','flap2_d','flap3_d','flap4_d','slat1_d','slat2_d','slat3_d','slat4_d'};

ReqMarg = TuningGoal.Margins('ud',6,45);
ReqContEffNoise = TuningGoal.Gain(noise_u_z_ddot,flap_d_slat_d,invWu*Vu*(1/Vn));
ReqContEffDist = TuningGoal.Gain({'dist_q_f1_ddot','dist_q_f2_ddot'},flap_d_slat_d,invWu*Vu*(1/Vd));
ReqAttenQF = TuningGoal.Gain({'dist_q_f1_ddot','dist_q_f2_ddot'},{'q_f1','q_f2'},Vp*(1/Vd));

% controller: one static 8x8 output-feedback gain (no dynamics)
tuneCont = tunableGain('Kafs',nud,nym);
tuneCL = lft(P,AnalysisPoint('ud',nud)*tuneCont*AnalysisPoint('ym',nym));

% Open a local pool only when random starts are actually requested and the
% Parallel Computing Toolbox is installed; without it systune ignores
% 'UseParallel' and evaluates the random starts serially (slower).
UseParallel = RS > 0 && ~isempty(ver('parallel'));
if UseParallel && isempty(gcp('nocreate')), parpool('Processes'); end
opt = systuneOptions('RandomStart',RS,'UseParallel',UseParallel,'Display','final');

% --- systune: the margin goal as the HARD constraint -----------------------
rng(SEED,'twister');
[CL_tuned,fSoft,gHard] = systune(tuneCL,[ReqAttenQF,ReqContEffNoise,ReqContEffDist],ReqMarg,opt);
disp(['Hard goal ReqMarg (must be <= 1) : ',num2str(gHard)])
disp(['Soft goals (smaller is better)   : ',num2str(fSoft)])

K_val = getBlockValue(CL_tuned,'Kafs');

% --- achieved margins at the plant input ----------------------------------
L = getLoopTransfer(CL_tuned,'ud',-1);
[~,MM] = diskmargin(L);
dm = arrayfun(@(m) m.DiskMargin,MM);
[dm_worst,i_worst] = min(dm(:));
gm_worst = 20*log10(MM(i_worst).GainMargin(2));
pm_worst = MM(i_worst).PhaseMargin(2);
disp(['Multiloop disk margin at ''ud'', worst of the ',num2str(length(V_inf)), ...
      ' synthesis velocities (V_inf = ',num2str(V_inf(i_worst)),' m/s): ', ...
      'disk margin ',num2str(dm_worst,'%.4f'),', gain margin ',num2str(gm_worst,'%.2f'), ...
      ' dB, phase margin ',num2str(pm_worst,'%.2f'),' deg'])

% --- verify: closed loop must be stable over the whole velocity range -----
V_check = linspace(20,160,32);
G = build_G_RectWing(V_check,1:8,1:8);
CL = feedback(G,-K_val);
for i_v = 1:length(V_check)
    assert(all(real(eig(CL(:,:,i_v))) < 0), ...
        ['closed loop unstable at V_inf = ',num2str(V_check(i_v)),' m/s'])
end
disp(['Closed loop stable for all V_inf in [',num2str(V_check(1)),', ',num2str(V_check(end)),'] m/s'])


% --- save to the repo data folder -----------------------------------------
RectWing_Cont = ss(K_val);
RectWing_Cont.InputName = P(:,:,1).OutputName(end-7:end);      % the eight u_zN_ddot
RectWing_Cont.OutputName = flap_d_slat_d(:);

script_dir = fileparts(mfilename('fullpath'));
ContPath = fullfile(script_dir,'data','RectWing_Cont_imu18_ail18.mat');
save(ContPath,'RectWing_Cont');

