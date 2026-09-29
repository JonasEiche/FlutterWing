%% R00: Synthesize the Goland AFS Demo Controller
%
% Regenerates data/Goland_Cont_nO4_lpf.mat - the pre-synthesized active
% flutter suppression controller that closes the loop in the tutorial
% (docs/TUTORIAL.md chapters 10 to 12, TUTORIAL.m sections (10) and (11))
% and in docs/figures/virtual_flight.gif.
% 
% SISO Goland case: 1 IMU, 1 trailing edge flap.
%
% Controller structure:
%   Goland_Cont_nO = nOcont(s) * G_lpf(s)
%   nOcont : tunable state space block, order nO = 4
%   G_lpf  : fixed first order low-pass for high frequency roll-off
%
% Tuning goals (soft, minimized):
%   1. ReqAttenQF      : gain dist_q_f{1,2}_ddot -> q_f{1,2}
%                        (suppress flutter mode response to modal forces)
%   2. ReqContEffNoise : weighted gain noise_u_z1_ddot -> flap1_d
%                        (limit actuator wear from sensor noise)
%   3. ReqContEffDist  : weighted gain dist_q_f{1,2}_ddot -> flap1_d
%                        (limit actuator activity during disturbance rejection)
% Hard constraint:
%   4. ReqMarg         : disk margins 6 dB / 45 deg at the plant input 'ud'
%
% The actuator effort weight W_theis concentrates the allowed control
% activity in the flutter frequency band. For the Goland wing the wind-off
% modes sit at 48 and 96 rad/s and the flutter mode at ~70 rad/s, so the
% corner frequencies are set to w1 = 20, w2 = 150 rad/s (the RectWing
% baseline uses 12 and 64 rad/s for its ~28 rad/s flutter mode).


clearvars
rng(20260704,'twister');            % reproducibility

num_modes = 5;
num_poles = 6;
imuIDX = 1;
ailIDX = 1;
modesIDX = [1,2];                   % target the first two structural modes
V_inf = [90,110,130,150,170,190];   % LPV synthesis grid (m/s), flutter onset ~169 m/s

P = build_P_Goland(V_inf,imuIDX,ailIDX,modesIDX);
nym = length(imuIDX);
nud = length(ailIDX);

% actuator effort weight (Theis bandpass-like filter) centered on the flutter band
w1 = 20;
w2 = 150;
s = tf('s');
invW_theis = ((s+0.01*w1)*(0.01*s+w2))/((s+w1)*(s+w2));

Vp=0.2;
Vn=0.1;
Vd=0.5;
Vu=0.5;
invWu = invW_theis;

ReqMarg = TuningGoal.Margins('ud',6,45);
ReqContEffNoise = TuningGoal.Gain({'noise_u_z1_ddot'},{'flap1_d'},invWu*Vu*(1/Vn));
ReqContEffDist = TuningGoal.Gain({'dist_q_f1_ddot','dist_q_f2_ddot'},{'flap1_d'},invWu*Vu*(1/Vd));
ReqAttenQF = TuningGoal.Gain({'dist_q_f1_ddot','dist_q_f2_ddot'},{'q_f1','q_f2'},Vp*(1/Vd));

% controller: tunable 4th order block + fixed first order low-pass roll-off
nO = 4;
w_lp = 300;                         % roll-off corner (rad/s), above the actuator bandwidth
G_lpf = w_lp/(s+w_lp);
tuneCont = tunableSS('nOcont',nO,1,1)*G_lpf;
tuneCL = lft(P,AnalysisPoint('ud',nud)*tuneCont*AnalysisPoint('ym',nym));

opt = systuneOptions('RandomStart',3,'UseParallel',false,'Display','final');
[CL_tuned,fSoft,gHard] = systune(tuneCL,[ReqAttenQF,ReqContEffNoise,ReqContEffDist],ReqMarg,opt);
disp(['Soft goals (smaller is better)  : ',num2str(fSoft)])
disp(['Hard goal  (must be <= 1)       : ',num2str(gHard)])

Goland_Cont_nO = getBlockValue(CL_tuned,'nOcont')*G_lpf;

% --- verify: closed loop must be stable over the whole velocity range -----
V_check = linspace(20,190,42);
G = build_G_Goland(V_check,imuIDX,ailIDX);
CL = feedback(G,-Goland_Cont_nO);
for i_v = 1:length(V_check)
    assert(all(real(eig(CL(:,:,i_v))) < 0), ...
        ['closed loop unstable at V_inf = ',num2str(V_check(i_v)),' m/s'])
end
disp(['Closed loop stable for all V_inf in [',num2str(V_check(1)),', ',num2str(V_check(end)),'] m/s'])

% --- save to the repo data folder -----------------------------------------
script_dir = fileparts(mfilename('fullpath'));
ContPath = fullfile(script_dir,'data','Goland_Cont_nO4_lpf.mat');
save(ContPath,'Goland_Cont_nO');
disp(['Saved controller to ',ContPath])
