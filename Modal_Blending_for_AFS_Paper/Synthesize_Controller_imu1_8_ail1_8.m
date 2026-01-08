%% ACC, V_inf=130, 8 IMUs, 8 AIL
disp('---------------------------------------------------------------------')
disp('                    ACC, V_inf=130, 8 IMUs, 8 AIL')
disp('---------------------------------------------------------------------')
%  -5-6-7-8-
% |         |
%  -1-2-3-4-
clearvars
num_modes = 5;
num_poles = 6;
V_inf_ref = 130;
V_inf = [90,100,110,120,130,140,150,160];
ailIDX = 1:8;                 % AIL used
imuIDX = 1:8;               % IMU (acc) used
modesIDX = [1,2];

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
Vu=1;
invWu = invW_theis;

nO=4;
ReqMarg = TuningGoal.Margins('ud',6,45);
% ReqMarg = TuningGoal.Margins('ud',6,30);
noise_u_z_ddot = {'noise_u_z1_ddot','noise_u_z2_ddot','noise_u_z3_ddot','noise_u_z4_ddot','noise_u_z5_ddot','noise_u_z6_ddot','noise_u_z7_ddot','noise_u_z8_ddot'};
flap_d_slat_d  = {'flap1_d','flap2_d','flap3_d','flap4_d','slat1_d','slat2_d','slat3_d','slat4_d'};
ReqContEffNoise = TuningGoal.Gain(noise_u_z_ddot,flap_d_slat_d,invWu*Vu*(1/Vn));
ReqContEffDist = TuningGoal.Gain({'dist_q_f1_ddot','dist_q_f2_ddot'},flap_d_slat_d,invWu*Vu*(1/Vd));
ReqAttenQF = TuningGoal.Gain({'dist_q_f1_ddot','dist_q_f2_ddot'},{'q_f1','q_f2'},Vp*(1/Vd));

rng(20250903,'twister'); 
opt = systuneOptions('RandomStart',7,'UseParallel',true);



%% H_2 Blending nO2
disp('>---------------> H_2 Blending nO=4       -  ACC, V_inf=130, 8 IMUs, 8 AIL')
%% cmG
% G = build_G_RectWing(V_inf_ref,imuIDX,ailIDX);
G = build_G_RectWing(90,imuIDX,ailIDX);
[V,DD] = eig(G.A);
flutIDX = find( (real(diag(DD)) > -1)  & (abs(imag(diag(DD))) < 40) & (abs(imag(diag(DD))) > 5) );
% critIDX = find( (real(diag(DD)) > -10)  & (abs(imag(diag(DD))) < 40) & (abs(imag(diag(DD))) > 5) );
% resIDX = find( (real(diag(DD)) > -10)  & (abs(imag(diag(DD))) < 120) & (abs(imag(diag(DD))) > 40) );
Am = DD;
Bm = V\G.B;
Cm = G.C*V;

Am_flut = Am(flutIDX,flutIDX);
Bm_flut = Bm(flutIDX,:);
Cm_flut = Cm(:,flutIDX);
cmG_flut = ss(Am_flut,Bm_flut,Cm_flut,0);

ky_H2 = [-0.099544537390600;-0.320113549928847;-0.532225871772634;-0.699695797290769;0.076443526375992;0.180003212532740;0.213141230616535;0.176367966086539];
ku_H2 = [0.044984777487334;0.311962575366850;0.447659003359420;0.560303866429342;0.115582305723231;0.268015425550478;0.385831765058990;0.390203827098918];

H2_NORM_cmG_flut = norm(cmG_flut,2);
H2_NORM_ky_cmG_flut_ku = norm(ky_H2'*cmG_flut*ku_H2,2);

tuneCont = ku_H2*tunableSS('nOcont',nO,1,1)*ky_H2';
tuneCL = lft(P,AnalysisPoint('ud',nud)*tuneCont*AnalysisPoint('ym',nym));
% wu wn & wu wd & wp wd
[CL_minQF_H2,fSoft_minQF_H2,gHard_minQF_H2] = systune(tuneCL,[ReqAttenQF,ReqContEffNoise,ReqContEffDist],ReqMarg,opt);
nOcont_H2 = getBlockValue(CL_minQF_H2,'nOcont');
Cont_minQF_H2 = ku_H2*nOcont_H2*ky_H2';



%% Modal Blending
disp('>---------------> Modal Blending nO=4     -  ACC, V_inf=130, 8 IMUs, 8 AIL')
%% cmG
G = build_G_RectWing(V_inf_ref,imuIDX,ailIDX);
% G = build_G_RectWing(90,imuIDX,ailIDX);
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

Cmr_flut = [real(Cm_flut(:,1)), imag(Cm_flut(:,1))];
Cmr_flut_inv=pinv(Cmr_flut);

Bmr_flut = [real(Bm_flut(1,:)); imag(Bm_flut(1,:))];
Bmr_flut_inv=pinv(Bmr_flut);

tuneCont = Bmr_flut_inv*tunableSS('nOcont',nO,2,2)*Cmr_flut_inv;
tuneCL = lft(P,AnalysisPoint('ud',nud)*tuneCont*AnalysisPoint('ym',nym));
% wu wn & wu wd & wp wd
[CL_minQF_MB,fSoft_minQF_MB,gHard_minQF_MB] = systune(tuneCL,[ReqAttenQF,ReqContEffNoise,ReqContEffDist],ReqMarg,opt);
nOcont_MB = getBlockValue(CL_minQF_MB,'nOcont');
Cont_minQF_MB = Bmr_flut_inv*nOcont_MB*Cmr_flut_inv;

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
data_dir = fullfile(script_dir, 'Data');
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end

ContPath = fullfile(data_dir, 'Controller_imu1_8_ail4_advanced_ReqMarg645_test.mat');
save(ContPath, 'Cont_minQF_H2','Cont_minQF_MB');