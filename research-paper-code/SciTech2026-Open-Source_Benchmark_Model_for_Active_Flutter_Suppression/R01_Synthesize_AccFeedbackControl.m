%R01_Synthesize_AccFeedbackControl  Static acceleration-feedback gain behind paper Figs. 12, 15 and 16.
%   Synthesizes a proportional output-feedback gain (IMUs 4 and 8 -> flap 4) at V_inf = 130 m/s with
%   systune: 6 dB / 45 deg margins at the plant input, disturbance-to-flap gain below 28.
%   Inputs:  none (generalized plant from build_P_RectWing). rng is seeded; systune results can
%            still differ slightly between MATLAB releases.
%   Writes:  Data/Controller_imu48_ail4_pof_test.mat (variable Cont_V130_2imu_POF, 1x2 static gain).
%            The _test suffix is on purpose: the shipped Data/Controller_imu48_ail4_pof.mat, which
%            fig12_VpzColormap_ol_pof and fig15_fig16_Vg_Vpz_ol_pof load, is not overwritten.
%            To plot your own result, point ContPath in those two fig scripts to the _test file.
%   Runtime: about a minute (estimate). Requires Control System Toolbox (systune).

%% ACC, V_inf=130, 2 IMUs, 1 AIL
disp('---------------------------------------------------------------------')
disp('                    ACC, V_inf=130, 2 IMUs, 1 AIL')
disp('---------------------------------------------------------------------')
%  -5-6-7-8-
% |         |
%  -1-2-3-4-
clearvars
num_modes = 5;
V_inf = 130;
ailIDX = 4;                 % AIL used
imuIDX = [4,8];               % IMU (acc) used
modesIDX = [1,2];

P = build_P_RectWing(V_inf,imuIDX,ailIDX,modesIDX);
nym = length(imuIDX);
nud = length(ailIDX);

disp('>---------------> Proportional Output Feedback  -  ACC, V_inf=130, 2 IMUs, 1 AIL')
tuneCont = tunableGain('Pcont',nud,nym);
tuneCL = lft(P,AnalysisPoint('ud',nud)*tuneCont*AnalysisPoint('ym',nym));
ReqMarg1 = TuningGoal.Margins('ud',6,45);
ReqAtten1 = TuningGoal.Gain({'dist_q_f1_ddot','dist_q_f2_ddot'},{'flap4_d'},28);     % constrains the largest singular value of the transfer matrix from inputname to outputname
rng(20250903,'twister'); 
% opt = systuneOptions('RandomStart',7,'UseParallel',true);
[CL,fSoft,gHard] = systune(tuneCL,[ReqAtten1],ReqMarg1);        % Hard Goals are not further optimized if met. They are treated as constraints!!!
PCont = getBlockValue(CL,'Pcont');
Cont_V130_2imu_POF = PCont;


script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
data_dir = fullfile(script_dir, 'Data');
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end

ContPath = fullfile(data_dir, 'Controller_imu48_ail4_pof_test.mat');
save(ContPath, 'Cont_V130_2imu_POF');