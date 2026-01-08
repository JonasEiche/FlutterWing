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