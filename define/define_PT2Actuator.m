function [G_act] = define_PT2Actuator(num_AIL, f0_Hz)
% define_PT2Actuator  Second-order (PT2) actuator models, one per control surface, channels flap1..flapN.
%   G_act = define_PT2Actuator(num_AIL) returns a 1 x num_AIL cell array of ss models.
%   G_act = define_PT2Actuator(num_AIL, f0_Hz) sets the bandwidth w0 = 2*pi*f0_Hz (default 32 Hz;
%   TUTORIAL.m section (8) also builds a 10 Hz one).
%   Used by the Goland entry points (build_G_Goland, build_P_Goland, build_LPV_*_Goland).
%   Same as define_RectWing_PT2Actuator except that every channel is named flapN
%   (the RectWing version names channels 5..8 slatN).
%
%   G_act{i}: input   flapN_d                          demanded deflection (rad)
%             outputs [flapN; flapN_dot; flapN_ddot]   deflection (rad), rate (rad/s),
%                     acceleration (rad/s^2) -> the deflection / rate / acceleration
%                     input blocks of build_ABCD_G / build_ABCD_P
%             states  [flapN; flapN_dot], the rate state is scaled by w0 for
%                     conditioning (StateScale below)
%
%   PT2 actuator model with K = 1, d = 0.9, w0 = 2*pi*f0_Hz rad/s (32 Hz by default):
%       G_act(s) = K*w0^2 / (s^2 + 2*d*w0*s + w0^2)     (flapN_d -> flapN)
%
% Jonas * July 2025
% _________________________________________________________________________
if nargin < 2, f0_Hz = 32; end
K=1; d=0.9; w0 = 2*pi*f0_Hz;

G_act = cell(1,num_AIL);
for i = 1:num_AIL
    Aact = [0,     1;
            -w0^2, -2*d*w0];
    Bact = [0; K*w0^2];
    Cact = [1,      0;
            0,      1;
            -w0^2,  -2*d*w0];
    Dact = [0; 0; K*w0^2];

    StateScale = [1; w0];
    G_act{i} = ss(diag(1./StateScale)*Aact*diag(StateScale),diag(1./StateScale)*Bact,Cact*diag(StateScale),Dact);
end

for i = 1:num_AIL
    G_act{i}.InputName = ['flap',num2str(i),'_d'];
    G_act{i}.OutputName = {['flap',num2str(i)],['flap',num2str(i),'_dot'],['flap',num2str(i),'_ddot']};
    G_act{i}.StateName = {['flap',num2str(i)],['flap',num2str(i),'_dot']};
end