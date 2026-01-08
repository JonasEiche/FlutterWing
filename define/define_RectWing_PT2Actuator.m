
function [G_act] = define_RectWing_PT2Actuator(num_AIL)
% ---- PT2 Actuators ------------------------------------------------------
% PT2 Actuator Model with K=1, d=1, w0 = 16 Hz = 16*2*pi rad/s
% G_act(s) = K*w0^2 / (s^2 + 2*d*w0*s + w0^2)
K=1; d=0.9; w0 = 32*2*pi;

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
    % StateScale = [1; 1];
    G_act{i} = ss(diag(1./StateScale)*Aact*diag(StateScale),diag(1./StateScale)*Bact,Cact*diag(StateScale),Dact);
end

for i = 1:4
    G_act{i}.InputName = ['flap',num2str(i),'_d'];
    G_act{i}.OutputName = {['flap',num2str(i)],['flap',num2str(i),'_dot'],['flap',num2str(i),'_ddot']};
    G_act{i}.StateName = {['flap',num2str(i)],['flap',num2str(i),'_dot']};
    G_act{4+i}.InputName = ['slat',num2str(i),'_d'];
    G_act{4+i}.OutputName = {['slat',num2str(i)],['slat',num2str(i),'_dot'],['slat',num2str(i),'_ddot']};
    G_act{4+i}.StateName = {['slat',num2str(i)],['slat',num2str(i),'_dot']};
end