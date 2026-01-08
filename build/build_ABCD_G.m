function [A_noS,B_noS,C_noS,D_noS] = build_ABCD_G(V_inf,Structure,Aero)

Kff = Structure.Kff;
Mff = Structure.Mff;
Dff = Structure.Dff;
Sfj = Structure.Sfj;
DRe_jf = Structure.DRe_jf;
DIm_jf = Structure.DIm_jf;
PHIgf = Structure.PHIgf;
DRe_jx = Structure.DRe_jx;
DIm_jx = Structure.DIm_jx;
PHIzg = Structure.PHIzg;

c_ref = Aero.c_ref;
rho = Aero.rho;
poles = Aero.poles;
Q0jj = Aero.Q0jj;
Q1jj = Aero.Q1jj;
QLpjj = Aero.QLpjj;

num_modes = size(Kff,1);
num_poles = length(poles);
num_panels = size(Q0jj,1);
num_AIL = size(DRe_jx,2);

%% build_A ----------------------------------------------------------------
% Right Realization with num_modes*num_poles aerodynamic lag states
num_x_L         = num_modes*num_poles;
% q_bar           = 0.5*rho*V_inf^2;

Mff_til         = Mff-0.5*rho*Sfj*Q1jj*0.5*c_ref*DIm_jf;

sumQLpjjKff_til = zeros(num_panels,num_modes);
for i = 1:num_poles
    sumQLpjjKff_til = sumQLpjjKff_til + QLpjj(:,:,i)*(DRe_jf-DIm_jf*poles(i)*2/c_ref);
end
Kff_til         = 0.5*rho*V_inf^2*Sfj*(Q0jj*DRe_jf+sumQLpjjKff_til)-Kff;

sumQLpjjDff_til = zeros(num_panels,num_modes);
for i = 1:num_poles
    sumQLpjjDff_til = sumQLpjjDff_til + QLpjj(:,:,i)*DIm_jf;
end
Dff_til         = 0.5*rho*Sfj*V_inf*( Q0jj*DIm_jf+Q1jj*0.5*c_ref*DRe_jf + sumQLpjjDff_til )-Dff;

% D_til           = repmat(speye(num_modes),1,num_poles);
D_til           = repmat(eye(num_modes),1,num_poles);
% R_til      = diag(sparse(reshape(repmat(-poles, num_modes,1),1,num_modes*num_poles)))*(2*V_inf/c_ref);
R_til      = diag(reshape(repmat(-poles, num_modes,1),1,num_modes*num_poles))*(2*V_inf/c_ref);

E_qg_til            = zeros(num_x_L,num_modes);
for i = 1:num_poles
    E_qg_til((i-1)*num_modes+1:i*num_modes,1:num_modes) = 0.5*rho*Sfj*V_inf^3*QLpjj(:,:,i)*(DIm_jf*(poles(i)*2/c_ref)^2 - DRe_jf*poles(i)*2/c_ref);
end

A_noS =               [zeros(num_modes), eye(num_modes),            zeros(num_modes,num_x_L);
                   Mff_til\Kff_til,  Mff_til\Dff_til,           Mff_til\D_til;
                   E_qg_til,         zeros(num_x_L,num_modes),  R_til];

%% build_B ----------------------------------------------------------------
E_ux_til            = zeros(num_x_L,num_AIL);
for i = 1:num_poles
    E_ux_til((i-1)*num_modes+1:i*num_modes,1:num_AIL) = 0.5*rho*Sfj*V_inf^3*QLpjj(:,:,i)*(DIm_jx*(poles(i)*2/c_ref)^2 - DRe_jx*poles(i)*2/c_ref);
end

sumQLpjjBjx = zeros(num_panels,num_AIL);
for i = 1:num_poles
    sumQLpjjBjx = sumQLpjjBjx + QLpjj(:,:,i)*(DRe_jx-DIm_jx*poles(i)*2/c_ref);
end
Bgx         = 0.5*rho*V_inf^2*Sfj*(Q0jj*DRe_jx+sumQLpjjBjx);

sumQLpjjBjx_dot = zeros(num_panels,num_AIL);
for i = 1:num_poles
    sumQLpjjBjx_dot = sumQLpjjBjx_dot + QLpjj(:,:,i)*DIm_jx;
end
Bgx_dot         = 0.5*rho*Sfj*V_inf*( Q0jj*DIm_jx+Q1jj*0.5*c_ref*DRe_jx + sumQLpjjBjx_dot );
% Bgx_ddot        = 0.5*rho*Sfj*Q1jj*0.5*c_ref*DIm_jx;
Bgx_ddot        = zeros(num_modes,num_AIL);             % TO AVOID FEEDTHROUGH FROM ACTUATOR DEMAND TO ACCELERATION SENSOR


B_noS =               [zeros(num_modes,3*num_AIL);
                   Mff_til\Bgx, Mff_til\Bgx_dot, Mff_til\Bgx_ddot;
                   E_ux_til, zeros(num_x_L,2*num_AIL)];




%% build_CD ---------------------------------------------------------------
D_til           = repmat(eye(num_modes),1,num_poles);
% D_til           = repmat(speye(num_modes),1,num_poles);

C_noS = [ PHIzg*PHIgf*(Mff_til\Kff_til), PHIzg*PHIgf*(Mff_til\Dff_til), PHIzg*PHIgf*(Mff_til\D_til) ];

D_noS = [ PHIzg*PHIgf*(Mff_til\Bgx),PHIzg*PHIgf*(Mff_til\Bgx_dot),PHIzg*PHIgf*(Mff_til\Bgx_ddot) ];