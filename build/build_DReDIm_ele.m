function [DRe_ele, DIm_ele] = build_DReDIm_ele(E1,E2,Ps)
% build_DReDIm_ele   Maps FEM structural deformation of bending-torsion
%                    beam onto the boundary conditions of the DLM.
%                    The 6 beam-element DOFs to Pistolesi Point (LP) downwash condition 
%
% INPUT
%   E1, E2   : physical coordinates of the FEM element nodes (column vectors)
%
%   Ps       : cell array, one entry per aerodynamic panel geometry in
%              structural coords.
%              For each panel the four corner nodes are stored as
% 
%              P{panel}{node} = [ x_pos;         1---4
%                                 y_pos;         |   |
%                                 z_pos ]        2---3
%
%              Downwash at the 3/4-point of each panel w_j depends on
%              structural deformation DRe*q_g and velocity DIm*q_g_dot
%              Only panels whose JP projects onto the current element
%              orthogonally contribute to DRe_ele / DIm_ele.
%
% OUTPUT
%   DRe_ele,  
%   DIm_ele  : num_panel × 6 matrix
%
%              w_j = DRe_ele           * q_ele 
%                  + DIm_ele * 1/V_inf * q_ele_dot
%
%              w_j = u_z/V_inf normalized downwash velocity at control point JP 
% 
% Sign convention
% _________________________________________________________________________
%   y_s = [0;1;0] is the direction of the beam in structural coord. system
%   z_s = [0;0;1] is the bending motion of the beam in structural coords.
%
%
% Jonas * July 2025
% _________________________________________________________________________
y_s = [0;1;0];
% y_s = (E2-E1)/norm(E2-E1);
z_s = [0;0;1];
L   = norm(E2-E1);
np = length(Ps);

DRe_ele = zeros(np, 6);
DIm_ele = zeros(np, 6); 
for i_p = 1:np
    JP1 = 0.75 * Ps{i_p}{2} + 0.25 * Ps{i_p}{1};
    JP2 = 0.75 * Ps{i_p}{3} + 0.25 * Ps{i_p}{4};
    JP  = 0.5  * (JP1 + JP2);
    N = cross((Ps{i_p}{2}-Ps{i_p}{1}),(Ps{i_p}{4}-Ps{i_p}{1}));
    N = N/norm(N);
    % c = norm(0.5*(Ps{p}{2}-Ps{p}{1}+Ps{p}{3}-Ps{p}{4}));
    S = Ps{i_p}{4}-Ps{i_p}{1} - (Ps{i_p}{2}-Ps{i_p}{1})/norm((Ps{i_p}{2}-Ps{i_p}{1}))^2*dot((Ps{i_p}{2}-Ps{i_p}{1}),(Ps{i_p}{4}-Ps{i_p}{1}));
    s = norm(S);
    S = S/s;

    t = dot(E2-E1,JP-E1)/norm(E2-E1)^2;
    if ~(0<t && t<=1) % skip panels whose LP is not inside this element; include last node exclude first
        continue
    end
    % ---------------------------------------------------------------------
    % Shape-function values in eps = 2·(y-y1)/L − 1
    % ---------------------------------------------------------------------
    eps = 2 * ( dot(E2-E1,JP-E1)/norm(E2-E1) ) / L - 1;                     % local coord. in [-1,1]
    % eps = 2*t-1

    S1 = 0.25 * ( 2 - 3*eps + eps^3 );
    S2 = (L/8) * ( 1 - eps - eps^2 + eps^3 );
    S3 = 0.5  * ( 1 - eps );
    S4 = 0.25 * ( 2 + 3*eps - eps^3 );
    S5 = (L/8) * ( -1 - eps + eps^2 + eps^3 );
    S6 = 0.5  * ( 1 + eps );

    DRe_ele(i_p, 3) =  S3 * dot(y_s,S);
    DRe_ele(i_p, 6) =  S6 * dot(y_s,S);
    
    DIm_ele(i_p, :) = [ S1 * dot(-N,z_s) , ...
                        S2 * dot(-N,z_s) , ...
                        S3 * dot(-N, cross(y_s,JP-E1)) , ...
                        S4 * dot(-N,z_s) , ...
                        S5 * dot(-N,z_s) , ...
                        S6 * dot(-N, cross(y_s,JP-E1)) ];
end
end