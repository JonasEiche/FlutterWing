function S_ele = build_S_ele(E1,E2,Ps)
% build_S_ele   Map aerodynamic panel forces applied at the
%               lifting-line point (LP) to the 6 beam-element DOFs
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
%              Lift is assumed to act at the quarter chord point’ of the panel:
%                   BVP1 = 0.25 P2 + 0.75 P1   (left edge)
%                   BVP2 = 0.25 P3 + 0.75 P4   (right edge)
%                   LP   = 0.5 (BVP1 + BVP2)
%
%              Only panels whose LP projects onto the current element
%              orthogonally contribute to S_ele.
%
% OUTPUT
%   S_ele    : 6 × num_panel matrix.  Each column holds the shape-function
%              weights that distribute aero-forces to the element DOFs:
%
%                       F_ele = S_ele * (q_bar Dp_j c s)
%
%              where  Dp_j  is the pressure difference coefficient at each aero panel,
%              q_bar  the dynamic pressure,   c  the panel chord,  and  s  the panel span (width).
% 
% Sign convention
% _________________________________________________________________________
%   y_s = [0;1;0] is the direction of the beam in structural coord. system
%   z_s = [0;0;1] is the bending motion of the beam in structural coords.
%
%
% Jonas * July 2025
% _________________________________________________________________________

y_s = [0;1;0];  % Hard Coded Direction of beam 
% y_s = (E2-E1)/norm(E2-E1);
z_s = [0;0;1];  % Hard Coded Direction of bending motion u_z
L   = norm(E2-E1);
np = length(Ps);
S_ele    = zeros(6, np);   

for i_p = 1:np
    BVP1 = 0.25 * Ps{i_p}{2} + 0.75 * Ps{i_p}{1};
    BVP2 = 0.25 * Ps{i_p}{3} + 0.75 * Ps{i_p}{4};
    LP   = 0.5  * (BVP1 + BVP2);
    N = cross((Ps{i_p}{2}-Ps{i_p}{1}),(Ps{i_p}{4}-Ps{i_p}{1}));
    N = N/norm(N);
    c = norm(0.5*(Ps{i_p}{2}-Ps{i_p}{1}+Ps{i_p}{3}-Ps{i_p}{4}));
    S = Ps{i_p}{4}-Ps{i_p}{1} - (Ps{i_p}{2}-Ps{i_p}{1})/norm((Ps{i_p}{2}-Ps{i_p}{1}))^2*dot((Ps{i_p}{2}-Ps{i_p}{1}),(Ps{i_p}{4}-Ps{i_p}{1}));
    s = norm(S);
    % S = S/s;
    a = s*c;

    t = dot(E2-E1,LP-E1)/norm(E2-E1)^2;
    if ~(0<t && t<=1) % skip panels whose LP is not inside this element; include last node exclude first
        continue
    end
    % ---------------------------------------------------------------------
    % Shape-function values in eps = 2·(y-y1)/L − 1
    % ---------------------------------------------------------------------
    eps = 2 * ( dot(E2-E1,LP-E1)/norm(E2-E1) ) / L - 1;                     % local coord. in [-1,1]
    % eps = 2*t-1
    S1 = 0.25 * ( 2 - 3*eps + eps^3 );
    S2 = (L/8) * ( 1 - eps - eps^2 + eps^3 );
    S3 = 0.5  * ( 1 - eps );
    S4 = 0.25 * ( 2 + 3*eps - eps^3 );
    S5 = (L/8) * (-1 - eps + eps^2 + eps^3 );
    S6 = 0.5  * ( 1 + eps );


    S_ele(:, i_p) = [  S1 * a * dot(N,z_s) ;
                       S2 * a * dot(N,z_s) ;
                       S3 * a * dot(y_s, cross(LP-E1,N) ) ;
                       S4 * a * dot(N,z_s) ;
                       S5 * a * dot(N,z_s) ;
                       S6 * a * dot(y_s, cross(LP-E1,N) ) ];
end

end