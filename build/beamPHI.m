function [PHIuzg, PHIpsiyg] = beamPHI(y_pos, E, ele)
% beamPHIztg   Shape-function matrices that interpolate the global
%              displacement vector q_g to arbitrary spanwise positions y.
%
%        u_z(y)       =  PHIzg(y,:) · q_g        (vertical bending deflection)
%        Psi_y(y)     =  PHItg(y,:) · q_g        (torsion angle about the y-axis)
%
% INPUT
%   y_pos   : column vector of spanwise stations where the response is needed
%
%   E       : cell array, one entry per element.  
%             For element e the nodal y-coordinates are stored as
%                     y1 = E{e}{1}(2),    y2 = E{e}{2}(2).
%
%   ele     : cell array, one entry per element.  
%             ele{e} is a 1 × 6 vector with the global DOF numbers
%             corresponding to element e.
%
% OUTPUT
%   PHIzg   : length(y_pos) × 3·(num_ele+1)   (maps to u_z )
%   PHItg   : length(y_pos) × 3·(num_ele+1)   (maps to Psi_y )
%
% Element formulation
% _________________________________________________________________________
% () Cubic Hermite interpolation for the bending DOFs (u_z, du_z/dy):
%       H_1(eps) = 0.25( 2 − 3eps + eps^3 )             * u_z1
%       H_2(eps) = 0.125L( 1 − eps − eps^2 + eps^3 )    * du_z1/dy
%       H_3(eps) = 0.25( 2 + 3eps − eps^3 )             * u_z2
%       H_4(eps) = 0.125L(−1 − eps + eps^2 + eps^3 )    * du_z2/dy
%
% () Linear interpolation for torsion DOFs (Psi_y):
%       N_1(eps) = 0.5(1 − eps) · Psi_y1 ,     N_2(eps) = 0.5(1 + eps) · Psi_y2
%
%   with eps \in [−1,1] given by  eps = 2(y − y_1)/L − 1.
%
%
% Jonas * July 2025
% _________________________________________________________________________

y_pos = y_pos(:);                         % works for row or column input

num_ele      = numel(E);
num_DoF_glob = 3 * (num_ele + 1);

PHIuzg = zeros(numel(y_pos), num_DoF_glob);
PHIpsiyg = zeros(numel(y_pos), num_DoF_glob);

for e = 1:num_ele
    % ---- element span and the y-stations that lie inside it -------------
    y1   = E{e}{1}(2);
    y2   = E{e}{2}(2);
    L    = y2 - y1;

    mask = (y_pos >= y1) & (y_pos <= y2);   % row selector
    if ~any(mask), continue, end            % skip empty ranges

    y_e  = y_pos(mask);                     % local slice
    eps  = 2 * (y_e - y1) / L - 1;          % eps-coordinates

    % ---- cubic Hermite shape functions (bending) ------------------------
    H1 = 0.25 * ( 2 - 3*eps + eps.^3 );
    H2 = (L/8) * ( 1 - eps - eps.^2 + eps.^3 );
    H3 = 0.25 * ( 2 + 3*eps - eps.^3 );
    H4 = (L/8) * (-1 - eps + eps.^2 + eps.^3 );

    % ---- linear shape functions (torsion) -------------------------------
    N1 = 0.5 * (1 - eps);
    N2 = 0.5 * (1 + eps);

    % ---- assemble rows into the global matrices -------------------------
    dof_w  = ele{e}([1 2 4 5]); 
    dof_t  = ele{e}([3 6]);

    PHIuzg(mask, dof_w) = [H1, H2, H3, H4];
    PHIpsiyg(mask, dof_t) = [N1, N2];
end
end
