function M_ele = build_M_ele(y1, y2, rho_bar, I_Tym, I_zm)
% build_M_ele  Element Mass Matrix of a bending–torsion beam with 6
%               DoFs:  - u_z 1        left node bending displacement
%                      - du_z/dy 1    left node bending slope
%                      - Psi_y 1      left node torsion angle
%                      - u_z 2        right node bending displacement
%                      - du_z/dy 2    right node bending slope
%                      - Psi_y 2      right node torsion angle
%
% INPUT
%   y1, y2   : nodal physical coordinates  (y2 > y1)
%   rho_bar  : mass per length              int_A rho dA
%   I_Tym    : polar mass moment            int_A rho r^2 dA
%   I_zm     : deviational (coupling) term  int_A rho x  dA 
%
% OUTPUT
%   M_ele    : 6×6 element mass matrix
%
% Jonas . July 2025
% _________________________________________________________________________


% ---- 4-point Gauss–Legendre rule
eps_GL = [-0.8611363115940526;
          -0.3399810435848563;
           0.3399810435848563;
           0.8611363115940526]; 
w_GL   = [ 0.3478548451374538;
           0.6521451548625461;
           0.6521451548625461;
           0.3478548451374538];
W  = diag(w_GL); 


% ---- eps Shape functions (cubic Hermite + linear Lagrange)

L = y2 - y1;
S = [ 0.25*( 2 - 3*eps_GL +   eps_GL.^3),   ...               % S1(eps)
      (L/8)*(1 -   eps_GL -   eps_GL.^2 + eps_GL.^3),   ...   % S2(eps)
      0.5 *(1 -   eps_GL),                ...                 % S3(eps)
      0.25*( 2 + 3*eps_GL -   eps_GL.^3),   ...               % S4(eps)
      (L/8)*(-1 -  eps_GL +   eps_GL.^2 + eps_GL.^3), ...     % S5(eps)
      0.5 *(1 +   eps_GL) ];                                  % S6(eps)
% S is a 4×6 matrix: rows = Gauss points, cols = DOFs

% ---- Galerkin (energy) integrals  ------------------------------------
%
% M_ij = int_L I_Tym    N_i(y) N_j(y)  dy
%      + int_L rho_bar  H_i(y) H_j(y)  dy
%      - int_L I_zm   ( N_i(y)H_j(y)+H_i(y)N_j(y) )  dy
%
%      negative sign on coupling term because positive torsion around
%      y-axis and positive leverarm in x-dir. leads to movement in negative
%      z-dir.
% 
Int = S.' * W * S;                                 % 6×6, symmetric

% ---- Coord. Trafo. J = L/2 multiplies every integral
J = L/2;

bend = [1 2 4 5]; 
tors = [3 6];

M_ele = zeros(6);

% ---- bending block
M_ele(bend,bend) = rho_bar * J * Int(bend,bend);

% ---- torsion block
M_ele(tors,tors) = I_Tym * J * Int(tors,tors);

% ---- coupling blocks
M_ele(tors,bend) = -I_zm * J * Int(tors,bend);
M_ele(bend,tors) =  M_ele(tors,bend).';        % enforce symmetry


end
