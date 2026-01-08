function K_ele = build_K_ele(y1, y2, EI_xxa, GI_Tya)
% build_K_ele   Element Stiffness Matrix of a bending–torsion beam with 6
%               DoFs:  - u_z 1        left node bending displacement
%                      - du_z/dy 1    left node bending slope
%                      - Psi_y 1      left node torsion angle
%                      - u_z 2        right node bending displacement
%                      - du_z/dy 2    right node bending slope
%                      - Psi_y 2      right node torsion angle
%
% INPUT
%   y1, y2   : nodal physical coordinates (y2 > y1)
%   EI_xxa   : bending stiffness  (E I_xxa)   I_xxa = int_A z^2 dA
%   GI_Tya   : torsional stiffness (G I_Tya)  I_Tya = int_A r^2 dA
%
% OUTPUT
%   K_ele    : 6×6 element stiffness matrix
%
% Jonas . July 2025
% _________________________________________________________________________

% ---- 4-point Gauss–Legendre rule on [-1,1]
eps_GL = [-0.8611363115940526;
          -0.3399810435848563;
           0.3399810435848563;
           0.8611363115940526];
w_GL   = [ 0.3478548451374538;
           0.6521451548625461;
           0.6521451548625461;
           0.3478548451374538];
W = diag(w_GL);

% ---- geometry-dependent scaling factors
L      = y2 - y1;                         % element length
fact_b = 8 / L^3;                         % bending   multiplier
fact_t = 2 / L;                           % torsional multiplier

% ---- eps-second derivatives for bending DOFs (4×6) ---------------------
d2S = [  1.5*eps_GL, ...                       % S1''(eps)
        (L/4)*(3*eps_GL - 1), ...              % S2''(eps)
         zeros(4,1), ...                       % S3''(eps) = 0 not used
        -1.5*eps_GL, ...                       % S4''(eps)
        (L/4)*(3*eps_GL + 1), ...              % S5''(eps)
         zeros(4,1) ];                         % S6''(eps) = 0 not used

% ---- eps-first derivatives for torsion DOFs (2×6) ----------------------
dS  = [ (-3 + 3*eps_GL.^2)/4, ...              % S1'(eps) not used
        (L/8)*(3*eps_GL.^2 - 2*eps_GL - 1), ...% S2'(eps) not used
        -0.5*ones(4,1), ...                    % S3'(eps)
        (3 - 3*eps_GL.^2)/4, ...               % S4'(eps) not used
        (L/8)*(3*eps_GL.^2 + 2*eps_GL - 1), ...% S5'(eps) not used
         0.5*ones(4,1) ];                      % S6'(eps)

% ---- Galerkin (energy) integrals  ------------------------------------
%
% K_ij = int_L E I_xxa H_i''(y) H_j''(y) dy
%      + int_L G I_Tya N_i'(y)  N_j'(y)  dy
%
% The factors fact_b, fact_t already include both J and the chain-rule
% powers (4/L^2 or 2/L), so we only need  dS  or d2S  and the weights.

K_ele = zeros(6);                             % full 6×6 matrix

bend = [1 2 4 5];
K_ele(bend,bend) = EI_xxa * fact_b * ( d2S(:,bend).' * W * d2S(:,bend) );

tors = [3 6];
K_ele(tors,tors) = GI_Tya * fact_t * ( dS(:,tors).'  * W * dS(:,tors)  );

end
