function [M,d] = h2_opt_input_siso(B,w)
%H2_OPT_INPUT_SISO  Real input direction maximising the energy injected into one complex pole.
%
%   [M,d] = h2_opt_input_siso(B,w)
%   NOTE: the direction d is the SECOND output; the MIMO variants
%   h2_opt_input_mimo / h2_opt_output_mimo return the direction FIRST.
%
%   Inputs
%   ------
%   B : (n x n_u) real    - plant input matrix
%   w : (n x 1)  complex  - left eigenvector of the pole
%
%   Outputs
%   -------
%   M : (n_u x n_u) real  - rank-2 matrix  M = a a^T + b b^T
%   d : (n_u x 1) real    - unit-norm direction that maximises |w^T B d|_2
%
%   The function follows the steps
%     1)  z = w^T B                       (row vector, 1 x n_u)
%     2)  a = Re(z)^T ,  b = Im(z)^T       (column vectors)
%     3)  M = a a^T + b b^T               (real, symmetric, ≤ rank-2)
%     4)  d = principal eigenvector(M).  (Rayleigh–Ritz theorem)
%
%   Reference
%   ---------
%   Derived from  |w^H Bd|_2  = d^T B^T Re(w w^T) B d
%                 with  w w^T  = (a+jb)(a+jb)^T
%
%   Scope
%   -----
%   Separate input/output method ("Method 2" in R6_siso_rhp_zeros.m and
%   R7_quadrature_mismatch.m of the AS2026 paper code). The joint H2-optimal
%   blending vectors reported in the papers come from proprietary code and
%   are hard-coded in those scripts; the h2_opt_* functions do not reproduce
%   them.
%

% ---------- Step 1: complex coupling row -------------------------------
z = (w') * B;                % 1 x n_u complex row (w' is Hermitian transpose)
% ---------- Step 2: split into real column vectors ----------------------
a = real(z).';               % n_u x 1
b = imag(z).';               % n_u x 1
% ---------- Step 3: build the  rank-2  matrix ---------------------------
M = a*a.' + b*b.';           % n_u x n_u, real symmetric
% ---------- Step 4: principal eigenvector (Rayleigh–Ritz) ---------------
% For numerical robustness we ask eig for only the largest eigenpair
if size(M,1) == 1            % scalar actuator case
    d = 1;                   % M is 1x1; the only direction is that actuator
else
    [v,~] = eigs(M,1,'la');  % 'la' = largest algebraic eigenvalue
    d = real(v);             % remove negligible imaginary round-off
    d = d / norm(d);         % ensure unit length
end
end
