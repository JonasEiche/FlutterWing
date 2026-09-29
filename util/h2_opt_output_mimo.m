function [H,H2] = h2_opt_output_mimo(C,V)
%H2_OPT_OUTPUT_MIMO  Real r-dimensional sensor projector maximising the output energy of r modes.
%
%   [H,H2] = h2_opt_output_mimo(C,V)
%   NOTE: the projector H is the FIRST output; the SISO variants
%   h2_opt_input_siso / h2_opt_output_siso return the direction SECOND.
%   Requires r = size(V,2) <= n_y (eigs(M,r) errors otherwise).
%
%   Inputs
%   ------
%     C : (n_y × n)   real     – plant output matrix  (y = C x)
%     V : (n   × r)   complex  – right eigenvectors of every mode of interest
%                                * For each real pole include its real v.
%                                * For a complex pair include BOTH  v  and  v_bar.
%
%   Outputs
%   -------
%     H  : (n_y × r)  real     – projector with orthonormal columns (H.'*H = I)
%     H2 : (1   × r)  real     – H2 norm of every column h_j, i.e.
%                                |h_j^T C (sI-A)^-1 B|_2 
%
%   Method
%   ------
%   For each mode  k
%       z_k = C v_k             (n_y x 1 complex)
%       a_k = Re(z_k),  b_k = Im(z_k)
%   Accumulate
%       M = sum_k (a_k a_k^T + b_k b_k^T)   (n_y × n_y, real, rank ≤ 2r)
%   The r principal eigenvectors of M give the optimal projector H.
%   Per-column energy:  H2(j) = sqrt(h_j^T M h_j).
%
%   To obtain "true" H^T energy  |h^T C v_k|_2/(−2sigma_k)  weight each rank-2 term
%   by 1/(−2*real(lambda_k))   inside the loop.
%
%   Scope
%   -----
%   Separate input/output method ("Method 2" in R6_siso_rhp_zeros.m and
%   R7_quadrature_mismatch.m of the AS2026 paper code). The joint H2-optimal
%   blending vectors reported in the papers come from proprietary code and
%   are hard-coded in those scripts; the h2_opt_* functions do not reproduce
%   them.
%
% -------------------------------------------------------------------------

[ny,n]  = size(C);
if size(V,1) ~= n
    error('Size mismatch: C is %d×%d but V is %d×%d.', ny,n, size(V,1), size(V,2));
end
r = size(V,2);

% ---------- aggregate observation matrix  M ------------------------------
M = zeros(ny,ny);
for k = 1:r
    z = C * V(:,k);            % (n_y × 1) complex
    % lambda_k = eigenvalue associated with V(:,k);
    % weight   = 1 / (-2*real(lambda_k));
    % a = weight * real(z);    b = weight * imag(z);
    a = real(z);
    b = imag(z);
    M = M + a*a.' + b*b.';     % rank-2 update
end

% ---------- principal  r  eigenvectors of  M ------------------------------
if size(M,1) == 1              % single sensor: only one possible direction
    H = 1;
elseif r == 1
    [vec,~] = eigs(M,1,'la');
    H = real(vec);
    H = H / norm(H);
else
    [Vh,~] = eigs(M,r,'la');
    H = real(Vh);
    [H,~] = qr(H,0);           % strict orthonormality
end

% ---------- H2 norm for each column h_j ----------------------------------
H2sq = diag(H.' * M * H);      % squared energy
H2   = sqrt(max(0,H2sq)).';    % row vector 1xr

end
