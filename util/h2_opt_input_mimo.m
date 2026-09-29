function [D,H2] = h2_opt_input_mimo(B,W)
%H2_OPT_INPUT_MIMO  Real r-dimensional input projector maximising the energy injected into r modes.
%
%   [D,H2] = h2_opt_input_mimo(B,W)
%   NOTE: the projector D is the FIRST output; the SISO variants
%   h2_opt_input_siso / h2_opt_output_siso return the direction SECOND.
%   Requires r = size(W,2) <= n_u (eigs(M,r) errors otherwise).
%
%   Inputs
%   ------
%     B : (n  x n_u)  real      — plant input matrix
%     W : (n  x r)    complex   — left eigenvectors of all modes of interest
%                                 * include both members of each complex pair
%                                 * include one real w for every real pole
%
%   Outputs
%   -------
%     D  : (n_u x r)  real      — optimal projector, D.'*D = eye(r)
%     H2 : (1   x r)  real      — H2 norm of each individual column d_j,
%                                 i.e.  |G_{d_j}|_2   (not squared)
%
%   Theory
%   ------
%   Aggregate symmetric matrix
%       M = sum_k  (Re(z_k) Re(z_k)^T + Im(z_k) Im(z_k)^T),
%       z_k = B^T w_k            (rank ≤ 2r)
%
%   Optimal D  = principal r eigenvectors of M.
%   Energy for each column d_j:
%       E_j = d_j^T M d_j   (squared H2 norm without sigma-weighting)
%       H2_j = sqrt(E_j).
%
%   If you want the *true* H2 energy   |w_k^H B d|_2 / (−2sigma_k),
%   multiply each rank-2 term by 1/(−2sigma_k) inside the loop.
%
%
% -------------------------------------------------------------------------

[n,nu] = size(B);
if size(W,1) ~= n
    error('Size mismatch: B is %dx%d but W is %dx%d.', n,nu, size(W,1), size(W,2));
end
r = size(W,2);

% ---------- build aggregate symmetric matrix M ---------------------------
M = zeros(nu,nu);
for k = 1:r
    z = B.' * W(:,k);              % (n_u x 1) complex
    % ---- OPTIONAL true H2 weighting ----
    % If you know the pole’s real part sigma_k < 0, multiply the next two
    % lines by 1/(-2*real(lambda_k)).
    a = real(z);
    b = imag(z);
    M = M + a*a.' + b*b.';         % rank-2 update
end

% ---------- principal r eigenvectors of  M -------------------------------
if size(M,1) == 1              % single actuator: only one possible direction
    D = 1;
elseif r == 1
    [vec,~] = eigs(M,1,'la');
    D = real(vec);
    D = D / norm(D);
else
    [V,~] = eigs(M,r,'la');        % columns orthonormal up to round-off
    D = real(V);
    [D,~] = qr(D,0);               % enforce strict orthonormality
end

% ---------- H2 norm for each column  -------------------------------------
H2sq = diag(D.' * M * D);          % squared “energy” of each d_j
H2   = sqrt(max(0,H2sq)).';        % row vector 1xr, guard round-off

end





