%% R04: SISO Blending Quadrature Diagnostics on a Canonical 2x2 Oscillator
%
% Self-contained analytical and numerical study proving how SISO blending
% of one oscillatory mode introduces an artificial finite zero whose
% location is governed entirely by the quadrature mismatch angle Delta
% between the blended input and output directions in the modal plane.
%
% Consider a single oscillatory mode with poles at s = sigma +/- j*omega,
% represented by the 2x2 state matrix A = sigma*I + omega*J where
% J = [0 1; -1 0] is the 90-degree rotation generator. The MIMO plant is
%
%       G(s) = C * (sI - A)^{-1} * B
%
% with B in R^{2 x n_u} and C in R^{n_y x 2}. SISO blending selects
% unit vectors k_y in R^{n_y} and k_u in R^{n_u} to form the scalar plant
%
%       g(s) = k_y^T * G(s) * k_u = c^T * (sI - A)^{-1} * b
%
% where c = C^T * k_y and b = B * k_u are the effective output and input
% directions in the 2D state space.
%
% In an orthonormal basis {q_1, q_2} of the modal plane
% (chosen so that A = sigma*I + omega*J in that basis), the directions
% c and b can be parameterised by angles theta_y and theta_u. The SISO
% transfer function then takes the form
%
%       g(s) = [alpha * (s - sigma) + omega * beta] / [(s-sigma)^2 + omega^2]
%
% where the two numerator coefficients are
%
%       alpha = c_hat^T * b_hat       (same-quadrature overlap, = cos(Delta))
%       beta  = c_hat^T * J * b_hat   (cross-quadrature overlap, = sin(Delta))
%
% with c_hat, b_hat being unit vectors in the modal plane and
%
%       Delta = theta_u - theta_y     (quadrature mismatch angle)
%
% The single finite zero of g(s) is therefore
%
%       z(Delta) = sigma - omega * beta/alpha = sigma - omega * tan(Delta)
%
%
% FOUR CANONICAL MISMATCH CASES
% -----------------------------
% The script highlights four geometrically significant configurations:
%
%   Delta = 0 deg:
%       Matched quadratures. alpha = 1, beta = 0. The zero sits at
%       z = sigma, i.e., exactly at the real part of the unstable pole.
%       Maximum direct coupling but the RHP zero severely limits
%       achievable bandwidth.
%
%   Delta = atan2d(sigma, omega):
%       Zero at the origin. The cross-quadrature component exactly
%       cancels the sigma-shift, placing z = 0. This is the boundary
%       between RHP and LHP zeros.
%
%   Delta = 90 deg:
%       Pure cross-quadrature coupling. alpha = 0, so no finite zero
%       exists (zero at infinity). The numerator is a pure constant
%       omega * beta. However, the direct damping authority (alpha) is
%       zero, making stabilisation difficult.
%
%   Delta = 180 deg:
%       Sign-reversed match. alpha = -1, beta = 0. The zero returns to
%       z = sigma (same as Delta = 0). The tan function has period 180 deg.
%
% PHYSICAL INTERPRETATION
% -----------------------
% The mismatch angle Delta quantifies how much the SISO blending rotates
% the actuator direction relative to the sensor direction within the
% oscillatory modal plane. The two quadrature components alpha and beta
% represent fundamentally different coupling mechanisms:
%
%   - alpha (same-quadrature): Direct damping-like coupling. Measures how
%     much the actuator pushes along the direction the sensor observes.
%     This is the "useful" coupling for active damping.
%
%   - beta (cross-quadrature): Stiffness-like coupling. Measures coupling
%     through the 90-degree rotated (quadrature) component. This shifts
%     the natural frequency but does not directly add damping.
%
% The zero z = sigma - omega * beta/alpha captures the trade-off: large
% cross-quadrature coupling (beta) relative to same-quadrature (alpha)
% pulls the zero deep into the LHP, but at the cost of reduced direct
% damping authority. The controller must navigate this trade-off.
%
% H2-OPTIMAL METHODS AND THEIR DEFICIENCY
% ----------------------------------------
% Two H2-optimal SISO blending methods are evaluated:
%
%   1. H2 Blending: Maximises the H2 norm of the SISO transfer
%      function G_flut(s)
%
%   2. H2-opt separate (h2_opt_output_siso + h2_opt_input_siso):
%      Independently maximises the H2 energy coupling from the mode to
%      outputs and from inputs to the mode, then combines.
%
% Both methods maximise modal energy transfer but are agnostic to zero
% placement. The script demonstrates that they land at unfavourable Delta
% values where the induced zero is deep in the RHP, explaining the
% performance gap of the rank-one cases in R06_Locked_Comparison.
%
% SCRIPT STRUCTURE
% ----------------
%   S1: Toy Plant Definition
%       2x2 unstable oscillator A = [2 12; -12 2] with 2 inputs / 2 outputs.
%       The plant has no finite transmission zeros (verified by tzero).
%
%   S2: Orthonormal Modal Basis
%       Eigendecomposition of A to extract sigma, omega, and eigenvectors.
%       Construction of orthonormal basis Q_x via Cholesky factorisation
%       of the Gram matrix (metricOrthonormalBasis), with right-handedness
%       enforced (makeRightHanded) so that beta = sin(Delta) has
%       consistent sign.
%
%   S3: Canonical Mismatch Cases
%       Tabulates the four canonical Delta values and their alpha, beta
%       coefficients. Provides analytical interpretation.
%
%   S4: H2-Optimal SISO Blending Vectors
%       Computes H2 Blending and H2-opt separate blending vectors,
%       projects them into the orthonormal modal basis to determine
%       their Delta values, and verifies the zero formula against tzero.
%
%   S5: Verify Zero Formula (compact sweep)
%       Sweeps Delta over a fine grid and verifies
%       z(Delta) = sigma - omega * tan(Delta) to machine precision
%       by comparing against tzero of the constructed SISO plant.
%       Also verifies the corollary: Delta = 0 => z = sigma.
%
%   S6: Systune Sweep over Delta
%       Sweeps Delta over [0, 175] deg in 2-degree steps with the H2
%       methods injected into the grid. At each Delta, constructs the
%       SISO plant, synthesises a 2nd-order controller via systune with
%       gain margin >= 6 dB, phase margin >= 45 deg (hard constraint)
%       and output-to-input gain < 10 (soft goal). Uses 3 random starts
%       with a fixed seed for reproducibility. Records soft and hard
%       goal values to characterise the performance landscape.
%
%   S7: Summary
%       Reports the best sweep point and the H2 methods' performance.
%       Prints the compact interpretation linking alpha, beta, z, and
%       controller performance.
%
%   S8: Save Results
%       Saves all sweep data and H2 method results to
%       Data/RHP_zeros_analysis.mat for use by fig04_fig05_SISO_Zero_Analysis.m.
%
% OUTPUTS
% -------
%   Data/RHP_zeros_analysis.mat containing:
%       - A, B, C: toy plant matrices
%       - sigma, omega: mode parameters (sigma=2, omega=12)
%       - Qx: orthonormal modal basis
%       - Delta_deg_sweep: sweep grid in degrees
%       - sweep_alpha, sweep_beta: numerator coefficients at each Delta
%       - sweep_z_formula, sweep_z_numeric: zero location (formula vs tzero)
%       - sweep_softGoal, sweep_hardGoal: systune performance metrics
%       - h2_names: {'H2 Blending', 'H2-opt separate'}
%       - h2_Delta_deg: raw Delta of H2 methods
%       - h2_Delta_display: Delta mapped to [0, 180) for display
%       - h2_sweep_idx: indices into sweep grid for H2 methods
%
%   Figures generated by fig04_fig05_SISO_Zero_Analysis.m (separate script):
%       - fig04_siso_numerator_geometry (alpha, beta, z vs Delta)
%       - fig05_siso_performance_vs_delta (soft goal vs Delta, z, alpha)
%
%
% DEPENDENCIES
% ------------
%   Toolboxes: Control System Toolbox (systune, TuningGoal, tunableSS, AnalysisPoint)
%   Functions: h2_opt_output_siso, h2_opt_input_siso
%
% SEE ALSO
%   R05_Quadrature_Mismatch.m  - Same phenomenon on full FlutterWing plant
%   fig04_fig05_SISO_Zero_Analysis.m  - Paper-ready figures from R04 results
%   R06_Locked_Comparison.m - 2x2 MIMO baseline on this oscillator
%
% -------------------------------------------------------------------------

clear; clc; close all;

%% S1: Toy Plant Definition
%
% 2x2 unstable oscillator representing a single flutter mode with two
% inputs and two outputs before loop closure.
% A = sigma*I + omega*J, where J is the 90 degree rotation generator.

A = [  2,  12;
     -12,   2];

B = [  1,   3;
      -3,   1];

C = [  1,   2;
      -2,   1];

G = ss(A, B, C, 0);
[ny, nu] = size(G);

% Verify that the original 2x2 plant has no finite zeros.
zG = tzero(G);
assert(isempty(zG), 'The original 2x2 toy plant should have no finite transmission zeros.');
fprintf('Original 2x2 plant has no finite zeros.\n');

%% S2: Orthonormal Modal Basis
%
% The angle Delta only has physical meaning if measured in a proper 2D modal
% coordinate system. First identify the 2D state-space modal plane, then
% build an orthonormal basis. Only in that basis do sensor and actuator
% directions become true geometric directions with a meaningful relative angle.

[V, Lambda] = eig(A);
lam = diag(Lambda);

% Pick the eigenvalue with positive imaginary part.
idx = find(imag(lam) > 0, 1, 'first');
assert(~isempty(idx), 'Expected one eigenvalue with positive imaginary part.');

lambda = lam(idx);
sigma  = real(lambda);
omega  = imag(lambda);
v      = V(:, idx);

% Left eigenvector (needed for h2_opt_input_siso).
W = inv(V)';
w = W(:, idx);

% Metric for orthonormality of the modal plane.
% Euclidean for this toy oscillator; replace by energy-consistent metric in general.
Mx = eye(2);

% Build orthonormal basis Qx = [q1 q2] of the 2D modal plane.
Qx = metricOrthonormalBasis([real(v), imag(v)], Mx);
Qx = makeRightHanded(Qx);

% Verify canonical oscillator form [sigma omega; -omega sigma].
A_modal = Qx' * A * Qx;
A_modal_expected = [sigma, omega; -omega, sigma];

fprintf('\nA expressed in the orthonormal modal basis:\n');
disp(A_modal);

assert(norm(A_modal - A_modal_expected, 'fro') < 1e-10, ...
    'Modal basis does not reproduce the canonical oscillator form.');

% Generator of 90 degree rotation in the orthonormal modal plane.
J = [0, 1; -1, 0];

%% S3: Canonical Mismatch Cases
%
% Delta = 0 deg:    matched quadratures, z = sigma (RHP)
% Delta = atan2d(sigma,omega):  zero at the origin
% Delta = 90 deg:   alpha = 0, zero at infinity (singular)
% Delta = 180 deg:  sign-reversed match, z = sigma again

Delta_zero_deg = atan2d(sigma, omega);

canonical_names = { ...
    'Matched quadratures', ...
    'Zero at origin', ...
    'Zero at infinity', ...
    'Sign-reversed match'};

canonical_Delta_deg = [0, Delta_zero_deg, 90, 180];

fprintf('\nCanonical cases:\n');
fprintf('%-24s %12s %12s %12s\n', 'Case', 'Delta [deg]', 'alpha', 'beta');
fprintf('%s\n', repmat('-', 1, 64));
for i = 1:numel(canonical_names)
    D = canonical_Delta_deg(i);
    a = cosd(D);
    b = sind(D);
    fprintf('%-24s %12.3f %12.3f %12.3f\n', canonical_names{i}, D, a, b);
end
fprintf('\nInterpretation:\n');
fprintf('  Delta = 0 deg      -> z = sigma = %.3f (RHP)\n', sigma);
fprintf('  Delta = %.3f deg -> z = 0\n', Delta_zero_deg);
fprintf('  Delta = 90 deg     -> no finite zero (zero at infinity)\n');
fprintf('  Delta = 180 deg    -> z = sigma = %.3f again\n', sigma);

%% S4: H2-Optimal SISO Blending Vectors
% Method 1: H2 Optimal Blending
% The code for the calculation of the H_2 optimal blending vectors is proprietary, hence hardcoded here: 
ky_h2 = [0.0825; 0.9966];
ku_h2 = [-1;0];

% Method 2: H2-optimal separate input/output (SISO)
[~, h] = h2_opt_output_siso(C, v);
[~, d] = h2_opt_input_siso(B, w);
ky_h2sep = -h;
ku_h2sep = d;

% Project blending vectors into orthonormal modal basis to find Delta
h2_names  = {'H2 Blending', 'H2-opt separate'};
h2_ky_all = {ky_h2, ky_h2sep};
h2_ku_all = {ku_h2, ku_h2sep};
h2_Delta_deg = nan(2, 1);
h2_alpha     = nan(2, 1);
h2_beta      = nan(2, 1);
h2_z_formula = nan(2, 1);
h2_z_tzero   = nan(2, 1);

fprintf('\n=== H2-Optimal SISO Blending Vectors ===\n');
fprintf('%-18s %12s %12s %12s %12s %12s\n', ...
    'Method', 'Delta [deg]', 'alpha', 'beta', 'z_formula', 'z_tzero');
fprintf('%s\n', repmat('-', 1, 82));

for i = 1:2
    ky_i = h2_ky_all{i};
    ku_i = h2_ku_all{i};

    % Project into modal basis
    c_state = C' * ky_i;
    b_state = B  * ku_i;
    c_modal = Qx' * Mx * c_state;
    b_modal = Qx' * Mx * b_state;

    theta_y_i = atan2(c_modal(2), c_modal(1));
    theta_u_i = atan2(b_modal(2), b_modal(1));
    Delta_i   = wrapToPiLocal(theta_u_i - theta_y_i);

    % Numerator coefficients from modal coordinates (unit-normalized)
    c_unit = c_modal / norm(c_modal);
    b_unit = b_modal / norm(b_modal);
    al = c_unit' * b_unit;
    be = c_unit' * J * b_unit;

    % Zero from formula
    if abs(al) > 1e-12
        zf = sigma - omega * be / al;
    else
        zf = sign(be) * Inf;
    end

    % Verify against tzero
    zt = real(tzero(ky_i' * G * ku_i));
    assert(abs(zf - zt) < 1e-6, 'Zero formula mismatch for %s', h2_names{i});

    h2_Delta_deg(i) = rad2deg(Delta_i);
    h2_alpha(i)     = al;
    h2_beta(i)      = be;
    h2_z_formula(i) = zf;
    h2_z_tzero(i)   = zt;

    fprintf('%-18s %12.1f %12.4f %12.4f %12.4f %12.4f\n', ...
        h2_names{i}, h2_Delta_deg(i), al, be, zf, zt);
end

%% S5: Verify Zero Formula (compact sweep)
%
% Verify z(Delta) = sigma - omega*tan(Delta) to machine precision.

Delta_verify = linspace(-1.2, 1.2, 51);
z_law = sigma - omega * tan(Delta_verify);
z_num = nan(size(Delta_verify));
theta_y_fixed = 0;

for i = 1:numel(Delta_verify)
    c_modal_i = [cos(theta_y_fixed); sin(theta_y_fixed)];
    b_modal_i = [cos(theta_y_fixed + Delta_verify(i)); sin(theta_y_fixed + Delta_verify(i))];

    c_state_i = Qx * c_modal_i;
    b_state_i = Qx * b_modal_i;

    ky_i = C' \ c_state_i;
    ku_i = B  \ b_state_i;

    z = tzero(ky_i' * G * ku_i);
    z_num(i) = real(z(1));
end

sweep_err = max(abs(z_num - z_law));
assert(sweep_err < 1e-10, 'Quadrature sweep formula mismatch exceeds tolerance');
fprintf('\nz(Delta) = sigma - omega*tan(Delta) verified over %d points (max err: %.1e)\n', ...
    numel(Delta_verify), sweep_err);

% Matched-quadrature corollary: Delta=0 => z=sigma
c_state_m = Qx * [1; 0];
b_state_m = Qx * [1; 0];
ky_m = C' \ c_state_m;
ku_m = B  \ b_state_m;
assert(abs(real(tzero(ky_m' * G * ku_m)) - sigma) < 1e-10, ...
    'Matched quadrature should give z=sigma');
fprintf('Corollary: Delta=0 => z = sigma = %.1f (verified)\n\n', sigma);

%% S6: Systune Sweep over Delta
%
% Sweep Delta over [0, 175] deg with 5-degree steps, plus the H2 methods'
% display-Deltas injected into the grid. This ensures the H2 methods go
% through the identical systune pipeline (same blending-vector construction,
% same random seed) for a fair comparison.
%
% H2 methods have negative Delta; we map to [0,180) via Delta+180.
% tan has period 180, so the induced zero is the same.

fprintf('=== Systune sweep over quadrature mismatch ===\n');

% Map H2 Delta into [0, 180) for inclusion in sweep grid
h2_Delta_display = h2_Delta_deg;
h2_Delta_display(h2_Delta_display < 0) = h2_Delta_display(h2_Delta_display < 0) + 180;

% Build sweep grid: regular 5-degree steps plus H2 display-Deltas
Delta_deg_sweep = unique(sort([0:2:175, h2_Delta_display(:)']));
nSweep = numel(Delta_deg_sweep);

% Identify which sweep indices correspond to H2 methods
[~, h2_sweep_idx] = ismember(h2_Delta_display, Delta_deg_sweep);

sweep_alpha     = nan(nSweep, 1);
sweep_beta      = nan(nSweep, 1);
sweep_z_formula = nan(nSweep, 1);
sweep_z_numeric = nan(nSweep, 1);
sweep_softGoal  = nan(nSweep, 1);
sweep_hardGoal  = nan(nSweep, 1);

nO = 2;  % 2nd order controller
ReqMarg = TuningGoal.Margins('ud', 6, 45);
ReqAtten = TuningGoal.Gain({'ym'}, {'ud'}, 10);

% Open a local pool if the Parallel Computing Toolbox is installed; without it
% systune ignores 'UseParallel' and evaluates the random starts serially (slower).
if ~isempty(ver('parallel')) && isempty(gcp('nocreate')), parpool('Processes'); end
opt = systuneOptions('RandomStart', 3, 'UseParallel', true, 'Display', 'off');

theta_y_syn = 0;  % fixed output angle
rngSeed = 20260305;     % fixed seed for reproducibility

fprintf('%-12s %12s %12s %12s %12s %12s %8s\n', ...
    'Delta [deg]', 'alpha', 'beta', 'z_formula', 'z_numeric', 'SoftGoal', 'Status');
fprintf('%s\n', repmat('-', 1, 84));

for k = 1:nSweep
    Ddeg = Delta_deg_sweep(k);
    D    = deg2rad(Ddeg);

    theta_u_k = theta_y_syn + D;

    % Unit directions in orthonormal modal basis
    c_modal = [cos(theta_y_syn); sin(theta_y_syn)];
    b_modal = [cos(theta_u_k);   sin(theta_u_k)];

    % Numerator coefficients
    sweep_alpha(k) = c_modal' * b_modal;
    sweep_beta(k)  = c_modal' * J * b_modal;

    % Exact zero formula
    if abs(sweep_alpha(k)) > 1e-12
        sweep_z_formula(k) = sigma - omega * sweep_beta(k) / sweep_alpha(k);
    else
        sweep_z_formula(k) = sign(sweep_beta(k)) * Inf;
    end

    % Convert to state coordinates and solve for blending vectors
    c_state = Qx * c_modal;
    b_state = Qx * b_modal;
    ky_k = C' \ c_state;
    ku_k = B  \ b_state;

    % Numeric zero verification
    g_k = minreal(ky_k' * G * ku_k, 1e-9);
    z_k = tzero(g_k);
    if isempty(z_k)
        sweep_z_numeric(k) = sign(sweep_beta(k)) * Inf;
    else
        sweep_z_numeric(k) = real(z_k(1));
    end

    % Systune with correct analysis points at the 2x2 plant I/O
    tuneCont_k = ku_k * tunableSS('cont', nO, 1, 1) * ky_k';
    tuneCL_k = feedback(G, AnalysisPoint('ud', nu) * tuneCont_k * AnalysisPoint('ym', ny));
    rng(rngSeed);
    [~, fSoft_k, gHard_k] = systune(tuneCL_k, ReqAtten, ReqMarg, opt);

    sweep_softGoal(k) = fSoft_k;
    sweep_hardGoal(k) = gHard_k;

    status_k = "OK";
    if gHard_k > 1, status_k = "FAIL"; end

    % Mark H2 methods in the table
    label = '';
    if k == h2_sweep_idx(1), label = '  <-- H2 Blending'; end
    if k == h2_sweep_idx(2), label = '  <-- H2-opt sep.'; end
    fprintf('%12.1f %12.4f %12.4f %12.4f %12.4f %12.4f %8s%s\n', ...
        Ddeg, sweep_alpha(k), sweep_beta(k), sweep_z_formula(k), ...
        sweep_z_numeric(k), fSoft_k, status_k, label);
end

% Verify zero formula over finite-zero sweep points
idxFinite = isfinite(sweep_z_formula);
maxZeroErr = max(abs(sweep_z_formula(idxFinite) - sweep_z_numeric(idxFinite)));
fprintf('\nMaximum |z_formula - z_numeric| over finite-zero cases: %.3e\n\n', maxZeroErr);

% Extract H2 method results from the sweep
h2_softGoal = sweep_softGoal(h2_sweep_idx);
h2_hardGoal = sweep_hardGoal(h2_sweep_idx);

%% S7: Summary

[bestSoft, bestIdx] = min(sweep_softGoal);
bestDelta = Delta_deg_sweep(bestIdx);

fprintf('=== Summary ===\n');
fprintf('Best SISO sweep:  Delta = %+.0f deg, alpha = %.4f, beta = %.4f, z = %.4f, SoftGoal = %.4f\n', ...
    bestDelta, sweep_alpha(bestIdx), sweep_beta(bestIdx), sweep_z_formula(bestIdx), bestSoft);
for i = 1:2
    idx_i = h2_sweep_idx(i);
    fprintf('%-18s display-Delta = %+.1f deg, z = %.4f, SoftGoal = %.4f\n', ...
        h2_names{i}, h2_Delta_display(i), sweep_z_formula(idx_i), sweep_softGoal(idx_i));
end


% See fig04_fig05_SISO_Zero_Analysis.m for plots

scriptDir = fileparts(mfilename('fullpath'));

%% S8: Save Results

dataDir = fullfile(scriptDir, 'Data');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

results.A = A;
results.B = B;
results.C = C;
results.sigma = sigma;
results.omega = omega;
results.Qx = Qx;
results.Delta_deg_sweep = Delta_deg_sweep;
results.sweep_alpha = sweep_alpha;
results.sweep_beta = sweep_beta;
results.sweep_z_formula = sweep_z_formula;
results.sweep_z_numeric = sweep_z_numeric;
results.sweep_softGoal = sweep_softGoal;
results.sweep_hardGoal = sweep_hardGoal;
results.h2_names = h2_names;
results.h2_Delta_deg = h2_Delta_deg;
results.h2_Delta_display = h2_Delta_display;
results.h2_sweep_idx = h2_sweep_idx;

save(fullfile(dataDir, 'RHP_zeros_analysis.mat'), 'results');
fprintf('Results saved to Data/RHP_zeros_analysis.mat\n');

%% ========================================================================
%  Helper Functions
%  ========================================================================

function Q = metricOrthonormalBasis(X, M)
% Build orthonormal basis Q spanning columns of X such that Q'*M*Q = I.
% For this toy example M = I; generalize to energy-consistent metric.
    G = X' * M * X;
    R = chol(G);
    Q = X / R;
end

function Q = makeRightHanded(Q)
% Ensure consistent orientation of the 2D orthonormal basis.
% The sign of beta = c^T J b depends on basis orientation.
    if det(Q) < 0
        Q(:,2) = -Q(:,2);
    end
end

function a = wrapToPiLocal(a)
    a = mod(a + pi, 2*pi) - pi;
end
