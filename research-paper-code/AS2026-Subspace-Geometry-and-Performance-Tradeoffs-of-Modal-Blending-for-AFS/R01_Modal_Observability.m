% Modal Observability/Controllability analysis using energy-normalized
% plant models WITH actuator dynamics across the velocity envelope.
%
% PURPOSE
%   Quantify how well the flutter and residual aeroelastic modes can be
%   observed (via IMU accelerometers) and controlled (via actuated control
%   surfaces) when the full actuator dynamics are included in the plant.
%   This complements R02, which strips actuator dynamics to isolate the
%   pure structural contribution.  Comparing both reveals whether actuator
%   bandwidth or coupling degrades modal observability/controllability.
%
%   Energy normalization (build_G_RectWing_EnergyNormalized) ensures that
%   the state basis weights strain- and kinetic-energy equally, giving
%   pole-vector norms a direct physical interpretation as modal energy
%   participation.
%
% METHODOLOGY
%   1. At each freestream velocity V_inf in [40, 160] m/s (32 points) the
%      open-loop plant G(s) is built and eigendecomposed.  Eigenvectors are
%      tracked across velocities using eigenshuffle to maintain consistent
%      pole identity despite mode veering.
%
%   2. Poles are classified at a reference velocity (~130 m/s, near flutter
%      onset) into:
%        - Flutter poles:   Re(lambda) > -1,  5 < |Im(lambda)| < 40
%        - Residual poles:  Re(lambda) > -10, 5 < |Im(lambda)| < 40,
%                           excluding flutter poles
%
%   3. For each conjugate pole pair (flutter, residual) the right
%      eigenvector v_R and left eigenvector w_L are balanced so that
%      ||v|| = ||w|| (geometric mean scaling), then decomposed into
%      real-form position/velocity components:
%        v_pos = sqrt(2) * real(v_bal),   v_vel = sqrt(2) * imag(v_bal)
%        w_pos = sqrt(2) * real(w_bal),   w_vel = sqrt(2) * imag(w_bal)
%      The sqrt(2) factor preserves the energy norm under the real
%      transformation.
%
%   4. Output pole vectors c = C * v and input pole vectors b = (w' * B)'
%      are formed.  These map the internal modal participation into the
%      physical sensor/actuator space.
%
% POLE-VECTOR SETS (9 total, covering 1D magnitude and 2D/4D MIMO)
%   1D sets (magnitude only, lose phase/direction information):
%     |flutter| 1D       - |C * v_bal_f|           / |(w_bal_f' * B)'|
%     |flutter pos| 1D   - |c_pos_f|               / |b_pos_f|
%     |flutter vel| 1D   - |c_vel_f|               / |b_vel_f|
%     |residual| 1D      - |C * v_bal_r|           / |(w_bal_r' * B)'|
%     |residual pos| 1D  - |c_pos_r|               / |b_pos_r|
%     |residual vel| 1D  - |c_vel_r|               / |b_vel_r|
%   2D sets (retain directional information within one mode):
%     flutter 2D         - [c_pos_f, c_vel_f]      / [b_pos_f, b_vel_f]
%     residual 2D        - [c_pos_r, c_vel_r]      / [b_pos_r, b_vel_r]
%   4D set (joint flutter + residual subspace):
%     flutter+res 4D     - [c_pos_f, c_vel_f, c_pos_r, c_vel_r] / same for B
%
%   The 1D sets test whether scalar (SISO) blending can capture sufficient
%   modal energy.  The 2D sets test independent observability within a
%   single mode's subspace.  The 4D set tests whether flutter and residual
%   modes can be simultaneously and independently observed/controlled.
%
% METRICS (computed for each set at each velocity)
%   gramDet (raw)        - Gram determinant det(M'*M); absolute volume of
%                          the parallelepiped spanned by pole-vector columns.
%   gramDet (normalized) - Normalized Gram determinant; 1 = perfectly
%                          orthogonal columns, 0 = linearly dependent.
%                          Measures independent measurability/controllability.
%   sigma_min            - Smallest singular value of the pole-vector matrix;
%                          worst-case modal gain in any direction.  A small
%                          sigma_min means at least one modal direction is
%                          poorly observable/controllable.
%   condition number     - sigma_max / sigma_min; quantifies how uniform
%                          the modal coupling is across directions.
%
% OUTPUTS
%   Data/Modal_Obsrv_Contr_data.mat containing:
%     V_inf_arr, set_names,
%     GramRaw_C, GramRaw_B    - raw Gram determinant (output/input)
%     GramNorm_C, GramNorm_B  - normalized Gram determinant
%     SigMin_C, SigMin_B      - minimum singular values
%     Cond_C, Cond_B          - condition numbers
%
%   Read by fig03_Sigmin_Velocity_Sweep.m (Fig. 3), fig07_Condition_Numbers.m
%   (Fig. 7) and print_Paper_Values.m.
%
%
% SEE ALSO
%   R02_Modal_Observability_NoActuators  - same analysis without actuators
%   build_G_RectWing_EnergyNormalized     - plant builder with energy norm
%   eigenshuffle                          - eigenvalue tracking
%   gramDet                               - Gram determinant computation
%   fig03_Sigmin_Velocity_Sweep, fig07_Condition_Numbers  - plotting scripts

clearvars

%% 1. Configuration
textwidth = 16.5730; % cm
figwidth = 0.47*textwidth; % cm
figheight = 0.7*figwidth;
docFontSize = 9; % pt
stdLineWidth = 1.2;

V_inf_arr = linspace(40,160,32);
nv = length(V_inf_arr);
imuIDX = 1:8;
ailIDX = 1:8;
ny = length(imuIDX);
nu = length(ailIDX);

set_names = {'|flutter| 1D', '|flutter pos| 1D', '|flutter vel| 1D', ...
             '|residual| 1D', '|residual pos| 1D', '|residual vel| 1D', ...
             'flutter 2D', 'residual 2D', 'flutter+res 4D'};
nSets = length(set_names);

%% 2. Pass 1 - Build Energy-Normalized Models & Eigendecomposition
G_ref = build_G_RectWing_EnergyNormalized(100, imuIDX, ailIDX);
nx = length(G_ref.A);

G_arr = ss(zeros(ny,nu,1,nv));
PHIV = zeros(nx,nx,nv);
PHIW = zeros(nx,nx,nv);
EV   = zeros(nx,nv);

for iv = 1:nv
    G = build_G_RectWing_EnergyNormalized(V_inf_arr(iv), imuIDX, ailIDX);
    G_arr(:,:,1,iv) = G;
    [V,DD] = eig(G.A);
    PHIV(:,:,iv) = V;
    EV(:,iv)     = diag(DD);
    if iv >= 2
        idx = eigenshuffle(EV(:,iv-1), PHIV(:,:,iv-1), EV(:,iv), PHIV(:,:,iv));
    else
        idx = 1:nx;
    end
    PHIV(:,:,iv) = PHIV(:,idx,iv);
    PHIW(:,:,iv) = inv(PHIV(:,:,iv));
    EV(:,iv)     = EV(idx,iv);
end

% Sort by real part at final velocity
[~, idx] = sort(EV(:,end), 'descend', 'ComparisonMethod', 'real');
EV   = EV(idx,:);
PHIV = PHIV(:,idx,:);
PHIW = PHIW(idx,:,:);

% Classify flutter/residual poles at reference velocity (~130 m/s)
[~, refIDX] = min(abs(V_inf_arr - 130));
EVref = EV(:,refIDX);
flutIDX = find((real(EVref) > -1) & (abs(imag(EVref)) < 40) & (abs(imag(EVref)) > 5));
critIDX = find((real(EVref) > -10) & (abs(imag(EVref)) < 40) & (abs(imag(EVref)) > 5));
resIDX = setdiff(critIDX, flutIDX, 'stable');

fprintf('Reference V_inf = %.1f m/s (index %d)\n', V_inf_arr(refIDX), refIDX);
fprintf('Flutter poles: ');
for i = 1:length(flutIDX)
    fprintf('%.2f %+.2fi  ', real(EVref(flutIDX(i))), imag(EVref(flutIDX(i))));
end
fprintf('\nResidual poles: ');
for i = 1:length(resIDX)
    fprintf('%.2f %+.2fi  ', real(EVref(resIDX(i))), imag(EVref(resIDX(i))));
end
fprintf('\n\n');

%% 3. Pass 2 - Pole-Vector Metrics Sweep
GramRaw_C = NaN(nSets, nv);  SigMin_C = NaN(nSets, nv);  Cond_C = NaN(nSets, nv);
GramRaw_B = NaN(nSets, nv);  SigMin_B = NaN(nSets, nv);  Cond_B = NaN(nSets, nv);
GramNorm_C = NaN(nSets, nv);
GramNorm_B = NaN(nSets, nv);

for iv = 1:nv
    if isempty(flutIDX) || isempty(resIDX)
        continue
    end

    C_iv = G_arr(:,:,1,iv).C;
    B_iv = G_arr(:,:,1,iv).B;

    % Flutter eigenvectors (first of conjugate pair)
    v_R_f = PHIV(:, flutIDX(1), iv);
    w_L_f = PHIW(flutIDX(1), :, iv)';

    % Residual eigenvectors
    v_R_r = PHIV(:, resIDX(1), iv);
    w_L_r = PHIW(resIDX(1), :, iv)';

    % Balance: equal norm for left/right
    alpha_f = sqrt(norm(w_L_f) / norm(v_R_f));
    v_bal_f = v_R_f * alpha_f;
    w_bal_f = w_L_f / alpha_f;

    alpha_r = sqrt(norm(w_L_r) / norm(v_R_r));
    v_bal_r = v_R_r * alpha_r;
    w_bal_r = w_L_r / alpha_r;

    % Real modal form: sqrt(2) * real/imag preserves energy norm
    v_pos_f = sqrt(2) * real(v_bal_f);
    v_vel_f = sqrt(2) * imag(v_bal_f);
    w_pos_f = sqrt(2) * real(w_bal_f);
    w_vel_f = sqrt(2) * imag(w_bal_f);

    v_pos_r = sqrt(2) * real(v_bal_r);
    v_vel_r = sqrt(2) * imag(v_bal_r);
    w_pos_r = sqrt(2) * real(w_bal_r);
    w_vel_r = sqrt(2) * imag(w_bal_r);

    % Output pole vectors: c = C * v
    c_pos_f = C_iv * v_pos_f;
    c_vel_f = C_iv * v_vel_f;
    c_pos_r = C_iv * v_pos_r;
    c_vel_r = C_iv * v_vel_r;

    % Input pole vectors: b = (w' * B)'
    b_pos_f = (w_pos_f' * B_iv)';
    b_vel_f = (w_vel_f' * B_iv)';
    b_pos_r = (w_pos_r' * B_iv)';
    b_vel_r = (w_vel_r' * B_iv)';

    % 9 output sets (ny x dim)
    C_sets = {abs(C_iv * v_bal_f), ...
              abs(c_pos_f), ...
              abs(c_vel_f), ...
              abs(C_iv * v_bal_r), ...
              abs(c_pos_r), ...
              abs(c_vel_r), ...
              [c_pos_f, c_vel_f], ...
              [c_pos_r, c_vel_r], ...
              [c_pos_f, c_vel_f, c_pos_r, c_vel_r]};

    % 9 input sets (nu x dim)
    B_sets = {abs((w_bal_f' * B_iv)'), ...
              abs(b_pos_f), ...
              abs(b_vel_f), ...
              abs((w_bal_r' * B_iv)'), ...
              abs(b_pos_r), ...
              abs(b_vel_r), ...
              [b_pos_f, b_vel_f], ...
              [b_pos_r, b_vel_r], ...
              [b_pos_f, b_vel_f, b_pos_r, b_vel_r]};

    % Compute metrics per set
    for i = 1:nSets
        GramRaw_C(i,iv) = gramDet(C_sets{i});
        GramNorm_C(i,iv) = gramDet(C_sets{i}, 'normalize', true);
        s_c = svd(C_sets{i});
        SigMin_C(i,iv) = s_c(end);
        Cond_C(i,iv) = s_c(1) / s_c(end);

        GramRaw_B(i,iv) = gramDet(B_sets{i});
        GramNorm_B(i,iv) = gramDet(B_sets{i}, 'normalize', true);
        s_b = svd(B_sets{i});
        SigMin_B(i,iv) = s_b(end);
        Cond_B(i,iv) = s_b(1) / s_b(end);
    end
end

%% 4. Save Data
script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
data_dir = fullfile(script_dir, 'Data');
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end
save(fullfile(data_dir, 'Modal_Obsrv_Contr_data.mat'), ...
    'V_inf_arr', 'set_names', ...
    'GramRaw_C', 'GramRaw_B', 'GramNorm_C', 'GramNorm_B', ...
    'SigMin_C', 'SigMin_B', 'Cond_C', 'Cond_B');

% Plots: fig03_Sigmin_Velocity_Sweep.m, fig07_Condition_Numbers.m

fprintf('Data saved to %s\n', fullfile(data_dir, 'Modal_Obsrv_Contr_data.mat'));
