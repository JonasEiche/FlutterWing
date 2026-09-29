% Geometric Mode Isolation Analysis: Full-State vs Real-Plant Measurability
%
% PURPOSE
%   Investigate whether flutter and residual aeroelastic modes can be
%   geometrically isolated in measurement space using orthogonal
%   projection, i.e., whether a blending vector can reject the residual
%   mode contribution while preserving the flutter mode signal.  This is
%   a pure linear algebra analysis (no controller synthesis) that exposes
%   the fundamental sensor-layout bottleneck for geometric mode isolation.
%
%   Two plants are compared at a single velocity (V_inf = 130 m/s, near
%   flutter onset):
%     1. Full-state plant:  C = I  (56x56) — artificial, sees all states
%     2. Real plant:        C = C_eng (8x56) — 8 distributed accelerometers
%   If geometric isolation works for (1) but fails for (2), the bottleneck
%   is the C-matrix projection from high-dimensional state space into the
%   low-dimensional measurement space, not the modal structure itself.
%
% METHODOLOGY
%   1-2  Build the energy-normalized plant at V_inf = 130 m/s.  Energy
%          normalization uses a diagonal similarity transform T derived from
%          stiffness (Kff), aerodynamic lag, and actuator energy metrics so
%          that pole-vector norms carry physical meaning (modal energy
%          participation).  The transform is constructed explicitly here
%          (not via build_G_RectWing_EnergyNormalized) because the script
%          also needs access to the intermediate aerodynamic matrices.
%
%   3    Eigendecompose A_eng and classify poles:
%           - Flutter poles:   Re(lambda) > -1,  5 < |Im(lambda)| < 40
%           - Residual poles:  Re(lambda) > -10, 5 < |Im(lambda)| < 40,
%                              excluding flutter poles
%         Left/right eigenvectors are balanced (||v|| = ||w||) via
%         geometric mean scaling.
%
%   4    Form real-modal pole vectors via sqrt(2)*real/imag decomposition:
%           v_pos = sqrt(2)*real(v_bal),  v_vel = sqrt(2)*imag(v_bal)
%         Output pole vectors c = C * [v_pos, v_vel] are computed for both
%         the full-state and real plants.  Input pole vectors b = w'*B are
%         computed once (shared B matrix).
%
%   5    SVD analysis: compute sigma_min, sigma_max, and condition number
%         of each pole-vector matrix {flutter 2D, residual 2D, combined 4D}
%         for both plants and for the input side.  sigma_min quantifies the
%         worst-case modal gain; a large condition number indicates that
%         some modal directions are much harder to observe/control.
%
%   6    Orthogonal projection analysis:
%           P_res = c_res * pinv(c_res)   — projector onto residual subspace
%           c_flut_proj = (I - P_res) * c_flut  — flutter after rejection
%         Signal loss = 1 - sigma_min(c_flut_proj) / sigma_min(c_flut).
%         A large signal loss means rejecting the residual mode unavoidably
%         destroys flutter information — geometric isolation is infeasible.
%
%   7    Principal angles between flutter and residual output subspaces:
%           theta_k = acos(svd(Q_flut' * Q_res))
%         where Q are orthonormal bases from QR decomposition.  Large
%         angles (>45 deg) indicate well-separated subspaces; small angles
%         indicate near-collinearity, meaning the modes project onto
%         overlapping directions in measurement space.
%
% POLE-VECTOR SETS (output side, for each plant)
%   flutter 2D      : C * [v_pos_f, v_vel_f]           (ny x 2)
%   residual 2D     : C * [v_pos_r, v_vel_r]           (ny x 2)
%   flutter+res 4D  : C * [v_pos_f, v_vel_f, v_pos_r, v_vel_r]  (ny x 4)
%
% POLE-VECTOR SETS (input side, shared)
%   flutter 2D      : [w_pos_f'*B; w_vel_f'*B]'        (nu x 2)
%   residual 2D     : [w_pos_r'*B; w_vel_r'*B]'        (nu x 2)
%   flutter+res 4D  : all four rows transposed           (nu x 4)
%
%
% OUTPUTS
%   Console tables in 5-8, and Data/Geometric_Isolation.mat (section 9:
%   the values of the paper's Table 4, read by print_Paper_Values.m).
%
% SEE ALSO
%   R01_Modal_Observability      - velocity sweep of sigma_min/gramDet
%                                  (with actuators)
%   R02_Modal_Observability_NoActuators - same without actuators
%   build_G_RectWing             - plant builder (before energy norm)
%   eigenshuffle                 - eigenvalue tracking (not used here;
%                                  single velocity point)
%   gramDet                      - Gram determinant (not used here;
%                                  SVD-based metrics only)

clearvars

%% 1. Configuration
V_inf = 130;
imuIDX = 1:8;
ailIDX = 1:8;
num_modes = 5;
num_poles = 6;

%% 2. Model Setup & Energy Normalization
% (energy metric and similarity transform built inline; see header)
[Structure, Aero] = define_RectWing_Structure_Aero(num_modes, num_poles);
G = build_G_RectWing(V_inf, imuIDX, ailIDX);

Kff = Structure.Kff;
Mff = Structure.Mff;
OMEGA = Structure.OMEGA;
Sfj = Structure.Sfj;
DRe_jx = Structure.DRe_jx;
DIm_jx = Structure.DIm_jx;
rho = Aero.rho;
c_ref = Aero.c_ref;
Q0jj = Aero.Q0jj;
Q1jj = Aero.Q1jj;
QLpjj = Aero.QLpjj;
poles_aero = Aero.poles;
q_bar_ref = 0.5 * rho * 100^2;  % V_inf_ref = 100 as in build_G_RectWing

nx = size(G.A, 1);

% Generalized Aerodynamic Forces from Actuators (for energy normalization)
num_AIL = size(DRe_jx, 2);
num_panels = size(Q0jj, 1);

sumQLpjjBjx = zeros(num_panels, num_AIL);
for i = 1:length(poles_aero)
    sumQLpjjBjx = sumQLpjjBjx + QLpjj(:,:,i)*(DRe_jx - DIm_jx*poles_aero(i)*2/c_ref);
end
Bgx = 0.5*rho*V_inf^2*Sfj*(Q0jj*DRe_jx + sumQLpjjBjx);

sumQLpjjBjx_dot = zeros(num_panels, num_AIL);
for i = 1:length(poles_aero)
    sumQLpjjBjx_dot = sumQLpjjBjx_dot + QLpjj(:,:,i)*DIm_jx;
end
Bgx_dot = 0.5*rho*Sfj*V_inf*(Q0jj*DIm_jx + Q1jj*0.5*c_ref*DRe_jx + sumQLpjjBjx_dot);

% Energy metric & similarity transform
n_aero = 2*num_modes + num_modes*num_poles;  % 40
n_act = nx - n_aero;                          % 16

w0_act = 32*2*pi;  % from define_RectWing_PT2Actuator
Kff_inv = inv(Kff);
T_diag_act_pos = zeros(length(ailIDX), 1);
T_diag_act_vel = zeros(length(ailIDX), 1);
for j = 1:length(ailIDX)
    F_delta     = Bgx(:, ailIDX(j));
    F_delta_dot = Bgx_dot(:, ailIDX(j));
    T_diag_act_pos(j) = sqrt(F_delta' * Kff_inv * F_delta);
    T_diag_act_vel(j) = sqrt(w0_act^2 * (F_delta_dot' * Kff_inv * F_delta_dot));
end
q_diag_act = reshape([T_diag_act_pos'; T_diag_act_vel'], [], 1).^2;

q_diag = [diag(Kff);
          diag(Kff);
          (q_bar_ref*c_ref^2)^2 * repmat(1./diag(Kff), num_poles, 1);
          q_diag_act];

% Normalize to dimensionless energy metric (consistent with build_P_RectWing scaling)
q_diag = q_diag / Kff(1,1);

T = diag(sqrt(q_diag));
T_inv = diag(1./sqrt(q_diag));

A_eng = T * G.A * T_inv;
B_eng = T * G.B;
C_eng = G.C * T_inv;

%% 3. Eigendecomposition & Classification
[V, DD] = eig(A_eng);
W = inv(V);
poles = diag(DD);

flutIDX = find((real(poles) > -1) & (abs(imag(poles)) < 40) & (abs(imag(poles)) > 5));
critIDX = find((real(poles) > -10) & (abs(imag(poles)) < 40) & (abs(imag(poles)) > 5));
resIDX = setdiff(critIDX, flutIDX, 'stable');

fprintf('Flutter poles: ');
for i = 1:length(flutIDX)
    fprintf('%.2f %+.2fi  ', real(poles(flutIDX(i))), imag(poles(flutIDX(i))));
end
fprintf('\nResidual poles: ');
for i = 1:length(resIDX)
    fprintf('%.2f %+.2fi  ', real(poles(resIDX(i))), imag(poles(resIDX(i))));
end
fprintf('\n\n');

%% 3b. Eigenvector Balancing
% Flutter pair
v_R_f = V(:, flutIDX(1));
w_L_f = W(flutIDX(1), :)';
alpha_f = sqrt(norm(w_L_f) / norm(v_R_f));
v_bal_f = v_R_f * alpha_f;
w_bal_f = w_L_f / alpha_f;

% Residual pair
v_R_r = V(:, resIDX(1));
w_L_r = W(resIDX(1), :)';
alpha_r = sqrt(norm(w_L_r) / norm(v_R_r));
v_bal_r = v_R_r * alpha_r;
w_bal_r = w_L_r / alpha_r;

%% 4. Real Modal Form & Pole Vectors
v_pos_f = sqrt(2) * real(v_bal_f);
v_vel_f = sqrt(2) * imag(v_bal_f);
w_pos_f = sqrt(2) * real(w_bal_f);
w_vel_f = sqrt(2) * imag(w_bal_f);

v_pos_r = sqrt(2) * real(v_bal_r);
v_vel_r = sqrt(2) * imag(v_bal_r);
w_pos_r = sqrt(2) * real(w_bal_r);
w_vel_r = sqrt(2) * imag(w_bal_r);

%% 4b. Output Pole Vectors for Two Plants
% Full-state feedback: C = I (artificial, sees all 56 states)
C_full = eye(nx);

% Real plant: C = C_eng (8 accelerometers)
C_real = C_eng;

% Compute output pole vectors for both plants
plants = struct('name', {'Full-state (C=I)', 'Real plant (8 IMUs)'}, ...
                'C', {C_full, C_real});

for p = 1:2
    Cp = plants(p).C;
    plants(p).c_flut = Cp * [v_pos_f, v_vel_f];        % ny x 2
    plants(p).c_res  = Cp * [v_pos_r, v_vel_r];        % ny x 2
    plants(p).c_comb = Cp * [v_pos_f, v_vel_f, v_pos_r, v_vel_r]; % ny x 4
end

% Input pole vectors (shared — same B for both plants)
b_flut = [w_pos_f' * B_eng; w_vel_f' * B_eng];  % 2 x nu
b_res  = [w_pos_r' * B_eng; w_vel_r' * B_eng];  % 2 x nu
b_comb = [b_flut; b_res];                         % 4 x nu

%% 5. SVD sigma_min Analysis
disp('========================================================================')
disp('  SVD ANALYSIS: Full-State vs Real Plant')
disp('========================================================================')
disp(' ')

set_labels = {'flutter 2D', 'residual 2D', 'flutter+res 4D'};

for p = 1:2
    fprintf('--- %s (ny = %d) ---\n', plants(p).name, size(plants(p).C, 1));
    fprintf('%-18s %4s %12s %12s %12s\n', 'Set', 'Dim', 'sigma_min', 'sigma_max', 'Cond');
    disp('------------------------------------------------------------')

    sets = {plants(p).c_flut, plants(p).c_res, plants(p).c_comb};
    for i = 1:3
        s = svd(sets{i});
        fprintf('%-18s %4d %12.4e %12.4e %12.2f\n', ...
            set_labels{i}, size(sets{i}, 2), s(end), s(1), s(1)/s(end));
    end
    fprintf('\n');
end

% Input side (shared)
fprintf('--- Input pole vectors (shared B) ---\n');
fprintf('%-18s %4s %12s %12s %12s\n', 'Set', 'Dim', 'sigma_min', 'sigma_max', 'Cond');
disp('------------------------------------------------------------')
b_sets = {b_flut', b_res', b_comb'};
for i = 1:3
    s = svd(b_sets{i});
    fprintf('%-18s %4d %12.4e %12.4e %12.2f\n', ...
        set_labels{i}, size(b_sets{i}, 2), s(end), s(1), s(1)/s(end));
end
fprintf('\n');

%% 6. Orthogonal Projection Analysis
disp('========================================================================')
disp('  ORTHOGONAL PROJECTION: Signal Loss After Residual Rejection')
disp('========================================================================')
disp(' ')

for p = 1:2
    c_f = plants(p).c_flut;
    c_r = plants(p).c_res;
    ny_p = size(plants(p).C, 1);

    % Projector onto residual output subspace
    P_res = c_r * pinv(c_r);

    % Flutter after rejecting residual
    c_flut_proj = (eye(ny_p) - P_res) * c_f;

    % SVD comparison
    s_orig = svd(c_f);
    s_proj = svd(c_flut_proj);

    signal_loss = 1 - s_proj(end) / s_orig(end);

    fprintf('--- %s (ny = %d) ---\n', plants(p).name, ny_p);
    fprintf('  sigma_min(c_flut):          %12.4e\n', s_orig(end));
    fprintf('  sigma_min(c_flut_proj):     %12.4e\n', s_proj(end));
    fprintf('  Signal loss:                %12.4f  (%.1f%%)\n', signal_loss, signal_loss*100);
    fprintf('  sigma_max(c_flut):          %12.4e\n', s_orig(1));
    fprintf('  sigma_max(c_flut_proj):     %12.4e\n', s_proj(1));
    fprintf('\n');
end

%% 7. Principal Angles Between Flutter and Residual Subspaces
disp('========================================================================')
disp('  PRINCIPAL ANGLES: Flutter vs Residual Output Subspaces')
disp('========================================================================')
disp(' ')

for p = 1:2
    c_f = plants(p).c_flut;
    c_r = plants(p).c_res;

    [Q_f, ~] = qr(c_f, 0);
    [Q_r, ~] = qr(c_r, 0);
    angles = acos(min(svd(Q_f' * Q_r), 1)) * 180/pi;

    fprintf('--- %s ---\n', plants(p).name);
    fprintf('  Principal angles: ');
    fprintf('%.2f°  ', angles);
    fprintf('\n');
    if all(angles > 45)
        fprintf('  => Well-separated subspaces (all angles > 45°)\n');
    else
        fprintf('  => Poorly separated (min angle = %.2f° < 45°)\n', min(angles));
    end
    fprintf('\n');
end

%% 8. Summary & Key Findings
disp('========================================================================')
disp('  SUMMARY')
disp('========================================================================')
disp(' ')

fprintf('%-25s %12s %12s %12s\n', '', 'sigma_min_f', 'sig_loss', 'min_angle');
disp('-------------------------------------------------------------')
for p = 1:2
    c_f = plants(p).c_flut;
    c_r = plants(p).c_res;
    ny_p = size(plants(p).C, 1);

    P_res = c_r * pinv(c_r);
    c_flut_proj = (eye(ny_p) - P_res) * c_f;
    s_orig = svd(c_f);
    s_proj = svd(c_flut_proj);
    signal_loss = 1 - s_proj(end) / s_orig(end);

    [Q_f, ~] = qr(c_f, 0);
    [Q_r, ~] = qr(c_r, 0);
    angles = acos(min(svd(Q_f' * Q_r), 1)) * 180/pi;

    fprintf('%-25s %12.4e %12.1f%% %11.2f°\n', ...
        plants(p).name, s_orig(end), signal_loss*100, min(angles));
end
disp(' ')

disp('KEY FINDING:')
disp('  C projects 56D state space into 8D measurement space.')
disp('  Flutter and residual modes that are well-separated in state space')
disp('  become nearly collinear in measurement space, making geometric')
disp('  mode isolation impractical with the available sensor set.')
disp(' ')

%% 9. Save results (Table 4 of the paper)
results = struct();
results.V_inf = V_inf;
results.flutter_pole = poles(flutIDX(1));
results.residual_pole = poles(resIDX(1));
plantKeys = {'full', 'real'};
for p = 1:2
    c_f = plants(p).c_flut;
    c_r = plants(p).c_res;
    ny_p = size(plants(p).C, 1);

    P_res = c_r * pinv(c_r);
    c_flut_proj = (eye(ny_p) - P_res) * c_f;
    s_orig = svd(c_f);
    s_proj = svd(c_flut_proj);

    [Q_f, ~] = qr(c_f, 0);
    [Q_r, ~] = qr(c_r, 0);
    angles = acos(min(svd(Q_f' * Q_r), 1)) * 180/pi;

    results.(plantKeys{p}) = struct( ...
        'name', plants(p).name, 'ny', ny_p, ...
        'sigminCflut', s_orig(end), ...
        'sigminCflutProj', s_proj(end), ...
        'signalLoss', 1 - s_proj(end)/s_orig(end), ...
        'minPrincipalAngleDeg', min(angles));
end

script_dir = fileparts(mfilename('fullpath'));
data_dir = fullfile(script_dir, 'Data');
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end
save(fullfile(data_dir, 'Geometric_Isolation.mat'), 'results');
disp('Results saved to Data/Geometric_Isolation.mat')
