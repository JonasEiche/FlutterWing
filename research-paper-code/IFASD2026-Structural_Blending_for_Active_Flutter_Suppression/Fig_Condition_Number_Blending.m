%% Condition Number of Structural Blending Matrix Blocks
%  Validates that pinv(PHIzf(:,1:2)) and pinv(PHIfx(1:2,:)) are
%  well-posed for the available sensor/actuator array.
%  Both matrices are V-independent (q_bar dropped from PHIfx).
clearvars
textwidth = 15.98; % cm (IFASD 2026 template: 455.24pt)
figwidth = 0.47*textwidth; % cm
figheight = 0.5*figwidth;
docFontSize = 9; % pt

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end

% ---- Load structural/aero parameters -----------------------------------
num_modes = 5;
num_poles = 6;

[Structure, Aero] = define_RectWing_Structure_Aero(num_modes, num_poles);
Sfj    = Structure.Sfj;
PHIgf  = Structure.PHIgf;
DRe_jx = Structure.DRe_jx;
DIm_jx = Structure.DIm_jx;
PHIzg  = Structure.PHIzg;

c_ref  = Aero.c_ref;
poles  = Aero.poles;
Q0jj   = Aero.Q0jj;
QLpjj  = Aero.QLpjj;

num_panels = size(Q0jj,1);
num_AIL = size(DRe_jx,2);

% ---- Output blending matrix: PHIzf(:,1:2) ------------------------------
%  Maps structural mode accelerations [q_f1_ddot; q_f2_ddot] to IMU
%  measurements u_z_ddot. Purely structural, velocity-independent.
PHIzf = PHIzg * PHIgf;           % 8 x 5
PHIzf_12 = PHIzf(:,1:2);         % 8 x 2
cond_output = cond(PHIzf_12);

% ---- Input blending matrix: PHIfx(1:2,:) --------------------------
%  Maps actuator deflections to generalized aerodynamic forces on modes 1&2.
%  Dynamic pressure q_bar = 0.5*rho*V^2 is a scalar factor that does not
%  change the direction of pinv(PHIfx) nor its condition number.
sumQLpjjBjx = zeros(num_panels, num_AIL);
for i = 1:num_poles
    sumQLpjjBjx = sumQLpjjBjx + QLpjj(:,:,i)*(DRe_jx - DIm_jx*poles(i)*2/c_ref);
end
PHIfx = Sfj*(Q0jj*DRe_jx + sumQLpjjBjx);   % without q_bar
PHIfx_12 = PHIfx(1:2,:);       % 2 x 8
cond_input = cond(PHIfx_12);

% ---- Singular values for additional insight -----------------------------
sv_output = svd(PHIzf_12);
sv_input  = svd(PHIfx_12);

% ---- Display results ----------------------------------------------------
fprintf('=== Condition Number Analysis ===\n');
fprintf('Output blending PHIzf(:,1:2):  cond = %.2f\n', cond_output);
fprintf('  sigma_max = %.4f,  sigma_min = %.4f\n', sv_output(1), sv_output(end));
fprintf('Input blending  PHIfx(1:2,:):  cond = %.2f\n', cond_input);
fprintf('  sigma_max = %.4f,  sigma_min = %.4f\n', sv_input(1), sv_input(end));

% ---- Plot: Bar chart of condition numbers and singular values -----------
fig = figure('Name','Condition Number of Blending Matrices');
set(fig, 'Units','centimeters', ...
         'Position', [7,7,figwidth,figheight], ...
         'defaultAxesFontSize', docFontSize, ...
         'defaultTextFontSize', docFontSize, ...
         'defaultLegendFontSize', docFontSize, ...
         'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex');
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperSize', [figwidth figheight]);
set(fig, 'PaperPosition', [0 0 figwidth figheight]);

t = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% Singular values of PHIzf(:,1:2)
ax1 = nexttile(t);
bar(ax1, sv_output, 'FaceColor', [0.3 0.3 0.3]);
xticklabels(ax1, {'$\sigma_1$','$\sigma_2$'});
ylabel(ax1, 'Singular Value', 'Interpreter', 'latex');
title(ax1, ['$\Phi_{zf}$, $n_f{=}2$, $\kappa = ' sprintf('%.1f', cond_output) '$'], ...
    'Interpreter', 'latex', 'FontSize', docFontSize);
grid(ax1, 'on');
set(ax1, 'Color', 'w');

% Singular values of PHIfx(1:2,:)
ax2 = nexttile(t);
bar(ax2, sv_input, 'FaceColor', [0.3 0.3 0.3]);
xticklabels(ax2, {'$\sigma_1$','$\sigma_2$'});
ylabel(ax2, 'Singular Value', 'Interpreter', 'latex');
title(ax2, ['$\Phi_{fx}$, $n_f{=}2$, $\kappa = ' sprintf('%.1f', cond_input) '$'], ...
    'Interpreter', 'latex', 'FontSize', docFontSize);
grid(ax2, 'on');
set(ax2, 'Color', 'w');

set(fig, 'Color', 'w');
FigureName = 'Fig14_Condition_Number_Blending_nf2';
fullFigurePath = fullfile(figures_dir, FigureName);
print(fig, fullFigurePath, '-dpdf', '-vector');
% print(fig, fullFigurePath, '-dmeta', '-vector');

%% ---- All five modes retained -------------------------------------------
PHIzf_15 = PHIzf(:,1:5);
PHIfx_15   = PHIfx(1:5,:);
cond_output_5 = cond(PHIzf_15);
cond_input_5  = cond(PHIfx_15);
sv_output_5   = svd(PHIzf_15);
sv_input_5    = svd(PHIfx_15);

fprintf('\n=== Condition Number — All 5 Modes ===\n');
fprintf('Output blending PHIzf(:,1:5):  cond = %.2f\n', cond_output_5);
fprintf('Input blending  PHIfx(1:5,:):  cond = %.2f\n', cond_input_5);

% ---- Plot: Bar chart of singular values (5 modes) ----------------------
fig2 = figure('Name','Condition Number of Blending Matrices (5 Modes)');
set(fig2, 'Units','centimeters', ...
         'Position', [7,7,figwidth,figheight], ...
         'defaultAxesFontSize', docFontSize, ...
         'defaultTextFontSize', docFontSize, ...
         'defaultLegendFontSize', docFontSize, ...
         'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex');
set(fig2, 'PaperUnits', 'centimeters');
set(fig2, 'PaperSize', [figwidth figheight]);
set(fig2, 'PaperPosition', [0 0 figwidth figheight]);

t2 = tiledlayout(fig2, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% Singular values of PHIzf(:,1:5)
ax1 = nexttile(t2);
bar(ax1, sv_output_5, 'FaceColor', [0.3 0.3 0.3]);
xticklabels(ax1, {'$\sigma_1$','$\sigma_2$','$\sigma_3$','$\sigma_4$','$\sigma_5$'});
ylabel(ax1, 'Singular Value', 'Interpreter', 'latex');
title(ax1, ['$\Phi_{zf}$, $n_f{=}5$, $\kappa = ' sprintf('%.1f', cond_output_5) '$'], ...
    'Interpreter', 'latex', 'FontSize', docFontSize);
grid(ax1, 'on');
set(ax1, 'Color', 'w');

% Singular values of PHIfx(1:5,:)
ax2 = nexttile(t2);
bar(ax2, sv_input_5, 'FaceColor', [0.3 0.3 0.3]);
xticklabels(ax2, {'$\sigma_1$','$\sigma_2$','$\sigma_3$','$\sigma_4$','$\sigma_5$'});
ylabel(ax2, 'Singular Value', 'Interpreter', 'latex');
title(ax2, ['$\Phi_{fx}$, $n_f{=}5$, $\kappa = ' sprintf('%.1f', cond_input_5) '$'], ...
    'Interpreter', 'latex', 'FontSize', docFontSize);
grid(ax2, 'on');
set(ax2, 'Color', 'w');

set(fig2, 'Color', 'w');
FigureName2 = 'Fig15_Condition_Number_Blending_nf5';
fullFigurePath2 = fullfile(figures_dir, FigureName2);
print(fig2, fullFigurePath2, '-dpdf', '-vector');
% print(fig2, fullFigurePath2, '-dmeta', '-vector');
