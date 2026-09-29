%% Structural Blending Vectors Visualization
%  Figure A: SB vs H2 blending vector entries (bar charts)
%  Figure B: Structural mode shapes with ky_SB overlay (geometric duality)
%  Produces: Figures/Fig12_BarChart_Blending_Vectors_SB_vs_H2.pdf
%            Figures/Fig13_Modeshape_Blending_Overlay.pdf


clearvars
set(groot, 'DefaultFigureRenderer', 'painters');
textwidth = 15.98; % cm (IFASD 2026 template: 455.24pt)
docFontSize = 9; % pt
stdLineWidth = 1.2;

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end

% ---- Model parameters ---------------------------------------------------
num_modes = 5;
num_poles = 6;

[Structure, Aero] = define_RectWing_Structure_Aero(num_modes, num_poles);
Sfj    = Structure.Sfj;
PHIgf  = Structure.PHIgf;
DRe_jx = Structure.DRe_jx;
DIm_jx = Structure.DIm_jx;
PHIzg  = Structure.PHIzg;

c_ref  = Aero.c_ref;
rho    = Aero.rho;
poles  = Aero.poles;
Q0jj   = Aero.Q0jj;
QLpjj  = Aero.QLpjj;

num_panels = size(Q0jj,1);
num_AIL = size(DRe_jx,2);

% ---- Structural Blending Vectors (V-independent) -----------------------
PHIzf = PHIzg * PHIgf;                      % 8 x 5
ky_SB = pinv(PHIzf(:,1:2))';                % 8 x 2

sumQLpjjBjx = zeros(num_panels, num_AIL);
for i = 1:num_poles
    sumQLpjjBjx = sumQLpjjBjx + QLpjj(:,:,i)*(DRe_jx - DIm_jx*poles(i)*2/c_ref);
end
Bgx_struct = Sfj*(Q0jj*DRe_jx + sumQLpjjBjx);
ku_SB = pinv(Bgx_struct(1:2,:));             % 8 x 2

% ---- H2 Optimal Blending Vectors (V_inf = 90 m/s, 8 IMUs, 8 AIL) --------
% The code for the calculation of the H_2 optimal blending vectors is proprietary, hence hardcoded here:
ky_H2 = [-0.099544537390600;-0.320113549928847;-0.532225871772634;-0.699695797290769;0.076443526375992;0.180003212532740;0.213141230616535;0.176367966086539];
ku_H2 = [0.044984777487334;0.311962575366850;0.447659003359420;0.560303866429342;0.115582305723231;0.268015425550478;0.385831765058990;0.390203827098918];

% ---- Labels and positions -----------------------------------------------
labels = {'F1','F2','F3','F4','S1','S2','S3','S4'};
s = 7.5;
num_ele = 16;
y_nodes = linspace(0, s, num_ele+1);
y_span  = [0.9375, 2.8125, 4.6875, 6.5625];
bending_dof = 1:3:3*(num_ele+1);
torsion_dof = 3:3:3*(num_ele+1);

% =========================================================================
%% Figure A: SB vs H2 blending vector entries (bar charts)
% =========================================================================
figwidth_A  = 0.9*textwidth;
figheight_A = 0.6*figwidth_A;

fig1 = figure('Name','Blending Vectors Comparison');
set(fig1, 'Units','centimeters', ...
    'Position', [5,5,figwidth_A,figheight_A], ...
    'defaultAxesFontSize', docFontSize, ...
    'defaultTextFontSize', docFontSize, ...
    'defaultLegendFontSize', docFontSize, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');
set(fig1, 'PaperUnits', 'centimeters');
set(fig1, 'PaperSize', [figwidth_A figheight_A]);
set(fig1, 'PaperPosition', [0 0 figwidth_A figheight_A]);

t1 = tiledlayout(fig1, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- (1,1) SB Output ky ---
ax = nexttile(t1);
b = bar(ax, ky_SB);
b(1).FaceColor = [0.2 0.2 0.2];
b(2).FaceColor = [0.65 0.65 0.65];
set(ax, 'XTick', 1:8, 'XTickLabel', labels);
ylabel(ax, 'Weight', 'Interpreter', 'latex');
title(ax, 'SB: $k_{y,\mathrm{SB}}$', 'Interpreter', 'latex', 'FontSize', docFontSize);
legend(ax, {'Mode 1','Mode 2'}, 'Interpreter', 'latex', ...
    'Location', 'best', 'FontSize', docFontSize-1);
grid(ax, 'on');
set(ax, 'Color', 'w');

% --- (1,2) H2 Output |ky| ---
ax = nexttile(t1);
bar(ax, abs(ky_H2), 'FaceColor', [0.2 0.2 0.2]);
set(ax, 'XTick', 1:8, 'XTickLabel', labels);
ylabel(ax, '$|$Weight$|$', 'Interpreter', 'latex');
title(ax, '$H_2$: $|k_{y,H_2}|$', 'Interpreter', 'latex', 'FontSize', docFontSize);
grid(ax, 'on');
set(ax, 'Color', 'w');

% --- (2,1) SB Input ku ---
ax = nexttile(t1);
b = bar(ax, ku_SB);
b(1).FaceColor = [0.2 0.2 0.2];
b(2).FaceColor = [0.65 0.65 0.65];
set(ax, 'XTick', 1:8, 'XTickLabel', labels);
xlabel(ax, 'Actuator', 'Interpreter', 'latex');
ylabel(ax, 'Weight', 'Interpreter', 'latex');
title(ax, 'SB: $k_{u,\mathrm{SB}}$', 'Interpreter', 'latex', 'FontSize', docFontSize);
legend(ax, {'Mode 1','Mode 2'}, 'Interpreter', 'latex', ...
    'Location', 'best', 'FontSize', docFontSize-1);
grid(ax, 'on');
set(ax, 'Color', 'w');

% --- (2,2) H2 Input |ku| ---
ax = nexttile(t1);
bar(ax, abs(ku_H2), 'FaceColor', [0.2 0.2 0.2]);
set(ax, 'XTick', 1:8, 'XTickLabel', labels);
xlabel(ax, 'Actuator', 'Interpreter', 'latex');
ylabel(ax, '$|$Weight$|$', 'Interpreter', 'latex');
title(ax, '$H_2$: $|k_{u,H_2}|$', 'Interpreter', 'latex', 'FontSize', docFontSize);
grid(ax, 'on');
set(ax, 'Color', 'w');

set(fig1, 'Color', 'w');
FigureName = 'Fig12_BarChart_Blending_Vectors_SB_vs_H2';
fullFigurePath = fullfile(figures_dir, FigureName);
print(fig1, fullFigurePath, '-dpdf', '-vector');
% print(fig1, fullFigurePath, '-dmeta', '-vector');

% =========================================================================
%% Figure B: Mode shape overlay with output blending vectors
% =========================================================================
figwidth_B  = 0.9*textwidth;
figheight_B = 0.4*figwidth_B;

fig2 = figure('Name','Mode Shape Overlay');
set(fig2, 'Units','centimeters', ...
    'Position', [5,5,figwidth_B,figheight_B], ...
    'defaultAxesFontSize', docFontSize, ...
    'defaultTextFontSize', docFontSize, ...
    'defaultLegendFontSize', docFontSize, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');
set(fig2, 'PaperUnits', 'centimeters');
set(fig2, 'PaperSize', [figwidth_B figheight_B]);
set(fig2, 'PaperPosition', [0 0 figwidth_B figheight_B]);

t2 = tiledlayout(fig2, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

mode1_bend = PHIgf(bending_dof, 1);
mode2_tors = PHIgf(torsion_dof, 2);

% --- Mode 1: Bending shape + ky_SB(:,1) ---
ax = nexttile(t2);
plot(ax, y_nodes, mode1_bend/max(abs(mode1_bend)), 'k-', 'LineWidth', stdLineWidth);
hold(ax, 'on');
% Flap IMUs (filled circles)
s1 = stem(ax, y_span, ky_SB(1:4,1)/max(abs(ky_SB(:,1))), 'filled');
s1.Color = [0.3 0.3 0.3]; s1.MarkerSize = 5; s1.LineWidth = 0.8;
% Slat IMUs (filled triangles)
s2 = stem(ax, y_span, ky_SB(5:8,1)/max(abs(ky_SB(:,1))), 'filled');
s2.Color = [0.6 0.6 0.6]; s2.Marker = '^'; s2.MarkerSize = 5;
s2.LineWidth = 0.8; s2.MarkerFaceColor = [0.6 0.6 0.6];
hold(ax, 'off');
xlabel(ax, '$y$ [m]', 'Interpreter', 'latex');
ylabel(ax, 'Norm.\ amplitude', 'Interpreter', 'latex');
title(ax, 'Bending (Mode 1)', ...
    'Interpreter', 'latex', 'FontSize', docFontSize);
legend(ax, {'$\phi_{f_1}^{u_z}(y)$', '$k_y$ flap', '$k_y$ slat'}, ...
    'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', docFontSize-1);
grid(ax, 'on');
set(ax, 'Color', 'w');

% --- Mode 2: Torsion shape + ky_SB(:,2) ---
ax = nexttile(t2);
plot(ax, y_nodes, mode2_tors/max(abs(mode2_tors)), 'k-', 'LineWidth', stdLineWidth);
hold(ax, 'on');
s1 = stem(ax, y_span, ky_SB(1:4,2)/max(abs(ky_SB(:,2))), 'filled');
s1.Color = [0.3 0.3 0.3]; s1.MarkerSize = 5; s1.LineWidth = 0.8;
s2 = stem(ax, y_span, ky_SB(5:8,2)/max(abs(ky_SB(:,2))), 'filled');
s2.Color = [0.6 0.6 0.6]; s2.Marker = '^'; s2.MarkerSize = 5;
s2.LineWidth = 0.8; s2.MarkerFaceColor = [0.6 0.6 0.6];
hold(ax, 'off');
xlabel(ax, '$y$ [m]', 'Interpreter', 'latex');
ylabel(ax, 'Norm.\ amplitude', 'Interpreter', 'latex');
title(ax, 'Torsion (Mode 2)', ...
    'Interpreter', 'latex', 'FontSize', docFontSize);
legend(ax, {'$\phi_{f_2}^{\psi_y}(y)$', '$k_y$ flap', '$k_y$ slat'}, ...
    'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', docFontSize-1);
grid(ax, 'on');
set(ax, 'Color', 'w');

set(fig2, 'Color', 'w');
FigureName = 'Fig13_Modeshape_Blending_Overlay';
fullFigurePath = fullfile(figures_dir, FigureName);
print(fig2, fullFigurePath, '-dpdf', '-vector');
% print(fig2, fullFigurePath, '-dmeta', '-vector');

