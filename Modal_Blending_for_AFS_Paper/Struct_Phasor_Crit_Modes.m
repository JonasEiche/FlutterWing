% Eigenvector analysis of the critical aeroelastic modes
clearvars
textwidth = 16.5730; % cm
figwidth = 0.50*textwidth; % cm
figheight = 1.0*figwidth;
docFontSize = 9; % pt
stdLineWidth = 1.2; 

V_inf_arr = [20,40,60,80,100,120,140,160];
nv = length(V_inf_arr);
imuIDX = 1:8;
ailIDX = 4;
ny=length(imuIDX);
nu=length(ailIDX);
G = build_G_RectWing(100,imuIDX,ailIDX);
nx=length(G.A);

G_arr = ss(zeros(ny,nu,1,nv));
PHIV = zeros(nx,nx,nv);
PHIW = zeros(nx,nx,nv);
EV   = zeros(nx,nv);

for iv = 1:nv
    V_inf = V_inf_arr(iv);
    G = build_G_RectWing(V_inf,imuIDX,ailIDX);
    G_arr(:,:,1,iv) = G;
    A = G.A; B=G.B; C=G.C; D=G.D;
    [V,DD] = eig(A);           % DD = V\G.A*V
    % W=inv(V);
    PHIV(:,:,iv)  =  V;
    EV(:,iv)      =  diag(DD);
    % track modes
    if iv>=2
	    idx = eigenshuffle(EV(:,iv-1),PHIV(:,:,iv-1),EV(:,iv),PHIV(:,:,iv));
    else
        idx = 1:size(EV,1);
    end
    PHIV(:,:,iv) = PHIV(:,idx,iv);
    PHIW(:,:,iv) = inv(PHIV(:,:,iv));
    EV(:,iv)     = EV(idx,iv);
end

[~, idx] = sort(EV(:,end),'descend','ComparisonMethod','real');
EV = EV(idx,:);
PHIV = PHIV(:,idx,:);   % PHIV(state,eigenmode,velocity)
PHIW = PHIW(idx,:,:);   % PHIW(eigenmode,state,velocity

EVval = EV(:,6);    % Evaluate criticallity at V_inf_arr(6) = 120 m/s
flutIDX = find( (real(EVval) > -1)  & (abs(imag(EVval)) < 40) & (abs(imag(EVval)) > 5) );
critIDX = find( (real(EVval) > -10)  & (abs(imag(EVval)) < 40) & (abs(imag(EVval)) > 5) );
resIDX = find( (real(EVval) > -10)  & (abs(imag(EVval)) < 120) & (abs(imag(EVval)) > 40) );

PHIV_flut = PHIV(:,flutIDX,:);
PHIW_flut = PHIW(flutIDX,:,:);
PHIV_crit = PHIV(:,critIDX,:);
PHIW_crit = PHIW(critIDX,:,:);
PHIV_res = PHIV(:,resIDX,:);
PHIW_res = PHIW(resIDX,:,:);

num_modes=5;

PHIV_crit_x12 = PHIV_crit(1:2,:,:);

PHIV_crit_x67 = PHIV_crit(num_modes+1:2,:,:);


% PHIV_crit_x12(:,:,6)    % PHIV(state,eigenmode,velocity)

PHIV_crit_x12_mag_v6 = abs(PHIV_crit_x12(:,[1,3],6));

PHIV_crit_x12_mag = abs(PHIV_crit_x12(:,[1,3],:));

PHIV_crit_x12_mag_lengths = vecnorm(PHIV_crit_x12_mag, 2, 1); % 2norm along dimension 1

PHIV_crit_x12_mag_norma = PHIV_crit_x12_mag./PHIV_crit_x12_mag_lengths;

PHIV_crit_x1_mag_m1 = squeeze(PHIV_crit_x12_mag_norma(1,1,:));
PHIV_crit_x2_mag_m1 = squeeze(PHIV_crit_x12_mag_norma(2,1,:));

PHIV_crit_x1_mag_m2 = squeeze(PHIV_crit_x12_mag_norma(1,2,:));
PHIV_crit_x2_mag_m2 = squeeze(PHIV_crit_x12_mag_norma(2,2,:));

% PHIV_crit_x1_mag_m1 = squeeze(PHIV_crit_x12_mag(1,1,:));
% PHIV_crit_x2_mag_m1 = squeeze(PHIV_crit_x12_mag(2,1,:));
% 
% PHIV_crit_x1_mag_m2 = squeeze(PHIV_crit_x12_mag(1,2,:));
% PHIV_crit_x2_mag_m2 = squeeze(PHIV_crit_x12_mag(2,2,:));

fig = figure('Name','Critical Modes Bending Torsion Magnitude over V_infty');
% t = tiledlayout(fig, 1, 1, 'Padding', 'tight', 'TileSpacing', 'none');
set(fig, 'Units','centimeters', ...
         'Position', [7,7,figwidth,figheight], ...
         'defaultAxesFontSize', docFontSize, ...
         'defaultTextFontSize', docFontSize, ...
         'defaultLegendFontSize', docFontSize, ...
         'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex');
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperSize', [figwidth figheight]);        % Set PDF page size to match figure
set(fig, 'PaperPosition', [0 0 figwidth figheight]); % Position plot to fill the PDF page exactly
% ax = nexttile;

quiver(zeros(nv,1), zeros(nv,1), PHIV_crit_x1_mag_m1, PHIV_crit_x2_mag_m1, 0, 'k-', 'LineWidth', stdLineWidth, 'MaxHeadSize', 0.1) % First group in blue
hold on
quiver(zeros(nv,1), zeros(nv,1), PHIV_crit_x1_mag_m2, PHIV_crit_x2_mag_m2, 0, 'k--', 'LineWidth', stdLineWidth, 'MaxHeadSize', 0.1) % Second group in red
% quiver(zeros(nv,1), zeros(nv,1), PHIV_crit_x12_mag_rshp(1,:), PHIV_crit_x12_mag_rshp(2,:), 0, 'LineWidth', 2, 'MaxHeadSize', 0.5);

% V_inf_arr = [20,40,60,80,100,120,140,160];

text(PHIV_crit_x1_mag_m1(1)*1.01, PHIV_crit_x2_mag_m1(1)*1.01, '20m/s', 'Interpreter', 'latex');
% text(PHIV_crit_x1_mag_m1(2)*1.05, PHIV_crit_x2_mag_m1(2)*1.05, '40m/s', 'Interpreter', 'latex');
text(PHIV_crit_x1_mag_m1(3)*1.01, PHIV_crit_x2_mag_m1(3)*1.01, '60m/s', 'Interpreter', 'latex');
text(PHIV_crit_x1_mag_m1(4)*1.01, PHIV_crit_x2_mag_m1(4)*1.01, '80m/s', 'Interpreter', 'latex');
text(PHIV_crit_x1_mag_m1(5)*1.01, PHIV_crit_x2_mag_m1(5)*1.01, '100m/s', 'Interpreter', 'latex');
text(PHIV_crit_x1_mag_m1(6)*1.01, PHIV_crit_x2_mag_m1(6)*1.01, '120m/s', 'Interpreter', 'latex');
text(PHIV_crit_x1_mag_m1(7)*1.01, PHIV_crit_x2_mag_m1(7)*1.01, '140m/s', 'Interpreter', 'latex');
text(PHIV_crit_x1_mag_m1(8)*1.01, PHIV_crit_x2_mag_m1(8)*1.01, '160m/s', 'Interpreter', 'latex');

text(PHIV_crit_x1_mag_m2(1)*1.01, PHIV_crit_x2_mag_m2(1)*1.01, '20m/s', 'Interpreter', 'latex');
% text(PHIV_crit_x1_mag_m2(2)*1.01, PHIV_crit_x2_mag_m2(2)*1.01, '40m/s', 'Interpreter', 'latex');
text(PHIV_crit_x1_mag_m2(3)*1.01, PHIV_crit_x2_mag_m2(3)*1.01, '60m/s', 'Interpreter', 'latex');
text(PHIV_crit_x1_mag_m2(4)*1.01, PHIV_crit_x2_mag_m2(4)*1.01, '80m/s', 'Interpreter', 'latex');
% text(PHIV_crit_x1_mag_m2(5)*1.10, PHIV_crit_x2_mag_m2(5)*1.10, '100m/s', 'Interpreter', 'latex');
% text(PHIV_crit_x1_mag_m2(6), PHIV_crit_x2_mag_m2(6), '120m/s');
text(PHIV_crit_x1_mag_m2(7)*1.01, PHIV_crit_x2_mag_m2(7)*1.01, '140m/s', 'Interpreter', 'latex');
text(PHIV_crit_x1_mag_m2(8)*1.01, PHIV_crit_x2_mag_m2(8)*1.01, '160m/s', 'Interpreter', 'latex');

axis equal;
% grid on;
xlabel('Bending', 'Interpreter', 'latex', 'FontSize', docFontSize);
ylabel('Torsion', 'Interpreter', 'latex', 'FontSize', docFontSize);
% title('Critical Eigenmodes over $V_\infty$', 'Interpreter', 'latex');
axis([0, 1.1, 0, 1.1])   % xmin xmax ymin ymax
set(fig, 'Color', 'w'); % Set figure background to white
ax = fig.CurrentAxes;
ax.TickLabelInterpreter = 'latex'; % Set tick labels to LaTeX
ax.FontSize = docFontSize;
set(ax, 'Color', 'w'); % Set axes background to white
% exportgraphics(f_phasor,'C:\Users\eich_j0\Documents\MATLAB\Flatter\Plots_for_Modal_Control_Paper\Figures\bending_torsion_phasor_critmodes.pdf','ContentType','vector', 'BackgroundColor', 'none');


figureName = 'bending_torsion_phasor_critmodes';
script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
fullFigurePath = fullfile(figures_dir, figureName);
print(fig, fullFigurePath, '-dpdf', '-vector');
print(fig, fullFigurePath, '-dmeta', '-vector');