% Eigenvector analysis of the critical aeroelastic modes
clearvars
textwidth = 16.5730; % cm
figwidth = 0.50*textwidth; % cm
figheight = 0.7*figwidth;
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

GRAMDET_PHIV_crit = zeros(nv,1);
GRAMDET_PHIV_flut = zeros(nv,1);
GRAMDET_CPHIV_crit = zeros(nv,1);
GRAMDET_CPHIV_flut = zeros(nv,1);
for iv = 1:nv
    GRAMDET_PHIV_crit(iv) = gramDet([PHIV_crit(:,1,iv), PHIV_crit(:,2,iv), PHIV_crit(:,3,iv), PHIV_crit(:,4,iv)],'normalize',true);
    GRAMDET_PHIV_flut(iv) = gramDet([PHIV_flut(:,1,iv), PHIV_flut(:,2,iv)],'normalize',true);
    C= G_arr(:,:,1,iv).C;
    Cphiv_crit = C*PHIV_crit(:,:,iv);
    Cphiv_flut = C*PHIV_flut(:,:,iv);
    GRAMDET_CPHIV_crit(iv) = gramDet([Cphiv_crit(:,1), Cphiv_crit(:,2), Cphiv_crit(:,3), Cphiv_crit(:,4)],'normalize',true);
    GRAMDET_CPHIV_flut(iv) = gramDet([Cphiv_flut(:,1), Cphiv_flut(:,2)],'normalize',true);
end

fig = figure('Name','Gram determinante of Fluttermode Output Projection');
t = tiledlayout(fig, 1, 1, 'Padding', 'none', 'TileSpacing', 'none');
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
ax = nexttile;
plot(V_inf_arr,GRAMDET_CPHIV_flut, 'k-', ...
    V_inf_arr,GRAMDET_CPHIV_crit, 'k--','LineWidth',stdLineWidth)

% h_title = title('Volume of Fluttermode Output Projection', 'Interpreter', 'latex');
h_xlabel = xlabel('Velocity $V_\infty$', 'Interpreter', 'latex');
h_ylabel = ylabel('Gram Determinante', 'Interpreter', 'latex');
h_legend = legend({'Fluttermode','Critical Modes'},'Location','southeast', 'Interpreter', 'latex');

% set(h_title, 'FontSize', docFontSize);
set(h_xlabel, 'FontSize', docFontSize);
set(h_ylabel, 'FontSize', docFontSize);
set(h_legend, 'FontSize', docFontSize);

grid on
grid minor;
set(fig, 'Color', 'w'); % Set figure background to white
ax.TickLabelInterpreter = 'latex'; % Set tick labels to LaTeX
set(ax, 'Color', 'w'); % Set axes background to white

figureName = 'gram_det_fluttermode_output';
script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
fullFigurePath = fullfile(figures_dir, figureName);
% print(fig, fullFigurePath, '-dpdf', '-vector');
% print(fig, fullFigurePath, '-dmeta', '-vector');
