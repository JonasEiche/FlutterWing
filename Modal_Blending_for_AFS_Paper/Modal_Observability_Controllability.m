% Eigenvector analysis of the critical aeroelastic modes Input Output
% Controllability / Observability
clearvars
textwidth = 16.5730; % cm
figwidth = 0.47*textwidth; % cm
figheight = 0.7*figwidth;
docFontSize = 9; % pt
stdLineWidth = 1.2; 

% V_inf_arr = [20,40,60,80,100,120,140,160];
V_inf_arr = linspace(20,160,32);
nv = length(V_inf_arr);
imuIDX = 1:8;
ailIDX = 1:8;
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

PHIW_normalized = PHIW ./ vecnorm(PHIW, 2, 2);

EVval = EV(:,6);    % Evaluate criticallity at V_inf_arr(6) = 120 m/s
flutIDX = find( (real(EVval) > -1)  & (abs(imag(EVval)) < 40) & (abs(imag(EVval)) > 5) );
critIDX = find( (real(EVval) > -10)  & (abs(imag(EVval)) < 40) & (abs(imag(EVval)) > 5) );
resIDX = find( (real(EVval) > -10)  & (abs(imag(EVval)) < 120) & (abs(imag(EVval)) > 40) );

PHIV_flut = PHIV(:,flutIDX,:);
PHIW_flut = PHIW(flutIDX,:,:);
PHIV_crit = PHIV(:,critIDX,:);
PHIW_crit = PHIW(critIDX,:,:);

PHIW_norm_flut = PHIW_normalized(flutIDX,:,:);

NORM_CPHIV_flut_Re = zeros(nv,1);
NORM_CPHIV_flut_Im = zeros(nv,1);
NORM_PHIWB_flut_Re = zeros(nv,1);
NORM_PHIWB_flut_Im = zeros(nv,1);
NORM_PHIWB_norm_flut_Re = zeros(nv,1);
NORM_PHIWB_norm_flut_Im = zeros(nv,1);
for iv = 1:nv
    C= G_arr(:,:,1,iv).C;
    Cphiv_crit = C*PHIV_crit(:,:,iv);
    Cphiv_flut = C*PHIV_flut(:,:,iv);

    NORM_CPHIV_flut_Re(iv) = norm(real(Cphiv_flut(:,1)));
    NORM_CPHIV_flut_Im(iv) = norm(imag(Cphiv_flut(:,1)));

    B= G_arr(:,:,1,iv).B;
    phiwB_crit = PHIW_crit(:,:,iv)*B;
    phiwB_flut = PHIW_flut(:,:,iv)*B;

    phiwB_norm_flut = PHIW_norm_flut(:,:,iv)*B;

    NORM_PHIWB_flut_Re(iv) = norm(real(phiwB_flut(1,:)));
    NORM_PHIWB_flut_Im(iv) = norm(imag(phiwB_flut(1,:)));

    NORM_PHIWB_norm_flut_Re(iv) = norm(real(phiwB_norm_flut(1,:)));
    NORM_PHIWB_norm_flut_Im(iv) = norm(imag(phiwB_norm_flut(1,:)));
end


fig = figure('Name','Flutter Modal Observability');
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
plot(V_inf_arr,NORM_CPHIV_flut_Re, 'k-', ...
     V_inf_arr,NORM_CPHIV_flut_Im, 'k--','LineWidth',stdLineWidth)

% h_title = title('Volume of Fluttermode Output Projection', 'Interpreter', 'latex');
h_xlabel = xlabel('Velocity $V_\infty$', 'Interpreter', 'latex');
h_ylabel = ylabel('Modal Observability', 'Interpreter', 'latex');
h_legend = legend({'Fluttermode Real','Fluttermode Imag'},'Location','best', 'Interpreter', 'latex');

% set(h_title, 'FontSize', docFontSize);
set(h_xlabel, 'FontSize', docFontSize);
set(h_ylabel, 'FontSize', docFontSize);
set(h_legend, 'FontSize', docFontSize);

grid on
grid minor;
set(fig, 'Color', 'w'); % Set figure background to white
ax.TickLabelInterpreter = 'latex'; % Set tick labels to LaTeX
set(ax, 'Color', 'w'); % Set axes background to white

ax.YMinorTick = 'on';

figureName = 'fluttermode_observability';
script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
fullFigurePath = fullfile(figures_dir, figureName);
print(fig, fullFigurePath, '-dpdf', '-vector');
print(fig, fullFigurePath, '-dmeta', '-vector');


% Controllability

fig2 = figure('Name','Flutter Modal Controllability');
t = tiledlayout(fig2, 1, 1, 'Padding', 'none', 'TileSpacing', 'none');
set(fig2, 'Units','centimeters', ...
         'Position', [7,7,figwidth,figheight], ...
         'defaultAxesFontSize', docFontSize, ...
         'defaultTextFontSize', docFontSize, ...
         'defaultLegendFontSize', docFontSize, ...
         'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex');
set(fig2, 'PaperUnits', 'centimeters');
set(fig2, 'PaperSize', [figwidth figheight]);        % Set PDF page size to match figure
set(fig2, 'PaperPosition', [0 0 figwidth figheight]); % Position plot to fill the PDF page exactly
ax = nexttile;
% plot(V_inf_arr,NORM_PHIWB_flut_Re, 'k-', ...
%      V_inf_arr,NORM_PHIWB_flut_Im, 'k--','LineWidth',stdLineWidth)
plot(V_inf_arr,NORM_PHIWB_norm_flut_Re, 'k-', ...
     V_inf_arr,NORM_PHIWB_norm_flut_Im, 'k--','LineWidth',stdLineWidth)


% h_title = title('Volume of Fluttermode Output Projection', 'Interpreter', 'latex');
h_xlabel = xlabel('Velocity $V_\infty$', 'Interpreter', 'latex');
h_ylabel = ylabel('Modal Controllability', 'Interpreter', 'latex');
h_legend = legend({'Fluttermode Real','Fluttermode Imag'},'Location','best', 'Interpreter', 'latex');

% set(h_title, 'FontSize', docFontSize);
set(h_xlabel, 'FontSize', docFontSize);
set(h_ylabel, 'FontSize', docFontSize);
set(h_legend, 'FontSize', docFontSize);

grid on
grid minor;
set(fig2, 'Color', 'w'); % Set figure background to white
ax.TickLabelInterpreter = 'latex'; % Set tick labels to LaTeX
set(ax, 'Color', 'w'); % Set axes background to white

ax.YMinorTick = 'on';

figureName = 'fluttermode_controllability';
script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
fullFigurePath = fullfile(figures_dir, figureName);
print(fig2, fullFigurePath, '-dpdf', '-vector');
print(fig2, fullFigurePath, '-dmeta', '-vector');