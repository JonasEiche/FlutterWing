clearvars
textwidth = 16.5730; % cm
figwidth = 0.70*textwidth; % cm
figheight = 0.5*figwidth;
docFontSize = 9; % pt
stdLineWidth = 1.2;

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end

V_inf=130;
modesIDX=[1,2];             % modes for performance output
ailIDX=4;                 % AIL used
imuIDX=[4,8];               % IMU (acc) used
num_modes=5;
G = build_G_RectWing(V_inf,imuIDX,ailIDX);
[V,DD] = eig(G.A);
flutIDX = find( (real(diag(DD)) > -1)  & (abs(imag(diag(DD))) < 40) & (abs(imag(diag(DD))) > 5) );
wi_flut=abs(imag(DD(flutIDX(1),flutIDX(1))));
V_flut=V(:,flutIDX);

q_f_flut=V_flut(1:num_modes,1);
q = @(t) real(q_f_flut).*cos(wi_flut*t)-imag(q_f_flut).*sin(wi_flut*t);
t = linspace(0,3*pi/wi_flut,36);
Q = q(t);


fig = figure('Name','Essential Eigenmode Contribution to Fluttermode');
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
plot(t,Q(1,:)/max(abs(Q(1,:))), 'k-', ...
    t,Q(2,:)/max(abs(Q(1,:))), 'k--','LineWidth',stdLineWidth)
% title('Essential Eigenmode Contribution to Fluttermode', 'Interpreter', 'latex', 'FontSize', docFontSize)
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', docFontSize) 
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', docFontSize) 
legend({'Bending Mode','Torsion Mode'},'Location','southeast', 'Interpreter', 'latex', 'FontSize', docFontSize)
grid on

axis([0, 0.35, -1, 1])   % xmin xmax ymin ymax
ax = fig.CurrentAxes; % Get current axes
ax.TickLabelInterpreter = 'latex'; % Set tick labels to LaTeX
ax.FontSize = docFontSize;
set(ax, 'Color', 'w'); % Set axes background to white
set(fig, 'Color', 'w'); % Set figure background to white

FigureName = 'essential_eigenmodes_contribution_fluttermode';
FigurePath = fullfile(figures_dir, FigureName);
print(fig, FigurePath, '-dpdf', '-vector');
print(fig, FigurePath, '-dmeta', '-vector');