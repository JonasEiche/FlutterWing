% Flutter mode bending/torsion composition over time
clearvars
textwidth = 15.98; % cm (IFASD 2026 template: 455.24pt)
figwidth = 0.7*textwidth; % cm
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
ailIDX=1:8;
imuIDX=1:8;
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
set(fig, 'PaperSize', [figwidth figheight]);
set(fig, 'PaperPosition', [0 0 figwidth figheight]);
normFactor = max(abs(Q(1,:)));
plot(t,Q(1,:)/normFactor, 'k-', ...
    t,Q(2,:)/normFactor, 'k--','LineWidth',stdLineWidth)
hold on
plot(t,Q(3,:)/normFactor, '-', 'Color',[0.5 0.5 0.5], 'LineWidth',stdLineWidth)
plot(t,Q(4,:)/normFactor, '--', 'Color',[0.5 0.5 0.5], 'LineWidth',stdLineWidth)
plot(t,Q(5,:)/normFactor, ':', 'Color',[0.5 0.5 0.5], 'LineWidth',stdLineWidth)
hold off
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', docFontSize)
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', docFontSize)
legend({'1st Bending','1st Torsion','2nd Torsion','2nd Bending','3rd Bending'}, ...
    'Location','southeast', 'Interpreter', 'latex', 'FontSize', docFontSize)
grid on

axis([0, 0.35, -1, 1])
ax = fig.CurrentAxes;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = docFontSize;
set(ax, 'Color', 'w');
set(fig, 'Color', 'w');

FigureName = 'Fig11_Flutter_Structural_Contributions';
FigurePath = fullfile(figures_dir, FigureName);
print(fig, FigurePath, '-dpdf', '-vector');
% print(fig, FigurePath, '-dmeta', '-vector');
