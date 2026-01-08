% Vizualisation of the Aerodynamic Panel Mesh with Control Surfaces and IMU
% Positions
clearvars
textwidth = 16.5730; % cm
figwidth = 0.5*textwidth; % cm
figheight = 0.7*figwidth;
docFontSize = 9; % pt
stdLineWidth = 1.2; 

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end

% ---- Geometry -----------------------------------------------------------
s = 7.5;                % semi span
c = 2;                  % root chord
c_ref = c;              % reference chord

np_s = 20;              % number panels spanwise
np_c = 10;              % number panels chordwise

sw = 0;                 % sweep angle at leading edge
dh = 0;                 % dihedral
tr = 1;                 % taper ratio
xf = 0.48;              % flex axis position relative to chord
[Pa,Ps] = build_PaPs(s,c,np_s,np_c,sw,dh,tr,xf);

%  21 22 23 24 25 26 27 28 29 30 31 32 ...
%  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20
cspanels_flap{1} = [1,2,3,4,5];
cspanels_flap{2} = [6,7,8,9,10];
cspanels_flap{3} = [11,12,13,14,15];            
cspanels_flap{4} = [16,17,18,19,20];           % panel numbers belonging to control surface cs

cspanels_slat{1} = 180+[1,2,3,4,5];
cspanels_slat{2} = 180+[6,7,8,9,10];
cspanels_slat{3} = 180+[11,12,13,14,15]; 
cspanels_slat{4} = 180+[16,17,18,19,20];  
num_FLAP = length(cspanels_flap);
num_SLAT = length(cspanels_slat);
num_AIL = length(cspanels_flap)+length(cspanels_slat);         % Number of ailerons



pos_IMU = { 0.5*(Pa{3}{1}+Pa{3}{4});
            0.5*(Pa{8}{1}+Pa{8}{4});
            0.5*(Pa{13}{1}+Pa{13}{4});
            0.5*(Pa{18}{1}+Pa{18}{4}) };


x_IMU = [ pos_IMU{1}(1); pos_IMU{1}(1)-8/10*c_ref;
          pos_IMU{2}(1); pos_IMU{2}(1)-8/10*c_ref;
          pos_IMU{3}(1); pos_IMU{3}(1)-8/10*c_ref;
          pos_IMU{4}(1); pos_IMU{4}(1)-8/10*c_ref ];

y_IMU = [ pos_IMU{1}(2); pos_IMU{1}(2);
          pos_IMU{2}(2); pos_IMU{2}(2);
          pos_IMU{3}(2); pos_IMU{3}(2);
          pos_IMU{4}(2); pos_IMU{4}(2) ];

z_IMU = [ pos_IMU{1}(3); pos_IMU{1}(3);
          pos_IMU{2}(3); pos_IMU{2}(3);
          pos_IMU{3}(3); pos_IMU{3}(3);
          pos_IMU{4}(3); pos_IMU{4}(3) ];


u_phi=[5;10;15;20];

Q=Pa;
for ics = 1:num_FLAP
    P_ics=Pa{cspanels_flap{ics}};
    P_ics_X = cellfun(@(x) x(1), P_ics);
    xROT = min(P_ics_X);
    for ip = cspanels_flap{ics}
        Z_CS = [(Pa{ip}{2}(1) - xROT), (Pa{ip}{1}(1) - xROT);
                (Pa{ip}{3}(1) - xROT), (Pa{ip}{4}(1) - xROT)].*u_phi(ics)/180*pi;
        Q{ip}{2}(3)=Pa{ip}{2}(3)+Z_CS(1,1);
        Q{ip}{1}(3)=Pa{ip}{1}(3)+Z_CS(1,2);
        Q{ip}{3}(3)=Pa{ip}{3}(3)+Z_CS(2,1);
        Q{ip}{4}(3)=Pa{ip}{4}(3)+Z_CS(2,2);
    end
end

for ics = 1:num_SLAT
    P_ics=Pa{cspanels_slat{ics}};
    P_ics_X = cellfun(@(x) x(1), P_ics);
    xROT = max(P_ics_X);
    for ip = cspanels_slat{ics}
        Z_CS = [(Pa{ip}{2}(1) - xROT), (Pa{ip}{1}(1) - xROT);
                (Pa{ip}{3}(1) - xROT), (Pa{ip}{4}(1) - xROT)].*u_phi(ics)/180*pi*(-1); 
        Q{ip}{2}(3)=Pa{ip}{2}(3)+Z_CS(1,1);
        Q{ip}{1}(3)=Pa{ip}{1}(3)+Z_CS(1,2);
        Q{ip}{3}(3)=Pa{ip}{3}(3)+Z_CS(2,1);
        Q{ip}{4}(3)=Pa{ip}{4}(3)+Z_CS(2,2);
    end
end

fig = figure('Name','Flatterfluegel');
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

for i = 1:length(Q)
    X = [Q{i}{2}(1), Q{i}{1}(1);
         Q{i}{3}(1), Q{i}{4}(1)];
    Y = [Q{i}{2}(2), Q{i}{1}(2);
         Q{i}{3}(2), Q{i}{4}(2)];
    Z = [Q{i}{2}(3), Q{i}{1}(3);
         Q{i}{3}(3), Q{i}{4}(3)];

    surf(X,Y,Z)
    hold on
end
% Plot red dots for IMUs
plot3(x_IMU, y_IMU, z_IMU,'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r');


axis equal
colormap("white"); % Apply colormap:  gray, bone, pink, parula

xlabel('$x_B$', 'Interpreter', 'latex', 'FontSize', docFontSize)
ylabel('$y_B$', 'Interpreter', 'latex', 'FontSize', docFontSize)
zlabel('$z_B$', 'Interpreter', 'latex', 'FontSize', docFontSize)

% axis off

% Remove ticks and tick labels
xticks([])      % Remove x-axis ticks
yticks([])      % Remove y-axis ticks
zticks([])      % Remove z-axis ticks
% 
xticklabels([]) % Remove x-axis tick labels
yticklabels([]) % Remove y-axis tick labels
zticklabels([]) % Remove z-axis tick labels

% colorbar
fig.CurrentAxes.ZDir = 'Reverse';
fig.CurrentAxes.XDir = 'Reverse';

fig.CurrentAxes.CameraPosition = [ -22.8970, -19.5183, -19.7860];

set(fig, 'Color', 'w'); % Set figure background to white

%ax = fig.CurrentAxes; % Get current axes
ax.TickLabelInterpreter = 'latex'; % Set tick labels to LaTeX
ax.FontSize = docFontSize;
set(ax, 'Color', 'w'); % Set axes background to white
axis off
set(fig, 'Color', 'w'); % Set figure background to white

FigureName = 'aero_panels_slat_flap_IMU';
fullVpzmap_PlotPath = fullfile(figures_dir, FigureName);
print(fig, fullVpzmap_PlotPath, '-dpdf', '-vector');
print(fig, fullVpzmap_PlotPath, '-dmeta', '-vector');





