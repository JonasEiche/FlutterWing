clearvars
textwidth = 16.5730; % cm
figwidth = 0.7*textwidth; % cm
figheight = 1.0*figwidth;
docFontSize = 9; % pt
stdLineWidth = 1.2; 

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);
figures_dir = fullfile(script_dir, 'Figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
data_dir = fullfile(script_dir, 'Data');

% % --- Produce Qjj DLM AIC Data --------------------------------------------
% 
% % ---- Geometry -----------------------------------------------------------
% s = 7.5;                % semi span
% c = 2;                  % root chord
% c_ref = c;              % reference chord
% 
% np_s = 20;              % number panels spanwise
% np_c = 10;              % number panels chordwise
% 
% sw = 0;                 % sweep angle at leading edge
% dh = 0;                 % dihedral
% tr = 1;                 % taper ratio
% xf = 0.48;              % flex axis position relative to chord
% [Pa,Ps] = build_PaPs(s,c,np_s,np_c,sw,dh,tr,xf);
% % ---- Unsteady Aerodynamics ----------------------------------------------
% Ma = 0.0;
% k_red = [0, 0.001, 0.01, 0.02, 0.05, 0.07, 0.1, 0.2, 0.5, 0.7, 0.9, 1.1];
% SYM =1;
% Qjj = build_Qjj(Ma,k_red,c_ref,Pa,SYM);
% [ny,nu,nk] = size(Qjj);
% save(fullfile(data_dir,'Qjj_npc10nps20_k_red12.mat'),'Qjj','k_red');
% 
% % --- Rational Function Approximation -------------------------------------
% num_poles = 6;
% [poles_rog,Q0jj_rog,Q1jj_rog,~,QLpjj_rog,D_rog,E_rog,R_rog] = rogersRFA_magW(k_red, Qjj, num_poles);
% save(fullfile(data_dir,'Q0jj_Q1jj_QLpjj_DER_P10x20_p6_kred1_1.mat'),'poles_rog','Q0jj_rog','Q1jj_rog','QLpjj_rog','D_rog','E_rog','R_rog');

load(fullfile(data_dir,'Qjj_npc10nps20_k_red12.mat'))
load(fullfile(data_dir,'Q0jj_Q1jj_QLpjj_DER_P10x20_p6_kred1_1.mat'))

% --- Evaluate RFA Fit Quality --------------------------------------------
[E_max_rog,E_rms_max_rog, E_rms_rog,H_rog] = evalRFA(k_red,Qjj, Q0jj_rog,Q1jj_rog,D_rog,E_rog,R_rog);

legends = {'Qjj(iw)','H rog'};

wi = k_red;
H_hat_cell = {H_rog};
rows = [1,150];
columns = [1,150];


set(0,'defaulttextinterpreter','latex')
title_str = sprintf('Frequenzgang H(i%d < iw < i%d)',min(wi), max(wi));
fig = figure('Name',title_str,'NumberTitle','off');
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
linestyles = {'-','--', '-.', ':'}; % Line styles for each sequence
colors = {'k',[0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], [0.6350, 0.0780, 0.1840]};      % Colors for each sequence


np=length(rows);
nq=length(columns);
for ip = 1:np
    p=rows(ip);
    for iq = 1:nq    
        q=columns(iq);
        
        i = (ip-1)*nq+iq;
        subplot(np,nq,i)
        plot(squeeze(complex(Qjj(p,q,:))),".",'MarkerSize',15)
        hold on
        
        for h = 1:numel(H_hat_cell)
            H_hat = H_hat_cell{h};
            plot(squeeze(complex(H_hat(p,q,:))),'LineStyle', linestyles{h}, 'Color', colors{h},'LineWidth',stdLineWidth)
        end
        axis equal
        grid on
        ax = fig.CurrentAxes; % Get current axes
        ax.TickLabelInterpreter = 'latex'; % Set tick labels to LaTeX
        ax.FontSize = docFontSize;
        set(ax, 'Color', 'w'); % Set axes background to white
        xlabel("Re(z)", 'Interpreter', 'latex', 'FontSize', docFontSize)
        ylabel("Im(z)", 'Interpreter', 'latex', 'FontSize', docFontSize)
        title(['$H(i \omega)$ ',num2str(p),',',num2str(q)], 'Interpreter', 'latex', 'FontSize', docFontSize)
        if ip==1 && iq==2
            legend(legends)
        end

    end
end
    
set(fig, 'Color', 'w'); % Set figure background to white


NyquistPlotName = 'RFA_Nyquist2x2';
fullVpzmap_PlotPath = fullfile(figures_dir, NyquistPlotName);
print(fig, fullVpzmap_PlotPath, '-dpdf', '-vector');
print(fig, fullVpzmap_PlotPath, '-dmeta', '-vector');
