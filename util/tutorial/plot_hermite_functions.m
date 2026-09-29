function fig = plot_hermite_functions(opts)
%PLOT_HERMITE_FUNCTIONS  The four cubic Hermite shape functions of the beam bending deflection (TUTORIAL.m section (2)).
%   fig = plot_hermite_functions()
%   fig = plot_hermite_functions('Name', 'fem_shape_functions')
%   Draws H_1 ... H_4 over the natural coordinate eps of one element, -1 <= eps <= 1: the
%   deflection functions of node 1 and node 2 and the two slope functions, drawn for an element
%   length of 2. The two torsion functions are the linear ones, 0.5*(1 -/+ eps).
%   Options (defaults)
%     'Name'   'Hermite shape functions'   figure name (the tutorial export uses the file stem)
%   fig      fw_figure of the wide size
%   See also beamPHI, build_K_ele, build_M_ele, fw_figure.

arguments
    opts.Name (1,:) char = 'Hermite shape functions'
end

S = fw_style();
eps_e = linspace(-1,1,201);
L_e   = 2;                                           % element length used for the two slope functions
H_shape = [0.25*(2 - 3*eps_e + eps_e.^3);            % H_1: nodal deflection u_z1
           0.125*L_e*(1 - eps_e - eps_e.^2 + eps_e.^3);   % H_2: nodal slope u_z1'
           0.25*(2 + 3*eps_e - eps_e.^3);            % H_3: nodal deflection u_z2
           0.125*L_e*(-1 - eps_e + eps_e.^2 + eps_e.^3)]; % H_4: nodal slope u_z2'

fig = fw_figure(S.size.wide(1), S.size.wide(2), 'Name', opts.Name);
ax  = axes(fig); hold(ax,'on')
plot(ax, eps_e, H_shape, 'LineWidth', S.lineWidth)
xlim(ax,[-1 1]); ylim(ax,[-0.80 1.05])              % room for the legend below the curves
xlabel(ax,'$\varepsilon$'); ylabel(ax,'shape function'); title(ax,'Cubic Hermite shape functions of the bending deflection')
legend(ax, {'$H_1$ ($u_{z1}$)','$H_2$ ($u_{z1}''$)','$H_3$ ($u_{z2}$)','$H_4$ ($u_{z2}''$)'}, ...
    'Location','south', 'NumColumns',2)
end
