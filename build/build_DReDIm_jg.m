function [DRe_jg, DIm_jg] = build_DReDIm_jg(E,ele,Ps)
% build_DReDIm_jg  Assemble the global boundary condition matrice that map
%            beam-element DOFs to normalized downwash at aero panel control points JP.
%
% INPUT
%   E       : cell array, one entry per element.  
%             For element e the nodal y-coordinates are stored as
%                     y1 = E{e}{1}(2),    y2 = E{e}{2}(2).
%
%   ele     : cell array, one entry per element.  
%             ele{e} is a 1 × 6 vector with the global DOF numbers
%             corresponding to element e.
%
%   Ps      : cell array of aerodynamic panels in structural coords
%
% OUTPUT
%   DRe_jg,  
%   DIm_jg  : num_panel × num_DoF_glob
%
%              w_j = DRe_jg           * q_g 
%                  + DIm_jg * 1/V_inf * q_g_dot
%
%              w_j = u_z/V_inf normalized downwash velocity at control point JP 
% 
% Sign convention
% _________________________________________________________________________
%   y_s = [0;1;0] is the direction of the beam in structural coord. system
%   z_s = [0;0;1] is the bending motion of the beam in structural coords.
%
%
% Jonas * July 2025
% _________________________________________________________________________

ne  = length(E);
np  = length(Ps);
ng  = 3 * (ne + 1);          % 3 DOF per node

DRe_jg = zeros(np, ng);
DIm_jg = zeros(np, ng);

for i_e = 1:ne
    E1 = E{i_e}{1};
    E2 = E{i_e}{2};

    [DRe_ele, DIm_ele] = build_DReDIm_ele(E1, E2, Ps);

    id = ele{i_e};                           % global DOF indices of element e
    DRe_jg(:, id) = DRe_jg(:, id) + DRe_ele;
    DIm_jg(:, id) = DIm_jg(:, id) + DIm_ele;
end