function Mgg = build_Mgg(E, ele, rho_bar, I_Tym, I_zm)
% build_Mgg  Assemble the global mass matrix of bending–torsion beam model.
%
%   Each two-node element contributes a 6 × 6 stiffness sub-matrix
%   (three DOF per node:    u_z      : bending deflection  
%                           du_z/dy  : rotation about the x-axis / bending slope
%                           Psi_y    : Torsion angle
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
%   rho_bar : mass per unit length      int_A rho dA      (constant)
%   I_Tym   : polar mass moment         int_A rho r^2 dA   (constant)
%   I_zm    : deviational coupling term int_A rho x  dA   (constant)
%
% OUTPUT
%   Mgg     : global  mass matrix of size 3 × (num_ele+1).
%
% Jonas . July 2025
% _________________________________________________________________________


num_ele      = numel(E);                  % number of elements
num_DoF_glob = 3 * (num_ele + 1);         % 3 DOF per node
Mgg          = zeros(num_DoF_glob);       % use sparse() for large problems

for e = 1:num_ele
    y1   = E{e}{1}(2);
    y2   = E{e}{2}(2);
    M_ele = build_M_ele(y1, y2, rho_bar, I_Tym, I_zm);

    id = ele{e};                         % 1×6 global DOF indices
    Mgg(id, id) = Mgg(id, id) + M_ele;
end
end
