function Kgg = build_Kgg(E, ele, EI_xxa, GI_Tya)
% build_Kgg  Assemble the global stiffness matrix of bending–torsion beam model.
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
%   EI_xxa  : bending stiffness  (E I_xxa)   I_xxa = int_A z^2 dA (constant)
%   GI_Tya  : torsional stiffness (G I_Tya)  I_Tya = int_A r^2 dA (constant)
%
% OUTPUT
%   Kgg     : global stiffness matrix (size 3 × (num_ele+1)).
%
% Jonas . July 2025
% _________________________________________________________________________

num_ele       = numel(E);              % number of elements
num_DoF_glob  = 3 * (num_ele + 1);     % total global DOF
Kgg           = zeros(num_DoF_glob);   % initialise (use sparse() if desired)

for e = 1:num_ele
    y1    = E{e}{1}(2);
    y2    = E{e}{2}(2);
    K_ele = build_K_ele(y1, y2, EI_xxa, GI_Tya);

    id = ele{e};                       % 1×6 global DOF indices
    Kgg(id, id) = Kgg(id, id) + K_ele;
end
end
