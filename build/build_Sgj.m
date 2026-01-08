function Sgj = build_Sgj(E,ele,Ps)
% build_Sgj  Assemble the global “splining’’ matrix that maps
%            aerodynamic panel forces to beam-element DOFs.
%
%            Fg_aero  =  q_bar  *  Sgj  *  DCp_j
%
%            where  q_bar   : dynamic pressure
%                   DCp_j   : vector of panel pressure coefficients
%                   Fg_aero : vector of nodal forces / moments on FEM DoFs
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
%   Sgj     : global mapping of aero panel pressure coeffs onto FEM 
%
% Jonas . July 2025
% _________________________________________________________________________

ne  = length(E);             % number of elements
np  = length(Ps);            % number of aerodynamic panels
ng  = 3 * (ne + 1);          % total structural DOFs

Sgj = zeros(ng, np);         % (use sparse for very large cases)

for i_e = 1:ne
    E1 = E{i_e}{1};
    E2 = E{i_e}{2};
    % --- element splining matrix, 6 × num_panel
    S_ele = build_S_ele(E1,E2,Ps);
    id = ele{i_e};
    Sgj(id, :) = Sgj(id, :) + S_ele;
end