function [E, ele] = build_E_y(Y_tip, Y_root, num_ele)
% build_E_y     Structural beam dimensions and global DoF assignment of
%               bending-torsion beam in y-direction with bending in z-direction
% INPUT
%   Y_tip     : beam tip y-coord
%   Y_root    : beam root y-coord
%   num_ele   : number of equidistant finite elements
%
% OUTPUT
%   E         : Cell array of FEM node positions in structural coords.
%               E{element}{node} = [ x_pos;
%                                    y_pos;
%                                    z_pos]
%   ele       : Element to Global DoF assembly table
%               ele{element}(DoF) = globalDoF
% 
% Sign convention
% _________________________________________________________________________
%   y_s = [0;1;0] is the direction of the beam in structural coord. system
%   z_s = [0;0;1] is the bending motion of the beam in structural coords.
%
%
% Jonas * July 2025
% _________________________________________________________________________
ele = cell(num_ele, 1);
E = cell(num_ele, 1);

Ey_ = linspace(Y_root,Y_tip,num_ele+1);
E1_y_ = Ey_(1:end-1);
E2_y_ = Ey_(2:end);
for e = 1:num_ele
    ele{e} = [1, 2, 3, 4, 5, 6]+3*(e-1);

    E{e}{1} = [0; E1_y_(e); 0];
    E{e}{2} = [0; E2_y_(e); 0];
end
end