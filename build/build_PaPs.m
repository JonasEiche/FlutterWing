function [Pa,Ps] = build_PaPs(s,c,np_s,np_c,sw,dh,tr,xf)
% build_PaPs      Setup aerodynamic and structural panel cell array for a
%                 planar (half-)wing with sweep, dihedral and taper. The 
%                 aerodynamic coordinate system x-axis is flow direction and 
%                 the xz-plane is the symmetry plane between left and right wing.
%               
%                 The structural coordinate system y-axis coincides with the wing
%                 flexural axis and x-axis is in plane with the wing.
% 
% INPUT
%   s           : semi span
%   c           : root chord
%   np_s        : number panels spanwise
%   np_c        : number panels chordwise
%   sw (deg)    : sweep angle at leading edge
%   dh (deg)    : dihedral
%   tr          : taper ratio
%   xf          : flex axis position relative to chord
%
% OUTPUT
%   Pa          : cell array, one entry per aerodynamic panel in aerodynamic/structural
%                 coordinate system with freestream in positive x-dir.
%                 For each panel the four corner nodes are stored as
% 
%                 P{panel}{node} = [ x_pos;         1---4
%                                    y_pos;         |   |
%                                    z_pos ]        2---3
%
%                 Panels are numbered rowwise against flow direction:
%                   5   6   7   8
%                   1   2   3   4


% CURRENTLY NOT ------------------------------------------------
%                 Panels are numbered rowwise in flow direction:
%                   1   2   3   4
%                   5   6   7   8
% CURRENTLY NOT ------------------------------------------------


%
% Jonas * July 2025
% _________________________________________________________________________
Y = linspace(0,s,np_s+1);
Z = linspace(0,tand(dh)*s,np_s+1);
C_y = linspace(1,tr,np_s+1)*c;
X_le = Y*tand(sw);
X_te = X_le+C_y;
X = zeros(np_c+1,np_s+1);
for i_y=1:np_s+1
    x_loc = linspace(X_le(i_y),X_te(i_y),np_c+1);
    X(:,i_y) = x_loc;
end
Y = repmat(Y,np_c+1,1);
Z = repmat(Z,np_c+1,1);

Pa = cell(np_c*np_s,1);
for i_x = 1:np_c
    for i_y = 1:np_s
        % i_p = (i_x-1)*np_s+i_y;
        i_p = (np_c-i_x)*np_s+i_y;
        Pa{i_p}{1} = [X(i_x,i_y);Y(i_x,i_y);Z(i_x,i_y)];
        Pa{i_p}{2} = [X(i_x+1,i_y);Y(i_x+1,i_y);Z(i_x+1,i_y)];
        Pa{i_p}{3} = [X(i_x+1,i_y+1);Y(i_x+1,i_y+1);Z(i_x+1,i_y+1)];
        Pa{i_p}{4} = [X(i_x,i_y+1);Y(i_x,i_y+1);Z(i_x,i_y+1)];
    end
end

z_s = [0;sind(dh);-cosd(dh)];
y_s = [tand(sw)*s+xf*c*(tr-1);s;tand(dh)*s];
s_s = norm(y_s);
y_s = y_s/norm(y_s);
x_s = cross(y_s, z_s);
x_s = x_s / norm(x_s);

R_sa = [x_s,y_s,z_s]'; % Rotation into Wing Plane

X = X-xf*c; % x_offset

x_p = X(:);
y_p = Y(:);
z_p = Z(:);
points = [x_p, y_p, z_p];
points_s = (R_sa * points')';
X_s = reshape(points_s(:,1), size(X));
Y_s = reshape(points_s(:,2), size(Y));
Z_s = reshape(points_s(:,3), size(Z));

assert(all(Z_s < 1e-12, 'all'),'R_sa transformiert nicht in Ebene!')
Z_s = zeros(size(Z_s));

Ps = cell(np_c*np_s,1);
for i_x = 1:np_c
    for i_y = 1:np_s
        % i_p = (i_x-1)*np_s+i_y;
        i_p = (np_c-i_x)*np_s+i_y;
        Ps{i_p}{1} = [X_s(i_x,i_y);Y_s(i_x,i_y);Z_s(i_x,i_y)];
        Ps{i_p}{2} = [X_s(i_x+1,i_y);Y_s(i_x+1,i_y);Z_s(i_x+1,i_y)];
        Ps{i_p}{3} = [X_s(i_x+1,i_y+1);Y_s(i_x+1,i_y+1);Z_s(i_x+1,i_y+1)];
        Ps{i_p}{4} = [X_s(i_x,i_y+1);Y_s(i_x,i_y+1);Z_s(i_x,i_y+1)];
    end
end


end