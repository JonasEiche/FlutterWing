function [DRe_jx,DIm_jx] = build_DReDIm_jx(cspanels,rot_axis,Pa)
% build_DReDIm_jx       Maps control surface deflections and deflection
%                       rates to aerodynamic panel boundary conditions on 
%                       downwash at control points JP
%
% INPUT
%   cspanels : panel numbers belonging to control surface cs
%              cspanels{cs}      = [p1,p2,...pcsn]
%   rot_axis : axis of rotation defined as line from (1) to (2)
%              rot_axis{cs}{1,2} = [x;y;z]_aero
%   Pa       : cell array, one entry per aerodynamic panel in aerodynamic
%              coordinate system with freestream in positive x-dir.
%              For each panel the four corner nodes are stored as
% 
%               P{panel}{node} = [ x_pos;         1---4
%                                  y_pos;         |   |
%                                  z_pos ]        2---3
% 
% OUTPUT
%   DRe_jx,        
%   DIm_jx        : [num_panel, num_cs] Matrix
%
%                   w_j = DRe_jx * u_x + 1/V_inf * DIm_jx * u_x_dot
% 
% Sign convention
% _________________________________________________________________________
% Freestream direction is positiv x-axis in aero coord. sys.
%
% Jonas * July 2025
% _________________________________________________________________________

np=length(Pa);
ncs=length(cspanels);

DRe_jx = zeros(np, ncs);
DIm_jx = zeros(np, ncs);
% HM_xj = zeros(ncs,np);      % Hinge Moment HMcs = q_bar * HM_xj * Cp_j
for i_cs = 1:ncs
    RP1 = rot_axis{i_cs}{1};
    RP2 = rot_axis{i_cs}{2};
    rot=(RP2-RP1);
    rot=rot/norm(rot);

    for i_p = cspanels{i_cs}
        JP1 = 0.75*Pa{i_p}{2} + 0.25*Pa{i_p}{1};
        JP2 = 0.75*Pa{i_p}{3} + 0.25*Pa{i_p}{4};
        JP  = 0.5*(JP1+JP2);
        N=cross((Pa{i_p}{2} - Pa{i_p}{1}), (Pa{i_p}{4}-Pa{i_p}{1}));
        N = N/norm(N);
        % S = [0;Pa{i_p}{4}(2:3) - Pa{i_p}{1}(2:3)];
        S = Pa{i_p}{4}-Pa{i_p}{1} - (Pa{i_p}{2}-Pa{i_p}{1})/norm((Pa{i_p}{2}-Pa{i_p}{1}))^2*dot((Pa{i_p}{2}-Pa{i_p}{1}),(Pa{i_p}{4}-Pa{i_p}{1}));
        s = norm(S);
        S = S/s;

        % c = norm(0.5*(Pa{i_p}{2}-Pa{i_p}{1}+Pa{i_p}{3}-Pa{i_p}{4}));
        % a = s*c;
        % BVP1 = 0.25 * Pa{i_p}{2} + 0.75 * Pa{i_p}{1};
        % BVP2 = 0.25 * Pa{i_p}{3} + 0.75 * Pa{i_p}{4};
        % LP   = 0.5  * (BVP1 + BVP2);
        % HM_xj(i_cs,i_p) = a * dot( rot, cross(LP-RP1,N) );


        DRe_jx(i_p,i_cs) = dot(rot, S);
        DIm_jx(i_p,i_cs) = dot(-N, cross(rot,JP-RP1));
    end
end

end









