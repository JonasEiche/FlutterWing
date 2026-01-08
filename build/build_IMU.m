function [PHIzg, PHItg] = build_IMU(x_pos, y_pos, E, ele)
% build_IMU   Maps structural FEM DoF accelerations to measured
%             acceleration at position x_pos, y_pos on structural cs.
%             Maps structural FEM DoF velocities to measured rotation rates
%             
%             Kinematic relation (small-angle assumption)
%             u_z_ddot  =  u_z_ddot(y) - x_pos * Psi_y_ddot(y) -->  PHIzg
%             Psi_y_dt  =  Psi_y_dot(y)                        -->  PHItg
%
% INPUT
%   x_pos   : vector with the IMU x-coordinates in struct. cs.  [n_IMU × 1]
%   y_pos   : vector with the IMU y-coordinates in struct. cs.  [n_IMU × 1]
%
%   E       : cell array, one entry per element.  
%             For element e the nodal y-coordinates are stored as
%                     y1 = E{e}{1}(2),    y2 = E{e}{2}(2).
%
%   ele     : cell array, one entry per element.  
%             ele{e} is a 1 × 6 vector with the global DOF numbers
%             corresponding to element e.
% OUTPUT
%   PHIzg   : n_IMU × n_DOF   – maps  q_g_ddot  to  u_z_ddot
%   PHItg   : n_IMU × n_DOF   – maps  q_g_dot   to  Psi_y_dt
% 
% Sign convention
% _________________________________________________________________________
%   y_s = [0;1;0] is the direction of the beam in structural coord. system
%   z_s = [0;0;1] is the bending motion of the beam in structural coords.
%
%
% Jonas * July 2025
% _________________________________________________________________________

% --- 0. Ensure column vectors -------------------------------------------
x_pos = x_pos(:);      % n × 1
y_pos = y_pos(:);

% --- 1. Shape-function matrices for u_z and Psi_y -----------------------
[PHIuzg, PHIpsiyg] = beamPHI(y_pos, E, ele);   % w(y) , ψ(y)

% --- 2. IMU kinematic coupling ------------------------------------------
PHIzg = PHIuzg - x_pos .* PHIpsiyg;
PHItg = PHIpsiyg;
end
