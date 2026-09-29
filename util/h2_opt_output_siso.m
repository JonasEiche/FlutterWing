function [M,h] = h2_opt_output_siso(C,v)
%H2_OPT_OUTPUT_SISO  Real sensor combination maximising the output energy of one complex pole.
%
%   [M,h] = h2_opt_output_siso(C,v)
%   NOTE: the direction h is the SECOND output; the MIMO variants
%   h2_opt_input_mimo / h2_opt_output_mimo return the direction FIRST.
%
%   Inputs
%   ------
%   C : (n_y x n)   real     – plant output matrix  (y_sensor = C x)
%   v : (n   x 1)   complex  – right eigenvector
%
%   Outputs
%   -------
%   M : (n_y x n_y) real  – rank-2 matrix  M = a a^T + b b^T
%   h : (n_y x 1)   real  – unit vector maximising |h^T C v|_2
%
%   Derivation
%   ----------
%     z = C v  (complex gains from mode --> each sensor)
%     a = Re(z),  b = Im(z)
%     |h^T C v|_2 = h^T (a a^T + b b^T) h  → Rayleigh quotient.
%
%



% ---------- Step 1: mode → sensor coupling vector -----------------------
z = C * v;                 % n_y x 1  complex column
% ---------- Step 2: real/imaginary split --------------------------------
a = real(z);               % n_y x 1
b = imag(z);               % n_y x 1
% ---------- Step 3: build the  rank-2  observation matrix ---------------
M = a*a.' + b*b.';         % real, symmetric
% ---------- Step 4: principal eigenvector (optimal sensor combo) --------
if size(M,1) == 1          % only one physical sensor
    h = 1;                 % trivial: that sensor is the only choice
else
    [vec,~] = eigs(M,1,'la');   % largest eigenvalue
    h = real(vec);              % remove negligible round-off
    h = h / norm(h);            % normalise
end
end
