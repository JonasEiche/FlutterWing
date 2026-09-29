function [EV, MS] = pkmethode(Structure,Aero, V_inf)
%PKMETHODE  p-k flutter analysis with a fresh DLM solve per iteration (slow).
%   [EV,MS] = pkmethode(Structure,Aero,V_inf)
%   Structure  needs Kff, Mff, Dff [n x n], Sfj [n x nj], DRe_jf, DIm_jf
%              [nj x n] and the DLM panel corners Pa, as set by
%              define_*_Structure_Aero or by sections (1)-(5) of
%              TUTORIAL.m; n = num_modes is taken from size(Kff,1)
%   Aero       needs c_ref (m) and rho (kg/m^3)
%   V_inf      [1 x n_vel] freestream velocities (m/s)
%   EV         [2*n x n_vel] eigenvalues (rad/s); row i is one mode followed
%              across velocities, a conjugate pair occupies two rows. Row
%              order = wind-off (k_red = 0) eigenvalues sorted by magnitude.
%   MS         [2*n x 2*n x n_vel] right eigenvectors of the first-order form
%              with state [q_f; q_f_dot]; MS(:,i,iv) belongs to EV(i,iv)
%   Fixed inside: Ma = 0 and SYM = 1 (must match the define_* files) and the
%   k_red convergence tolerance tol = 1e-3 (absolute).
%   For every mode and velocity the reduced frequency is iterated to the
%   fixed point k_red = |Im(lambda)|*c_ref/(2*V_inf). Every iterate calls
%   build_Qjj, i.e. one DLM solve: the Goland wing on 11 velocities takes
%   about 9 s (R2026a, 2025 laptop), the RectWing minutes.
%   The while loop has no iteration cap: if k_red does not converge for a
%   mode the call does not return. Mode identity between iterates and
%   velocities is kept with eigenshuffle (util/eigenshuffle).
%   Example (Goland wing, TUTORIAL.m section (7)):
%     V_inf = linspace(90,190,11);
%     [EV,MS] = pkmethode(Structure,Aero,V_inf);
%     Vg_plot(V_inf, EV, 'open loop, p-k method', 'Band', [3 20]);

nv = length(V_inf);
Ma = 0.0;
SYM =1;


Kff = Structure.Kff;
Mff = Structure.Mff;
Dff = Structure.Dff;
Sfj = Structure.Sfj;
DRe_jf = Structure.DRe_jf;
DIm_jf = Structure.DIm_jf;

Pa = Structure.Pa;
c_ref = Aero.c_ref;
rho = Aero.rho;
num_modes = size(Kff,1);
tol = 0.001;

Q_windoff = [zeros(num_modes, num_modes) eye(num_modes, num_modes); -Mff\(Kff) -Mff\(Dff)];
[modeshapes, evals] = eig(Q_windoff);    % wind off: k_red=0   -->  omega_windoff = imag(evals)
evals = diag(evals);
[evals, idx] = sort(evals);
modeshapes = modeshapes(:,idx);
dim_A = length(evals);              % dim_A = 2*num_modes

% MS : [2*num_modes,2*num_modes,num_vel] eigenvector x mode x velocity
% EV : [2*num_modes,num_vel]             eigenvalue for every mode x velocity
% NK : [2*num_modes,num_vel]             number of k_red iterations per mode x velocity (internal, not returned)

MS_i_v = repmat(modeshapes,[1 1 dim_A]);
EV_i_v = repmat(evals,[1 dim_A]);
MS = zeros(dim_A,dim_A,nv);            % complex conjugate modes are tracked as two separate rows
EV = zeros(dim_A,nv);
NK = zeros(dim_A,nv);
for i_v = 1:nv
V_inf_i = V_inf(i_v);
disp(['PK Iteration V_inf:  ',num2str(V_inf_i)])
    for i_m = 1:dim_A
        k_red = 0;

        modeshapes = MS_i_v(:,:,i_m);
        evals = EV_i_v(:,i_m);

%         alternative relative-tolerance criterion, not used:
%         if k_red<1.0;
%         flag = abs(abs(imag(evals(i_m)))*c_ref/V_inf-k_red) > tol;
%         elseif k_red>=1.0;
%         flag = abs(abs(imag(evals(i_m)))*c_ref/V_inf-k_red) > tol*k_red;
%         end
        nk = 0;
        while abs(abs(imag(evals(i_m)))*0.5*c_ref/V_inf_i-k_red) > tol
            nk=nk+1;
            k_red = abs(imag(evals(i_m)))*0.5*c_ref/V_inf_i;              % abs() keeps k_red >= 0; k_red = 0 (static mode) breaks below
            if k_red == 0, break, end   % static/divergence: skip oscillatory aero update (avoids /omega=0), keep last valid eig
            Qjj = build_Qjj(Ma,k_red,c_ref,Pa,SYM);
            QjjRe = real(Qjj);
            QjjIm = imag(Qjj);

            q_bar = 0.5*rho*V_inf_i^2;
            omega = 2*V_inf_i/c_ref*k_red;
            QggRe_til = q_bar*Sfj*(QjjRe*DRe_jf - QjjIm/V_inf_i*DIm_jf*omega);
            QggIm_til = q_bar*Sfj*(QjjIm*DRe_jf/omega+QjjRe/V_inf_i*DIm_jf);
            A = [zeros(num_modes, num_modes)     ,       eye(num_modes, num_modes); 
                 Mff\(QggRe_til-Kff) ,       Mff\(QggIm_til-Dff)];

            [modeshapes1, evals1] = eig(A);
            evals1 = diag(evals1);
            idx = eigenshuffle(evals, modeshapes, evals1, modeshapes1); % mode tracking to chose the eigenvalue belonging to the current mode
            modeshapes = modeshapes1(:,idx);
            evals = evals1(idx);
        end
        MS_i_v(:,:,i_m) = modeshapes;
        EV_i_v(:,i_m) = evals;
        MS(:,i_m,i_v) = modeshapes(:,i_m);
        EV(i_m,i_v) = evals(i_m);
        NK(i_m,i_v) = nk;
    end
end

end
