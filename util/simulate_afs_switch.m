function sim = simulate_afs_switch(G,K,opts)
%SIMULATE_AFS_SWITCH  Flutter grows with the controller off, then the controller switches on: one time simulation, frame by frame.
%   sim = simulate_afs_switch(G, K)
%   sim = simulate_afs_switch(G, K, Name=Value)
%   The wing starts in its flutter mode shape at ONE velocity above the flutter
%   speed, diverges open loop for GrowthCycles flutter cycles, then the controller
%   closes the loop at frame i_on and the same model runs on for ControlCycles
%   cycles. Both segments are stepped with expm(A*dt), so the trajectory is an
%   exact solution of the LTI model sampled FramesPerCycle times per flutter
%   cycle. Hand sim.q_f, sim.delta, sim.acc and sim.active to animate_wing.
%
%   G   ss at ONE velocity from build_G_RectWing(V, imuIDX, ailIDX) or
%       build_G_Goland(V, 1, 1): the scaled plant with the PT2 actuators inside,
%       inputs flapN_d / slatN_d (demanded deflection, rad), outputs u_zN_ddot
%       (scaled IMU acceleration), D = 0 (build_ABCD_G drops the apparent-mass
%       feedthrough on purpose). The scaled plant is the right input: q_f and the
%       actuator angles are physical there, and the shipped controllers were
%       tuned on the scaled measurement.
%   K   controller in the u = K*y convention (the same as lft(P,K)): a numeric
%       gain or an ss/tf/zpk, closed as CL = feedback(G, -K, FeedIn, FeedOut).
%       K = [] simulates the open loop only (ControlCycles is ignored).
%
%   Options (defaults)
%   'FeedIn' 1:nu, 'FeedOut' 1:ny   inputs of G the controller drives and outputs
%                    it reads: indices or channel names (e.g. 'flap4_d',
%                    {'u_z4_ddot','u_z8_ddot'})
%   'GrowthCycles' 6, 'ControlCycles' 8, 'FramesPerCycle' 12
%                    dt = 2*pi/omega/FramesPerCycle, omega = imag of the flutter pole
%   'Pluck' 'flutter'   real part of the eigenvector of the least stable
%                    oscillatory pole of G.A, phase-aligned so that the largest
%                    q_f component is real and positive (the actuators start at
%                    rest); or a full state vector x0; or a q_f0 vector
%                    [num_modes x 1] with every other state zero
%   'Amplitude' 0.05    (m) max |q_f(0)|. PHIgf columns are max-normalised, so
%                    q_f1 is about the tip bending displacement
%   'MinImag' 3      (rad/s) poles with |imag| <= MinImag are ignored in the
%                    flutter-mode search and for the closed-loop decay pole
%   'Verbose' true   three console lines: flutter mode, timeline, closed-loop decay
%
%   sim   .t (s), .dt, .Nt, .i_on (GrowthCycles*FramesPerCycle + 1; Nt+1 when
%         K = []), .active (logical 1 x Nt)
%         .q_f [num_modes x Nt] (m); .delta [nu x Nt] (rad) actuator angles in
%         G.InputName order; .u [nu x Nt] (rad) commands, zero before i_on;
%         .acc [ny x Nt] IMU outputs in the scaled units of G
%         .x [n_x + nK x Nt] states of G, then of K; .StateName
%         .lambda, .f_Hz, .sigma_ol   flutter pole; .growthPerCycle = exp(sigma_ol*2*pi/omega)
%         .sigma_cl, .decayPerCycle   least damped oscillatory closed-loop pole (NaN when K = [])
%         .residual   max |q_f| over the last cycle divided by max |q_f| over the
%                     cycle before i_on: what is left after the catch (measured)
%         .CL, .K (as ss), .FeedIn, .FeedOut (resolved indices)
%   Errors: G is an ss array or has D ~= 0; no unstable oscillatory pole ('no
%   flutter at this velocity'); closed loop unstable ('controller does not
%   stabilise this velocity'). States are found by name (q_fN, flapN / slatN);
%   without names the build_G_* layout is assumed and a warning is issued.
%
%   Example, the README hero scene (QUICKSTART.m section 6 and docs/figures/make_readme_figures.m
%   copy these lines):
%     [Structure, ~] = define_RectWing_Structure_Aero(5, 6);
%     G = build_G_RectWing(108, 1:8, 1:8);            % 3.7 m/s above the 104.3 m/s onset
%     load RectWing_Cont_imu18_ail18.mat              % RectWing_Cont: the shipped static gain, 8 IMUs -> 8 surfaces (R01)
%     sim = simulate_afs_switch(G, RectWing_Cont, 'GrowthCycles', 6, 'ControlCycles', 8, 'FramesPerCycle', 12);
%     animate_wing(Structure, sim.q_f, 'Surfaces', sim.delta, 'IMUAccel', sim.acc, 'Active', sim.active, ...
%         'Skin', 'foil', 'NACA', '2408', 'FlapChord', 0.25, 'SlatChord', 0.20, ...
%         'TipAmplitude', 0.07, 'SurfaceAmplitude', deg2rad([33 25]), 'Output', 'figure');
%   Runtime: 0.05 s for 168 frames (model build excluded), measured on R2026a.

arguments
    G
    K
    opts.FeedIn = []                                    % [] -> 1:nu
    opts.FeedOut = []                                   % [] -> 1:ny
    opts.GrowthCycles (1,1) double {mustBePositive,mustBeInteger} = 6
    opts.ControlCycles (1,1) double {mustBeNonnegative,mustBeInteger} = 8
    opts.FramesPerCycle (1,1) double {mustBePositive,mustBeInteger} = 12
    opts.Pluck = 'flutter'
    opts.Amplitude (1,1) double {mustBePositive} = 0.05
    opts.MinImag (1,1) double {mustBeNonnegative} = 3
    opts.Verbose (1,1) logical = true
end

% --- plant checks --------------------------------------------------------
assert(isa(G,'ss'), 'simulate_afs_switch: G must be an ss, got %s.', class(G));
assert(size(G,3) == 1, ['simulate_afs_switch: G must be a single model, not ' ...
    'an ss array (%d models). Pick one velocity, e.g. G(:,:,i_v).'], size(G,3));
assert(all(G.D(:) == 0), ['simulate_afs_switch: G.D must be exactly zero ' ...
    '(max|D| = %g). build_ABCD_G drops the control-surface apparent mass on ' ...
    'purpose so that flapN_d -> u_zN_ddot has no feedthrough; a nonzero D ' ...
    'means this is not a build_G_* plant.'], max(abs(G.D(:))));

n_x = order(G);
nu  = size(G,2);
ny  = size(G,1);

FeedIn  = resolve_channels(opts.FeedIn,  G.InputName,  nu, 'FeedIn',  'input');
FeedOut = resolve_channels(opts.FeedOut, G.OutputName, ny, 'FeedOut', 'output');

% --- state bookkeeping: q_f rows and actuator deflection rows ------------
[qfRows,actRows] = locate_states(G,n_x,nu);
num_modes = numel(qfRows);

% --- flutter mode: least stable oscillatory pole of G.A ------------------
[Phi_G,e_G] = eig(full(G.A));
e_G  = diag(e_G);
cand = find(imag(e_G) > opts.MinImag);          % positive-frequency branch
assert(~isempty(cand) && max(real(e_G(cand))) > 0, ['simulate_afs_switch: ' ...
    'no flutter at this velocity - G.A has no oscillatory pole (imag > %g ' ...
    'rad/s) with a positive real part. Rebuild G above the flutter speed ' ...
    '(RectWing about 104.5 m/s, Goland about 175.6 m/s).'], opts.MinImag);
[~,imax]  = max(real(e_G(cand)));
lambda    = e_G(cand(imax));
v_flutter = Phi_G(:,cand(imax));
omega     = imag(lambda);
sigma_ol  = real(lambda);

fpc = opts.FramesPerCycle;
dt  = 2*pi/omega/fpc;

% --- initial condition ---------------------------------------------------
x0 = zeros(n_x,1);
if (ischar(opts.Pluck) || isstring(opts.Pluck)) && isscalar(string(opts.Pluck))
    assert(strcmpi(opts.Pluck,'flutter'), ['simulate_afs_switch: Pluck must ' ...
        'be ''flutter'', a full state vector [%d x 1] or a modal vector ' ...
        '[%d x 1], got ''%s''.'], n_x, num_modes, char(opts.Pluck));
    % phase align so that the largest q_f entry becomes real and positive,
    % then take the real part (the physical snapshot of the mode)
    [~,imq] = max(abs(v_flutter(qfRows)));
    v_flutter = v_flutter * exp(-1i*angle(v_flutter(qfRows(imq))));
    x0 = real(v_flutter);
else
    assert(isnumeric(opts.Pluck) && isvector(opts.Pluck) && all(isfinite(opts.Pluck)), ...
        ['simulate_afs_switch: Pluck must be ''flutter'' or a finite numeric ' ...
        'vector of length %d (full state) or %d (modal).'], n_x, num_modes);
    p = double(opts.Pluck(:));
    if numel(p) == n_x
        x0 = p;
    elseif numel(p) == num_modes
        x0(qfRows) = p;
    else
        error(['simulate_afs_switch: Pluck vector has length %d; expected %d ' ...
            '(full state x0) or %d (modal q_f0).'], numel(p), n_x, num_modes);
    end
end
qf0max = max(abs(x0(qfRows)));
assert(qf0max > 0, ['simulate_afs_switch: the initial condition has no modal ' ...
    'displacement (all q_f entries zero), so it cannot be scaled to ' ...
    'Amplitude = %g m.'], opts.Amplitude);
x0 = (x0/qf0max) * opts.Amplitude;      % divide first: max|q_f(0)| == Amplitude exactly

% --- controller and closed loop -----------------------------------------
Ng = opts.GrowthCycles*fpc;
if isempty(K)
    Kss = [];  CL = [];  nK = 0;  Nc = 0;
    sigma_cl = NaN;
else
    Kss = to_ss(K,numel(FeedIn),numel(FeedOut));
    CL  = feedback(G,-Kss,FeedIn,FeedOut);      % u = +K*y
    nK  = order(Kss);
    assert(order(CL) == n_x+nK && isequal(CL.StateName(1:n_x),G.StateName), ...
        ['simulate_afs_switch: the closed-loop state vector does not start ' ...
        'with the states of G (feedback reordered the interconnection); the ' ...
        'state bookkeeping below would be wrong.']);
    assert(max(real(eig(full(CL.A)))) < 0, ['simulate_afs_switch: controller ' ...
        'does not stabilise this velocity - max real part of eig(CL.A) = ' ...
        '%+g. Synthesise K for this velocity, or check the u = K*y sign ' ...
        'convention and the FeedIn/FeedOut channels.'], ...
        max(real(eig(full(CL.A)))));
    e_CL  = eig(full(CL.A));
    oscCL = e_CL(abs(imag(e_CL)) > opts.MinImag);
    if isempty(oscCL)
        sigma_cl = NaN;
    else
        sigma_cl = max(real(oscCL));
    end
    Nc = opts.ControlCycles*fpc;
end
Nt   = Ng + Nc;
i_on = Ng + 1;

% --- one continuous simulation: pluck -> diverge -> control -------------
Phi_ol = expm(full(G.A)*dt);                    % controller off (command = 0)
X = zeros(n_x+nK,Nt);
X(1:n_x,1) = x0;
for k = 2:Ng
    X(1:n_x,k) = Phi_ol*X(1:n_x,k-1);
end
if Nc > 0
    Phi_cl = expm(full(CL.A)*dt);               % controller on
    xc = [X(1:n_x,Ng); zeros(nK,1)];            % K starts from rest
    X(:,i_on) = Phi_cl*xc;
    for k = i_on+1:Nt
        X(:,k) = Phi_cl*X(:,k-1);
    end
end

active = false(1,Nt);
active(i_on:end) = true;

% --- outputs -------------------------------------------------------------
acc = G.C*X(1:n_x,:);                           % D = 0
u = zeros(nu,Nt);
if Nc > 0
    u(FeedIn,active) = Kss.C*X(n_x+1:end,active) + Kss.D*acc(FeedOut,active);
end

if isnan(sigma_cl)
    decayPerCycle = NaN;
else
    decayPerCycle = exp(sigma_cl*2*pi/omega);
end
if Nc > 0
    q_peak = max(abs(X(qfRows,max(1,i_on-fpc):i_on-1)),[],'all');
    q_end  = max(abs(X(qfRows,Nt-fpc+1:Nt)),[],'all');
    residual = q_end/q_peak;
else
    residual = NaN;
end

sim = struct();
sim.t      = (0:Nt-1)*dt;
sim.dt     = dt;
sim.Nt     = Nt;
sim.i_on   = i_on;
sim.active = active;
sim.q_f    = X(qfRows,:);
sim.delta  = X(actRows,:);
sim.u      = u;
sim.acc    = acc;
sim.x      = X;
if isempty(CL)
    sim.StateName = G.StateName;
else
    sim.StateName = CL.StateName;
end
sim.lambda         = lambda;
sim.f_Hz           = omega/(2*pi);
sim.sigma_ol       = sigma_ol;
sim.sigma_cl       = sigma_cl;
sim.growthPerCycle = exp(sigma_ol*2*pi/omega);
sim.decayPerCycle  = decayPerCycle;
sim.residual       = residual;
sim.CL             = CL;
sim.K              = Kss;
sim.FeedIn         = FeedIn;
sim.FeedOut        = FeedOut;

% --- console summary -----------------------------------------------------
if opts.Verbose
    fprintf('flutter mode: f = %.2f Hz, sigma = %+.3f 1/s, growth x%.3f per cycle\n', ...
        sim.f_Hz, sigma_ol, sim.growthPerCycle);
    if Nc > 0
        fprintf(['open loop %d cycles (%d frames), controller on at frame %d, ' ...
            'closed loop %d cycles (%d frames), dt = %.1f ms\n'], ...
            opts.GrowthCycles, Ng, i_on, opts.ControlCycles, Nc, dt*1e3);
        fprintf(['closed loop: least-damped oscillatory pole sigma = %+.3f 1/s, ' ...
            'decay x%.3f per cycle; residual after %d cycles: %.0f %% of the peak\n'], ...
            sigma_cl, decayPerCycle, opts.ControlCycles, 100*residual);
    else
        fprintf('open loop %d cycles (%d frames), no controller, dt = %.1f ms\n', ...
            opts.GrowthCycles, Ng, dt*1e3);
    end
end
end

% =========================================================================
function idx = resolve_channels(sel,names,n,optName,what)
% Channel selection by index or by name, always returned as a row of indices.
if isempty(sel) && ~ischar(sel)
    idx = 1:n;
    return
end
if isnumeric(sel) || islogical(sel)
    if islogical(sel)
        assert(numel(sel) == n, ['simulate_afs_switch: logical %s must have ' ...
            '%d entries (one per %s of G).'], optName, n, what);
        idx = find(sel(:)');
    else
        idx = double(sel(:)');
    end
    assert(all(idx == round(idx)) && all(idx >= 1 & idx <= n), ...
        ['simulate_afs_switch: %s indices must be integers in 1..%d ' ...
        '(G has %d %ss).'], optName, n, n, what);
else
    assert(ischar(sel) || isstring(sel) || iscellstr(sel), ...
        ['simulate_afs_switch: %s must be indices or channel names ' ...
        '(char, string or cellstr).'], optName);
    sel = cellstr(sel);
    idx = zeros(1,numel(sel));
    for i = 1:numel(sel)
        hit = find(strcmp(names,sel{i}));
        assert(~isempty(hit), ['simulate_afs_switch: %s channel ''%s'' is not ' ...
            'an %s name of G. Available: %s.'], optName, sel{i}, what, ...
            strjoin(names(:)',', '));
        assert(isscalar(hit), ['simulate_afs_switch: %s channel ''%s'' is not ' ...
            'unique in G (%d matches).'], optName, sel{i}, numel(hit));
        idx(i) = hit;
    end
end
assert(numel(unique(idx)) == numel(idx), ...
    'simulate_afs_switch: %s lists the same channel twice.', optName);
end

% =========================================================================
function [qfRows,actRows] = locate_states(G,n_x,nu)
% Row indices of the modal displacements q_f1..q_fN and of the actuator
% deflection states, found by name; positions are the documented fallback.
sn = G.StateName;
if isempty(sn)
    sn = repmat({''},n_x,1);
end
in = G.InputName;
if isempty(in)
    in = repmat({''},nu,1);
end

hit = ~cellfun(@isempty,regexp(sn,'^q_f\d+$','once'));
qfRows = find(hit(:)');
[~,ord] = sort(cellfun(@(s) str2double(s(4:end)),sn(qfRows)));
qfRows  = qfRows(ord);                          % q_f1, q_f2, ... numeric order

actRows = zeros(1,nu);
for j = 1:nu
    nm = in{j};
    if numel(nm) > 2 && strcmp(nm(end-1:end),'_d')
        h = find(strcmp(sn,nm(1:end-2)));
        if isscalar(h)
            actRows(j) = h;
        end
    end
end

if isempty(qfRows) || any(actRows == 0) || numel(unique(actRows)) < nu
    % unnamed or renamed model: fall back to the build_G_* state layout,
    % x = [q_f(1:nm); q_f_dot(1:nm); aero_lag(nm*6); (angle,rate) per surface]
    num_modes = (n_x - 2*nu)/8;                 % num_poles = 6
    assert(num_modes == round(num_modes) && num_modes > 0, ...
        ['simulate_afs_switch: cannot locate the q_f and actuator states of ' ...
        'G. The state names do not match build_G_* (q_f1.., flapN/slatN) and ' ...
        'the order %d is not 8*num_modes + 2*%d either.'], n_x, nu);
    warning('simulate_afs_switch:stateNames', ...
        ['G has no usable q_f / actuator state names; assuming the ' ...
        'build_G_* layout: q_f = states 1:%d, actuator deflections = states ' ...
        '%d:2:%d (num_poles = 6).'], num_modes, n_x-2*nu+1, n_x-1);
    qfRows  = 1:num_modes;
    actRows = (n_x-2*nu) + 2*(0:nu-1) + 1;      % (angle,rate) pairs, in input order
end
end

% =========================================================================
function Kss = to_ss(K,nu_K,ny_K)
% Controller as an ss, with a dimension check against FeedIn/FeedOut.
if isnumeric(K)
    assert(ismatrix(K) && all(isfinite(K(:))), ...
        'simulate_afs_switch: a numeric K must be a finite [%d x %d] matrix.', nu_K, ny_K);
    Kss = ss(double(K));
else
    assert(isa(K,'DynamicSystem'), ['simulate_afs_switch: K must be a numeric ' ...
        'gain or an ss/tf/zpk model, got %s.'], class(K));
    assert(size(K,3) == 1, ['simulate_afs_switch: K must be a single model, ' ...
        'not an array (%d models).'], size(K,3));
    Kss = ss(K);
    assert(Kss.Ts == 0, ['simulate_afs_switch: K must be continuous time ' ...
        '(Ts = %g); this simulation steps the continuous model with expm.'], Kss.Ts);
end
assert(size(Kss,1) == nu_K && size(Kss,2) == ny_K, ['simulate_afs_switch: K ' ...
    'is [%d x %d] but u = K*y needs [%d x %d]: %d FeedIn channel(s) and %d ' ...
    'FeedOut channel(s).'], size(Kss,1), size(Kss,2), nu_K, ny_K, nu_K, ny_K);
end
