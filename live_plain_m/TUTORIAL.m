%% The Goland wing from beam element to virtual flight
% Build the Goland wing from beam elements to a controlled aeroelastic model, with equations and
% experiments. Each section (n) belongs to chapter n of |docs/TUTORIAL.md|, which has the
% figures, equations and derivations.
%
% To run it, use MATLAB R2023a or newer with Control System Toolbox. Set the Current Folder to
% the repository root and run |startup.m|. Section (0) is explanatory; then run sections (1) to
% (12) in order. See |docs/CONVENTIONS.md| for the coordinate frames, signs, numbering and
% scaling.

clearvars

%% (0) What is flutter and why this repository
% The textbook definition of flutter: Flutter denotes a self-excited aeroelastic instability
% arising from the interaction of unsteady aerodynamics with structural modes, often involving
% coupling between bending and torsional dynamics.
%
% Flutter occurs when two or more structural modes couple through the aerodynamic forces induced
% by their oscillation. When the net work produced by the airflow over one structural
% oscillatory cycle becomes positive, i.e. the airflow pumps energy into the structure, the
% oscillation grows unbounded. On the Goland wing, the benchmark case treated in this tutorial,
% this happens at about 165 m/s. Bending and twist couple: twist changes the lift, and lift
% drives bending.
%
% As modern aircraft design trends increasingly favor lightweight and flexible structures, the
% susceptibility to the flutter phenomenon emerges as a fundamental constraint. Conventional
% approaches mitigate this risk through conservative structural design. However,active control
% techniques offer the potential to extend operational boundaries by artificially increasing
% modal damping.
%
% **The V-g plot**
%
% Everything this tutorial computes ends up in one type of plot. The linear model of the wing at
% one airspeed $V_\infty$ has eigenvalues $\lambda = \sigma + \mathrm{i}\omega$, one pair per
% retained structural mode. Each pair is drawn as a frequency and a damping,
%
% $$f = \frac{|\lambda|}{2\pi}, \qquad g = -100\,\frac{\mathrm{Re}\,\lambda}{|\lambda|} ,$$
%
% which is what |Vg_plot| draws: $f$ in Hz on top, $g$ in percent below (the damping ratio
% $\zeta$ in percent, positive when the motion decays). Displayed over the relevant range of
% airspeeds. Here only the bending and torsion modes are displayed for clarity. Read the plot
% from left to right. At low airspeed the two modes sit at their wind-off frequencies, 7.67 Hz
% for bending and 15.27 Hz for torsion. As $V_\infty$ grows, the aerodynamic forces pull the
% torsion frequency down and push the bending frequency up. Below the flutter speed the air
% takes energy out of the motion, at the flutter speed the system is marginally stable.
%
% **Why this repository exists**
%
% When I started my PhD in active control of aeroelastic systems some years ago, I struggled to
% find a plant model to try my controller synthesis ideas on. The industrial models I
% encountered had hundreds of poles and layers of SIMULINK blocks, backed by NASTRAN models
% tuned against CFD and test data. I could work with their linearized $A$, $B$, $C$, $D$
% matrices, but following each matrix back to the physics was difficult.
%
% Two-dimensional airfoil models were easier to inspect, but the results typically do not
% generalize to the industrial setting. Hence I (somewhat foolishly) decided to build my own
% model: simple enough to follow along and understand the physics, yet with the same modeling
% stages used in industrial preliminary loads analysis.
%
% The result is this repository. The code is self-contained MATLAB, with short functions that
% can be inspected and modified individually.
%
% **What the model assumes**
%
% The results depend on the following assumptions.
%
% - **A slender wing as a beam.** The structure is a one-dimensional cantilever beam along the
%   flexural axis with bending and torsion. Chordwise deformation does not exist.
% - **Flat-plate potential flow.** The unsteady air loads come from the planar doublet-lattice
%   method (DLM), a frequency-domain solution of the linearised potential-flow equations: thin
%   flat lifting surface, small motions, inviscid and irrotational flow. The Mach number is
%   fixed per build and is zero (incompressible flow) for both wings of this repository.
% - **A few modes.** The structure is reduced to its first five vibration modes before any
%   aerodynamics is attached; flutter of a clean wing lives in the lowest ones.
% - **Coupling by virtual work.** Panel forces reach the beam, and beam motion reaches the
%   panels, through the beam's own shape functions along a rigid bar perpendicular to the beam.
%   No surface spline, no interpolation layer.
% - **Linear and time-invariant per airspeed.** The result is one state-space model per
%   freestream velocity, stored as an |ss| array, or the same family as a velocity-scheduled
%   |lpvss| for simulation.
% - **Second-order actuators.** Each control surface is driven by a second-order servo model
%   without deflection or rate limits.
%
% **Source:** the whole pipeline is |define/define_Goland_Structure_Aero.m| feeding
% |build/build_G_Goland.m|.

%% (1) The Goland wing in numbers
% The Goland wing is the textbook flutter case: a rectangular, untwisted, uniform cantilever
% half-wing published by Goland in 1945 with an analytical flutter speed. Its geometry, its
% section properties per unit length and the discretization are all in
% |define/define_Goland_Structure_Aero.m|. You may change any one of them and watch the V-g plot
% react.
%
% - $\rho$ = 1.020 kg/m^3: air density
% - $s$ = 6.096 m: semi span, root to tip
% - $c$ = 1.8288 m: chord, also the reference chord $c_\mathrm{ref}$ of the reduced frequency
% - $x_f$ = 0.33: flexural axis position, as a fraction of the chord from the leading edge
% - $\zeta$ = 0: structural damping ratio: the benchmark is undamped
% - $x_m$ = 0.43: mass axis position, chord fraction from the leading edge
% - $\bar\rho$ = 35.71 kg/m: mass per unit length, $\int_A \rho\,\mathrm{d}A$
% - $I_{Tym}$ = 8.64 kg m: polar mass moment of inertia per unit length about the flexural axis,
%   $\int_A \rho r^2\,\mathrm{d}A$
% - $EI$ = 9.77 MN m^2: flexural rigidity
% - $GJ$ = 0.99 MN m^2: torsional rigidity
%
% One derived quantity couples bending and torsion in the mass matrix: the static mass moment
% per unit length about the flexural axis,
%
% $$I_{zm} = \int_A \rho\, x\,\mathrm{d}A = -\bar\rho\,(x_m - x_f)\,c = -6.53\ \mathrm{kg},$$
%
% negative because the structural $x$ axis points toward the leading edge and the mass axis lies
% behind the flexural axis. Move the mass axis forward of the flexural axis and $I_{zm}$ changes
% sign; chapter 6 does exactly that.
%
% The Doublet Lattice Method is used to calculate complex valued Aerodynamic Influence
% Coefficient Matrices with $n_c \times n_s$ = 5 x 10: doublet-lattice panels. For the beam FEM
% only $n_{ele}$=10 beam finite elements along the span are sufficient to model the structure of
% the wing. The unsteady aerodynamics are computed at 14 reduced frequencies from 0 to 1.1
% (chapter 3). A reduced frequency $k$ at a velocity $V_\infty$ is a physical frequency
% $f = k \cdot 2V_\infty/(2\pi c_\mathrm{ref})$, so the largest sample covers 19.1 Hz at 100 m/s
% and 36 Hz at 190 m/s, comfortably above the two modes that flutter.
%
% Two coordinate frames appear in the code, and both are drawn at the root of the wing below.
% The aerodynamic frame |Pa| has its origin at the root leading edge, $x_a$ pointing downstream,
% $y_a$ along the span and $z_a$ up. The doublet-lattice solver and the control-surface hinges
% live in it. The structural frame |Ps| has its origin on the flexural axis at the root, $x_s$
% pointing toward the leading edge, $y_s$ along the beam and $z_s$ pointing down, as beam theory
% likes it. The finite elements, the coupling and the accelerometer live in it. For an unswept
% wing the two are related by $x_s = x_f c - x_a$, $y_s = y_a$, $z_s = -z_a$.
%
% **Source:** |define/define_Goland_Structure_Aero.m|; |build/build_PaPs.m| (the panel corners
% in both frames).

rho = 1.020;                % air density  (kg/m^3)
s = 6.096;                  % semi span (m)
c = 1.8288;                 % root chord (m)
c_ref = c;                  % reference chord (m)
sw = 0;                     % sweep angle at leading edge (deg)
dh = 0;                     % dihedral (deg)
tr = 1;                     % taper ratio

xf = 0.33;                  % flexural axis position relative to chord: xf=0.5 is mid chord flexural axis
dampRatio = 0.00; %0.01;    % structural damping estimate: Dff = Mff*diag(2*dampRatio*OMEGA)
xm = 0.43;                  % mass axis position relative to chord

rho_bar = 35.71;            % mass per unit length (kg/m) : int_A rho dA
I_Tym = 8.64;               % polar mass moment of inertia per unit length (kg*m) : int_A rho r^2 dA
I_zm = rho_bar*-(xm-xf)*c;  % static mass moment per unit length (kg) : int_A rho x dA  (bending-torsion coupling; negative because the mass axis lies behind the flexural axis)
EI_xxa = 9.77e6;            % flexural rigidity (N*m^2)   : E(N/m^2) * int_A z^2 dA
GI_Tya = 0.99e6;            % torsional rigidity (N*m^2)  : G(N/m^2) * int_A r^2 dA


np_s = 10;                  % number of DLM panels spanwise
np_c = 5;                   % number of DLM panels chordwise

num_ele = 10;               % number of FEM elements
num_modes = 5;              % number of structural modes to keep

[Pa,Ps] = build_PaPs(s,c,np_s,np_c,sw,dh,tr,xf);   % panel corners in aero (Pa) and structural (Ps) coordinates
cspanels = {[8,9,10]};                             % the flap: outer three panels of the trailing-edge row

% the wing in structural coordinates, seen from aft, above and towards the root: the flap, the
% accelerometer, the two chordwise axes of the benchmark (the mass axis in coral: its position
% behind the flexural axis is what makes this wing flutter) and both frame triads at their origins
x_mass = -(xm-xf)*c;                                % mass axis: (xm-xf)*c behind the flexural axis
wing_scene(struct('Pa',{Pa},'Ps',{Ps},'cspanels',{cspanels}), 'Surfaces', 1, 'IMUs', 1, ...
    'Axes', {'flexural axis $0.33c$', 0, []; 'mass axis $0.43c$', x_mass, fw_style().coral}, ...
    'Triad', 'both', 'Freestream', true, 'IMUNumbers', false, 'Name', 'goland_planform');

%% (2) The structure: a bending-torsion beam FEM
% The structure is a cantilever beam along the flexural axis, split into ten elements of equal
% length. Each node carries three degrees of freedom: the vertical deflection $u_z$, its slope
% along the span $u_z' = \partial u_z/\partial y$, and the twist angle $\psi_y$ about the beam
% axis.
%
% The element matrices are the energies written in these coordinates. Stiffness is the strain
% energy of bending and torsion,
%
% $$K_{ij} = \int_L EI\, H_i''\, H_j''\,\mathrm{d}y + \int_L GJ\, N_i'\, N_j'\,\mathrm{d}y ,$$
%
% with $H_i$ the Hermite functions of the bending degrees of freedom and $N_i$ the linear
% functions of the twist. Mass is the kinetic energy. In the general case bending and torsion
% DOFs are coupled:
%
% $$M_{ij} = \int_L \bar\rho\, H_i H_j\,\mathrm{d}y + \int_L I_{Tym}\, N_i N_j\,\mathrm{d}y - \int_L I_{zm}\,\big(N_i H_j + H_i N_j\big)\,\mathrm{d}y .$$
%
% The last term is the coupling. A point of the section at chordwise position $x$ moves
% vertically by $u_z - x\,\psi_y$. A section whose mass centre is off the flexural axis
% therefore gains kinetic energy from the product of deflection rate and twist rate, weighted by
% the static mass moment $I_{zm}$. The minus sign is that kinematic relation: positive twist
% about the $y$ axis moves a point on the positive $x$ side (toward the leading edge) in the
% negative $z$ direction, which is up. This term makes the wind-off modes a mixture of bending
% and twist, and later lets the airflow feed one motion through the other.
%
% Elements are assembled into the global matrices $\mathbf K_{gg}$ and $\mathbf M_{gg}$ (the
% subscript $g$ marks the physical FEM degrees of freedom), the three root degrees of freedom
% are clamped, and the eigenvalue problem
%
% $$\mathbf K_{gg}\,\boldsymbol\phi = \omega^2\,\mathbf M_{gg}\,\boldsymbol\phi$$
%
% gives the wind-off modes: 7.67 Hz, a bending mode, and 15.27 Hz, a torsion mode. The mode
% shapes are the columns of $\boldsymbol\Phi_{gf}$ (from physical DOFs $g$ to modal coordinates
% $f$). The modal coordinates $\mathbf q_f$, with
% $\mathbf q_g = \boldsymbol\Phi_{gf}\,\mathbf q_f$, are all the structure keeps from here on.
% In modal coordinates the mass and stiffness matrices are diagonal:
%
% $$\mathbf M_{ff} = \boldsymbol\Phi_{gf}^\top \mathbf M_{gg}\, \boldsymbol\Phi_{gf}, \qquad \mathbf K_{ff} = \mathbf M_{ff}\,\mathrm{diag}(\omega_i^2), \qquad \mathbf D_{ff} = \mathbf M_{ff}\,\mathrm{diag}(2\zeta\omega_i) ,$$
%
% where $\mathbf D_{ff}$ is a modal damping matrix built from the damping ratio $\zeta$, zero
% for this wing.
%
% **Source:** |build/build_E_y.m|, |build_K_ele.m|, |build_M_ele.m|, |build_Kgg.m|,
% |build_Mgg.m|, |build_PHIgf.m|, |beamPHI.m|.

[E, ele] = build_E_y(s, 0, num_ele);
Mgg = build_Mgg(E, ele, rho_bar, I_Tym, I_zm);
Kgg = build_Kgg(E, ele, EI_xxa, GI_Tya);

% clamp the root node (build_PHIgf drops its three DOFs), keep the first num_modes modes
[PHIgf,OMEGA] = build_PHIgf(Kgg,Mgg,num_modes);      % mode shapes (largest entry +1) and eigenfrequencies (rad/s)
disp(['First  Structural Mode Eigenfrequency in Hz:  ',num2str(OMEGA(1)/(2*pi))])
disp(['Second Structural Mode Eigenfrequency in Hz:  ',num2str(OMEGA(2)/(2*pi))])
Mff = PHIgf'*Mgg*PHIgf;
max_offdiag = max(max(abs(Mff - diag(diag(Mff)))));  % test diagonality
assert(max_offdiag < 1e-8,"Generalized Mass Matrix Mff is not diagonal")
Mff = diag(diag(Mff));
Kff = Mff*diag(OMEGA.^2);
Dff = Mff*diag(2*dampRatio*OMEGA);

% the four bending shape functions of one beam element (the formulas are in build/beamPHI.m)
plot_hermite_functions('Name', 'fem_shape_functions');

% the first two wind-off mode shapes, interpolated with the same shape functions
plot_mode_shapes(E, ele, PHIgf, OMEGA, 1:2, 'Name', 'goland_modes');

%% (3) The aerodynamics: doublet lattice
% The air loads come from the doublet-lattice method (DLM), the workhorse of industrial flutter
% analysis since Albano and Rodden published it in 1969. The DLM is a frequency domain method.
% It maps harmonically oscillating wing surfaces at a frequency $\omega$ in a freestream
% $V_\infty$, to harmonically oscillating pressure difference coefficients at the same frequency
% but with varying gain and phase per panel. The gain and phase shift is expressed as one
% complex number per possible transmission path per frequency. Each panel's flow condition
% transmits to each panel's pressure coefficient. Since the boundary condition at each panel's
% control point influences mostly its own pressure difference the resulting $n_p$ x $n_p$
% matrices are diagonally dominant. The natural frequency variable is the reduced frequency
%
% $$k = \frac{\omega\, c_\mathrm{ref}}{2\,V_\infty},$$
%
% the phase the oscillation advances while the flow crosses half a chord. $k \to 0$ is
% quasi-steady.
%
% The boundary condition is that the flow stays tangent to the moving surface. On panel $j$ it
% is imposed at the control point as a normalised downwash $w_j = u_{z,j}/V_\infty$, the
% vertical velocity of the surface there divided by the freestream. The DLM returns the pressure
% difference across each panel as a pressure coefficient. The map from downwash to pressure is
% the aerodynamic influence coefficient (AIC) matrix,
%
% $$\Delta \mathbf c_p(k) = \mathbf Q_{jj}(k)\, \mathbf w(k),$$
%
% one complex $50 \times 50$ matrix per reduced frequency (the subscript $j$ marks the panels).
% Complex, because at $k > 0$ the pressure lags the downwash: a real part in phase with the
% motion and an imaginary part in phase with its rate.
%
% **Source:** |build/build_PaPs.m| (panel corners, numbering), |build/build_Qjj.m| (geometry,
% symmetry, the calls into |DLMpro/|).

Ma = 0.0;               % mach number
SYM =1;                 % symmetry flag xz-plane for second semi-span
V_inf_ref = 100;        % reference freestream velocity

% must start at 0: the RFA of section (5) takes Qjj(0) over exactly
k_red = [0, 0.001, 0.01, 0.02, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 1.1];

max_omega_Hz = max(k_red)*2*V_inf_ref/c_ref /(2*pi);
disp(['At ',num2str(V_inf_ref),' m/s the maximal frequency covered by the RFA is:  ',num2str(max_omega_Hz),' Hz'])

Qjj = build_Qjj(Ma,k_red,c_ref,Pa,SYM);

% the panel grid in aerodynamic coordinates, numbered from the trailing-edge row inboard; panel 8
% (the inboard flap panel) carries the doublet line at 1/4 chord and the control point at 3/4 chord
plot_dlm_panel(Pa, Ps, 8, 'Name', 'goland_panels');

% two entries of the AIC matrix over the 14 reduced frequencies
plot_aic_entries(k_red, Qjj, [1 1; 1 25], 'Name', 'aic_nyquist');

%% (4) Coupling the two grids and projecting onto modes
% The structure is defined in FEM deflections, the aerodynamics in panel pressures and downwash.
% Every point of the wing surface at chordwise position $x$ (structural frame) and span station
% $y$ "rides" on a rigid bar perpendicular to the beam, so its vertical position is
%
% $$z(x, y) = u_z(y) - x\,\psi_y(y):$$
%
% the beam's deflection plus the twist times the lever arm.
%
% **Pressure to force.** A pressure coefficient $\Delta c_{p,j}$ on panel $j$ is a force
% $\bar q\, a_j\, \Delta c_{p,j}$ at its load point, with $a_j$ the panel area and
% $\bar q = \tfrac12\rho V_\infty^2$ the dynamic pressure. The force acts on the beam station
% under the load point. The shape functions evaluated there determine how it splits over the
% nodal deflections and slopes, and the lever arm to the flexural axis defines the moment that
% splits over the nodal twists. Collected for all panels,
%
% $$\mathbf F_g = \bar q\, \mathbf S_{gj}\, \Delta \mathbf c_p ,$$
%
% with $\mathbf S_{gj}$ a $33 \times 50$ matrix (physical DOFs $g$ from panels $j$).
%
% **Motion to downwash.** The downwash at the control point of panel $j$ has two parts. The
% steady part is the angle the surface presents to the flow, which for a flat plate on a
% twisting beam is the twist at that station. The unsteady part is the vertical velocity of the
% control point itself, $\dot u_z - x_j \dot\psi_y$, divided by $V_\infty$:
%
% $$\mathbf w = \mathbf D^{Re}_{jg}\, \mathbf q_g + \frac{1}{V_\infty}\, \mathbf D^{Im}_{jg}\, \dot{\mathbf q}_g .$$
%
% The names come from the frequency domain: for harmonic motion
% $\dot{\mathbf q}_g = \mathrm{i}\omega\, \mathbf q_g$, so the downwash is
% $(\mathbf D^{Re}_{jg} + \mathrm{i}\,\tfrac{\omega}{V_\infty}\mathbf D^{Im}_{jg})\,\mathbf q_g$,
% a real and an imaginary part.
%
% **Why this is consistent.** Both directions use the same six shape functions of the same
% element: the force distribution is the transpose of the displacement interpolation. The work
% done by the panel pressures on a virtual displacement of the surface then equals the work done
% by the nodal forces on the corresponding virtual nodal displacement, so no energy is created
% or lost at the interface. This is the principle of virtual work. It is why the code needs no
% surface spline (the SPLINE cards of a NASTRAN deck): the beam's own basis is the
% interpolation.
%
% **Onto the modes.** Everything is now projected onto the five mode shapes of chapter 2,
%
% $$\mathbf S_{fj} = \boldsymbol\Phi_{gf}^\top \mathbf S_{gj}, \qquad \mathbf D^{Re}_{jf} = \mathbf D^{Re}_{jg}\, \boldsymbol\Phi_{gf}, \qquad \mathbf D^{Im}_{jf} = \mathbf D^{Im}_{jg}\, \boldsymbol\Phi_{gf},$$
%
% which leaves $5 \times 50$ and $50 \times 5$ matrices. Together with $\mathbf M_{ff}$,
% $\mathbf K_{ff}$, $\mathbf D_{ff}$ and the mode shapes they fill the |Structure| dictionary;
% $\mathbf Q_{jj}(k)$, the reference chord and the air density go into |Aero|. Every |build_*|
% function downstream reads only these two structs, which is what makes the pipeline swappable:
% a bigger wing is a bigger pair of dictionaries.
%
% **Source:** |build/build_S_ele.m|, |build_Sgj.m|, |build_DReDIm_ele.m|, |build_DReDIm_jg.m|;
% the dictionary layout is documented in the header of |define/define_Goland_Structure_Aero.m|.

show_image('DLM_FEM_Coupling.png');               % drawn by docs/figures/make_coupling_diagram.py

% coupling matrices between aero (DLM) and structure (FEM)
Sgj = build_Sgj(E,ele,Ps);
[DRe_jg, DIm_jg] = build_DReDIm_jg(E,ele,Ps);

Sfj = PHIgf'*Sgj;
DRe_jf = DRe_jg*PHIgf;
DIm_jf = DIm_jg*PHIgf;

% save to the Structure & Aero dictionaries consumed by all build_* functions
Structure.OMEGA = OMEGA;
Structure.E = E;
Structure.ele = ele;
Structure.Pa = Pa;
Structure.Ps = Ps;
Structure.cspanels = cspanels;
Structure.Sgj = Sgj;
Structure.DRe_jg = DRe_jg;
Structure.DIm_jg = DIm_jg;
Structure.Mgg = Mgg;
Structure.Kgg = Kgg;

Structure.Kff = Kff;
Structure.Mff = Mff;
Structure.Dff = Dff;
Structure.Sfj = Sfj;
Structure.DRe_jf = DRe_jf;
Structure.DIm_jf = DIm_jf;
Structure.PHIgf = PHIgf;

Aero.c_ref = c_ref;
Aero.rho = rho;

%% (5) From frequency to time: Roger's RFA
% The doublet lattice gives the air loads as a table: one matrix per sampled reduced frequency,
% valid for harmonic motion only. A state-space model needs the loads as a function of the
% Laplace variable, so that any motion, growing, decaying or transient, can be simulated.
% Roger's rational function approximation (RFA) is the classic bridge. It fits every entry of
% the table with
%
% $$\hat{\mathbf Q}(k) = \mathbf Q_0 + \mathrm{i}k\, \mathbf Q_1 + \sum_{r=1}^{n_p} \mathbf A_r\, \frac{\mathrm{i}k}{\mathrm{i}k + p_r} ,$$
%
% a form with a physical reading. $\mathbf Q_0$ is the quasi-steady load, the table's entry at
% $k = 0$, taken over exactly. $\mathrm{i}k\,\mathbf Q_1$ grows with frequency: the part of the
% load in phase with the surface velocity. Each lag term is a first-order filter with a real
% pole $p_r$; together they give the load a memory of the motion a moment ago, which is the
% wake. The poles are not fitted but placed,
%
% $$p_r = \frac{k_{max}}{n_p,\ n_p - 1,\ \dots,\ 1}, \qquad n_p = 6,\ k_{max} = 1.1: \quad p_r = 0.18,\ 0.22,\ 0.28,\ 0.37,\ 0.55,\ 1.1 ,$$
%
% which keeps the fit linear in the unknowns $\mathbf Q_1$ and $\mathbf A_r$. At each sampled
% $k \ne 0$ the real and the imaginary part of one matrix entry give two real equations,
%
% $$\begin{bmatrix} 0 & \dfrac{k^2}{k^2 + p_r^2} \\ k & \dfrac{k\,p_r}{k^2 + p_r^2} \end{bmatrix} \begin{bmatrix} Q_{1} \\ A_{r} \end{bmatrix} = \begin{bmatrix} \mathrm{Re}\,(Q - Q_0) \\ \mathrm{Im}\,(Q - Q_0) \end{bmatrix} \quad (r = 1 \dots n_p \text{ side by side}),$$
%
% so the 13 nonzero samples give 26 equations for 7 unknowns per entry, solved in the
% least-squares sense.
%
% Modal projection reduces the state count. Realised in panel space, each lag term would need
% one state per panel, $50 \times 6 = 300$ aerodynamic states. Projected onto the modes first,
% as chapter 6 does, each lag term needs one state per mode: $5 \times 6 = 30$. The whole Goland
% wing then has 10 structural and 30 aerodynamic states.
%
% |evalRFA| re-evaluates $\hat{\mathbf Q}$ at the 14 samples and reports three absolute errors:
% the largest deviation of any entry at any frequency, 0.0106; the largest root-mean-square
% deviation over frequency of any entry, 0.0052; and the root-mean-square deviation over all
% entries and frequencies, 0.00088, against entries of order one.
%
% **Experiment: a cheaper fit**
%
% Take only six of the fourteen reduced frequencies and two poles instead of six, and fit again.
% Evaluated on all fourteen samples the errors grow to five to seven times the full fit.
%
% As expected a less expressive fit is less accurate; however, what is accurate enough? Chapter
% 7 answers this question by comparing the cheap aerodynamics in |Aero_c| to the p-k method's
% flutter result as baseline.
%
% **Source:** |rfa/rogersRFA_magW.m|, |util/evalRFA.m|, |util/fiterrs.m|, |util/mimo_nyquist.m|.

num_poles = 6;          % number of RFA poles
[poles,Q0jj,Q1jj,~,QLpjj,D_rog,E_rog,R_rog] = rogersRFA_magW(k_red, Qjj, num_poles);
Aero.poles = poles;
Aero.Q0jj = Q0jj;
Aero.Q1jj = Q1jj;
Aero.QLpjj = QLpjj;

% Evaluate RFA Fit Quality
[E_max_rog,E_rms_max_rog, E_rms_rog,H_rog] = evalRFA(k_red,Qjj, Q0jj,Q1jj,D_rog,E_rog,R_rog);
disp(['Maximal difference between Qjj and H_hat at any panel any frequency:  ',num2str(E_max_rog)])
disp(['Maximal rms of error over all frequencies at any panel             :  ',num2str(E_rms_max_rog)])
disp(['Rms of error over all panels and frequencies                       :  ',num2str(E_rms_rog)])
mimo_nyquist(k_red, Qjj, {H_rog}, {'DLM samples', 'Roger fit'}, [1,25], [1,25], ...
    'LegendTile', 'north', 'Name', 'rfa_nyquist')

% experiment: how much fit does the RFA really need?
idx_c = [1 2 5 7 10 14];                         % 6 of the 14 reduced frequencies: 0, 0.001, 0.05, 0.1, 0.4, 1.1
k_red_c = k_red(idx_c);  Qjj_c = Qjj(:,:,idx_c); % reuse the DLM samples: no new solve
num_poles_c = 2;
[poles_c,Q0jj_c,Q1jj_c,~,QLpjj_c,D_c,E_c,R_c] = rogersRFA_magW(k_red_c, Qjj_c, num_poles_c);
[E_max_c,E_rms_max_c,E_rms_c,H_c] = evalRFA(k_red,Qjj,Q0jj_c,Q1jj_c,D_c,E_c,R_c);   % errors on the full 14-point grid
disp(['Coarse fit, maximal difference between Qjj and H_hat at any panel any frequency:  ',num2str(E_max_c)])
disp(['Coarse fit, maximal rms of error over all frequencies at any panel             :  ',num2str(E_rms_max_c)])
disp(['Coarse fit, rms of error over all panels and frequencies                       :  ',num2str(E_rms_c)])
mimo_nyquist(k_red, Qjj, {H_rog, H_c}, {'DLM samples','6 poles, 14 $k$','2 poles, 6 $k$'}, [1,25], [1,25], ...
    'LegendTile', 'north', 'Name', 'rfa_experiment')

% keep the coarse aerodynamics for the flutter comparison of section (7)
Aero_c = Aero;  Aero_c.poles = poles_c;  Aero_c.Q0jj = Q0jj_c;  Aero_c.Q1jj = Q1jj_c;  Aero_c.QLpjj = QLpjj_c;

%% (6) Flutter analysis I: the p-method, the flutter mode
% With the rational fit in hand the aeroelastic equations become an ordinary linear system. The
% state vector stacks the modal coordinates, their rates and the aerodynamic lag states,
%
% $$\mathbf x = \begin{bmatrix} \mathbf q_f \\ \dot{\mathbf q}_f \\ \mathbf x_L \end{bmatrix}, \qquad \dot{\mathbf x} = \mathbf A(V_\infty)\,\mathbf x,$$
%
% 40 states for the bare wing (10 structural, 30 aerodynamic), and $\mathbf A$ depends on the
% airspeed through the dynamic pressure and through the lag poles. The p-method is then one line
% per velocity: the eigenvalues of $\mathbf A(V_\infty)$ are the aeroelastic modes, their real
% parts the growth rates, their imaginary parts the frequencies. |build_ABCD_G| assembles
% $\mathbf A$ (and the $\mathbf B$, $\mathbf C$, $\mathbf D$ of chapter 8, zero placeholders for
% now). |getEigenvalueModeshape| calls |eig| at each of the 42 velocities and matches the
% eigenvalues from one velocity to the next by comparing eigenvectors, so that each row of its
% output is one mode followed across the grid. |Vg_plot| draws the rows between 3 and 20 Hz, the
% band of the two modes that matter, and finds the crossing.
%
% The comparison below follows Table 2 of Murua, Palacios and Graham (2010), including its
% attribution of 175.6 m/s to Goland. These are numerical benchmarks with different aerodynamic
% formulations and discretizations, rather than measurements of one physical wing.
%
% - Goland (1945), as reported by Murua et al., analytical: 175.6 m/s
% - Wang et al. (2006), ZAERO (panel method): 174.3 m/s
% - Wang et al. (2006), UVLM: 163.8 m/s
% - Murua et al. (2010), UVLM: 165 m/s
% - Murua et al. (2010), UVLM with RFA: 177 m/s
% - this model, DLM, Roger RFA, p-method: 164.8 m/s
%
% References: J. Murua, R. Palacios and J. M. R. Graham, [*Modeling of Nonlinear Flexible
% Aircraft Dynamics Including Free-Wake
% Effects*](https://openresearch.surrey.ac.uk/esploro/outputs/conferencePresentation/Modeling-of-Nonlinear-Flexible-Aircraft-Dynamics/99511820602346),
% AIAA 2010-8226, Table 2; Z. Wang, P. C. Chen, D. D. Liu, D. T. Mook and M. J. Patil, [*Time
% Domain Nonlinear Aeroelastic Analysis for HALE Wings*](https://doi.org/10.2514/6.2006-1640),
% AIAA 2006-1640. The original benchmark is M. Goland, *The Flutter of a Uniform Cantilever
% Wing*, Journal of Applied Mechanics 12(4), A197-A208 (1945).
%
% This model is close to the two UVLM results. Chapter 7 provides a separate check on its
% rational fit: the p-k method, which uses the doublet-lattice table directly, gives a flutter
% speed within 0.3 m/s. That comparison tests the fit on this grid; a convergence study would
% also vary the beam mesh, retained modes and panel grid.
%
% The flutter mode itself is the eigenvector of the unstable pole. At 180 m/s the pole sits at
% 10.9 Hz with a growth of 42 % per cycle. Its eigenvector describes the amplitudes and phase
% relation of bending and twist. Their coupled motion lets the aerodynamic forces do positive
% net work over a cycle. Modal signs depend on the eigenvector convention; interpret the phase
% through the reconstructed physical motion.
%
% **Experiment: move the mass axis**
%
% The mass coupling of chapter 2 ties bending to twist inside the structure. Move the mass axis
% from 0.43 to 0.30 of the chord, ahead of the flexural axis at 0.33, and $I_{zm}$ changes sign.
% Only the mass matrix changes, so the rebuild is cheap: new $\mathbf M_{gg}$, new modes, new
% projection, the aerodynamics untouched.
%
% With the mass axis at 0.30 c there is no flutter: the damping of both modes rises with
% velocity instead of turning down. This is mass balancing, the classical passive cure for
% flutter. With the mass centre ahead of the flexural axis, the inertia of an accelerating
% section twists it the other way, and the airflow can no longer feed the twist through the
% bending. Try the other knobs while the script is open: sweep |xm|, give the wing structural
% damping through |dampRatio|, or stiffen it through |GI_Tya|.
%
% **Source:** |build/build_ABCD_G.m|, |util/getEigenvalueModeshape.m|, |util/eigenshuffle/|,
% |util/Vg_plot.m|, |util/animate_wing.m|.

V_inf = linspace(20,190,42);                     % the velocity grid of build_G_Goland and of the README

num_panels=size(Q0jj,1);
num_dof=size(Kgg,1);
num_AIL = 1;
num_IMU = 1;
Structure.DRe_jx = zeros(num_panels,num_AIL);   % placeholder - built in section (8)
Structure.DIm_jx = zeros(num_panels,num_AIL);   % placeholder - built in section (8)
Structure.PHIzg = zeros(num_IMU,num_dof);       % placeholder - built in section (8)

ny = num_IMU;
nu = 3*num_AIL;     % inputs: flap deflection, rate and acceleration
nv = length(V_inf);
G = ss(zeros(ny,nu,nv));
for i_v = 1:length(V_inf)
    V_inf_iv = V_inf(i_v);
    [A_noS,B_noS,C_noS,D_noS] = build_ABCD_G(V_inf_iv,Structure,Aero);
    G_noAct = ss(A_noS,B_noS,C_noS,D_noS);
    G(:,:,i_v) = G_noAct;
end

[EV_G_p, MS_G_p] = getEigenvalueModeshape(G,num_modes);
[~, cr_p] = Vg_plot(V_inf, EV_G_p, 'open loop', 'Band', [3 20], 'Name', 'vg_p');

% flutter speeds reported in the literature for the Goland wing:
%   Table 2 of Murua et al. (2010), cited in docs/TUTORIAL.md chapter 6:
%   Goland (reported there) 175.6 m/s, Wang (ZAERO) 174.3 m/s, Wang (UVLM) 163.8 m/s,
%   Murua (UVLM) 165 m/s, Murua (UVLM + RFA) 177 m/s

V_anim = 180;                                    % 15 m/s past the onset
A_anim = build_ABCD_G(V_anim,Structure,Aero);
[v_anim,e_anim] = eig(A_anim);  e_anim = diag(e_anim);
cand = find(real(e_anim) > 0 & imag(e_anim) > 3);   % unstable oscillatory poles
[~,imax] = max(real(e_anim(cand)));  lambda_f = e_anim(cand(imax));
q_f_vis = v_anim(1:num_modes, cand(imax));       % complex modal eigenvector: bending and torsion with their phase
fprintf('Flutter mode at %g m/s: %.2f Hz, amplitude x%.3f per cycle\n', V_anim, imag(lambda_f)/(2*pi), exp(real(lambda_f)*2*pi/imag(lambda_f)))
fprintf('Torsion/bending amplitude ratio %.3f, phase %+.1f deg\n', abs(q_f_vis(2))/abs(q_f_vis(1)), angle(q_f_vis(2)/q_f_vis(1))*180/pi)
animate_wing(Structure, q_f_vis, 'Cycles', [3 1], 'FramesPerCycle', 14, 'TipAmplitude', 0.12, 'Skin', 'flat');
% The fluttermodes animation Replay is not available inside the live script. Please run above
% command in the Command Window to open the interactive Figure.

% experiment: move the mass axis ahead of the flexural axis (mass balancing).
% Only the structure changes - the aerodynamics (DLM + RFA) are untouched.
xm_exp = 0.30;                          % mass axis ahead of the flexural axis
I_zm_exp = rho_bar*-(xm_exp-xf)*c;
Mgg_exp = build_Mgg(E, ele, rho_bar, I_Tym, I_zm_exp);
[PHIgf_exp,OMEGA_exp] = build_PHIgf(Kgg,Mgg_exp,num_modes);

Structure_exp = Structure;
Structure_exp.Mgg = Mgg_exp;
Structure_exp.PHIgf = PHIgf_exp;
Structure_exp.OMEGA = OMEGA_exp;
Mff_exp = PHIgf_exp'*Mgg_exp*PHIgf_exp;
Mff_exp = diag(diag(Mff_exp));
Structure_exp.Mff = Mff_exp;
Structure_exp.Kff = Mff_exp*diag(OMEGA_exp.^2);
Structure_exp.Dff = Mff_exp*diag(2*dampRatio*OMEGA_exp);
Structure_exp.Sfj = PHIgf_exp'*Sgj;
Structure_exp.DRe_jf = DRe_jg*PHIgf_exp;
Structure_exp.DIm_jf = DIm_jg*PHIgf_exp;

G_exp = ss(zeros(ny,nu,nv));
for i_v = 1:length(V_inf)
    [A_e,B_e,C_e,D_e] = build_ABCD_G(V_inf(i_v),Structure_exp,Aero);
    G_exp(:,:,i_v) = ss(A_e,B_e,C_e,D_e);
end
[EV_exp, ~] = getEigenvalueModeshape(G_exp,num_modes);
[~, cr_exp] = Vg_plot({V_inf, V_inf}, {EV_G_p, EV_exp}, {'mass axis at $0.43c$', 'mass axis at $0.30c$'}, ...
    'Band', [3 20], 'Name', 'vg_mass_axis');

%% (7) Flutter analysis II: the p-k method as reference
% The p-method is fast because it does not run the doublet lattice method after the initial fit.
% The classical way, avoiding the time-domain altogether, is the p-k method. Start from the
% modal equation of motion with the aerodynamic force on the right,
%
% $$\mathbf M_{ff}\,\ddot{\mathbf q}_f + \mathbf D_{ff}\,\dot{\mathbf q}_f + \mathbf K_{ff}\,\mathbf q_f = \bar q\,\mathbf Q_{ff}(k)\,\mathbf q_f, \qquad \mathbf Q_{ff}(k) = \mathbf S_{fj}\,\mathbf Q_{jj}(k)\,\Big(\mathbf D^{Re}_{jf} + \mathrm{i}\,\tfrac{2k}{c_\mathrm{ref}}\,\mathbf D^{Im}_{jf}\Big),$$
%
% the AIC table sandwiched between the coupling matrices. Its real part acts as a stiffness and
% its imaginary part, for harmonic motion at $\omega$, as a damping:
% $\mathrm{i}\,\bar q\,\mathbf Q_{Im}\,\mathbf q_f = \bar q\,\mathbf Q_{Im}\,\dot{\mathbf q}_f/\omega$.
% Seeking $\mathbf q_f \propto e^{pt}$ gives the aeroelastic eigenproblem
%
% $$\Big[\,\mathbf M_{ff}\,p^2 + \big(\mathbf D_{ff} - \bar q\,\mathbf Q_{Im}/\omega\big)\,p + \big(\mathbf K_{ff} - \bar q\,\mathbf Q_{Re}\big)\Big]\,\mathbf x = \mathbf 0 .$$
%
% The "argument" is however circular: $\mathbf Q_{ff}$ depends on the reduced frequency,
% $k = \omega\,c_\mathrm{ref}/(2V_\infty)$ depends on the frequency of the mode, and that
% frequency is the answer. The solution is iterative. Guess $k$, evaluate the aerodynamics there
% (one doublet-lattice solve), solve the eigenproblem, read the mode's $p$, set
% $k_{n+1} = |\mathrm{Im}\,p_n|\,c_\mathrm{ref}/(2V_\infty)$, and repeat until $k$ stops moving.
% Because the eigenvalues come out of |eig| in no particular order, the mode is followed through
% the iterations and across the velocities with |eigenshuffle|, the same matching the p-method
% uses.
%
% The two flutter methods must agree: the p-k method puts the crossing at 164.5 m/s and 11.34
% Hz, the p-method at 164.8 m/s and 11.30 Hz, 0.3 m/s apart. That agreement is the first thing
% to verify on any new wing: a fit that has gone wrong shows up as a gap between the two curves.
%
% **Experiment: the cheap fit meets the referee**
%
% Chapter 5 kept the two-pole, six-sample fit in |Aero_c|. Run the p-method with it (the p-k
% result does not depend on any fit, so it needs no rerun) and compare all three.
%
% The coarse aerodynamics move the flutter speed from 164.8 to 162.8 m/s, a change of 2 m/s or
% 1.2 %, and the flutter frequency from 11.30 to 11.13 Hz, while the six-pole fit reproduces the
% p-k reference to 0.3 m/s. For this wing the cheap fit would still be a usable estimate.
%
% **Source:** |util/pkmethode.m| (Ma = 0 and SYM = 1 fixed inside; must match the define file),
% |util/eigenshuffle/|.

V_inf_pk = linspace(90,190,11);                  % coarser grid: every point costs one DLM solve per iterate
[EV_G_pk, MS_G_pk] = pkmethode(Structure,Aero, V_inf_pk);   % about 10 s: one DLM solve per iterate
[~, cr_pk] = Vg_plot({V_inf, V_inf_pk}, {EV_G_p, EV_G_pk}, ...
    {'p method', 'p-k method'}, 'Band', [3 20], 'Name', 'vg_p_vs_pk');

% experiment: the same p-method sweep on the coarse RFA of section (5)
G_c = ss(zeros(ny,nu,nv));
for i_v = 1:length(V_inf)
    [A_c,B_c,C_c,D_c2] = build_ABCD_G(V_inf(i_v),Structure,Aero_c);
    G_c(:,:,i_v) = ss(A_c,B_c,C_c,D_c2);
end
[EV_c, ~] = getEigenvalueModeshape(G_c,num_modes);
[~, cr_c] = Vg_plot({V_inf, V_inf, V_inf_pk}, {EV_G_p, EV_c, EV_G_pk}, ...
    {'p method, 6 poles / 14 $k$', 'p method, 2 poles / 6 $k$', 'p-k reference'}, 'Band', [3 20], ...
    'Name', 'vg_rfa_coarse');

%% (8) Control surface, accelerometer, actuator
% **Flap.** A control surface is a list of panels plus a hinge line, given by two points in
% aerodynamic coordinates. For the benchmark these are the outer three panels of the
% trailing-edge row and the line from the inboard leading-edge corner of panel 8
% ($\mathbf{RP}_1$) to the outboard one of panel 10 ($\mathbf{RP}_2$), the grey surface of the
% chapter 1 figure. Deflecting it by an angle $u_x$ tilts those panels and moves their control
% points, so the flap enters the aerodynamics as a localized downwash boundary condition:
%
% $$\mathbf w = \mathbf D^{Re}_{jx}\,u_x + \frac{1}{V_\infty}\,\mathbf D^{Im}_{jx}\,\dot u_x, \qquad D^{Re}_{jx} = \mathbf r\cdot\mathbf S_j, \quad D^{Im}_{jx} = -\mathbf N_j\cdot\big(\mathbf r\times(\mathbf j - \mathbf{RP}_1)\big),$$
%
% with $\mathbf r$ the unit hinge vector, $\mathbf S_j$ the panel's spanwise direction,
% $\mathbf N_j$ its normal and $\mathbf j$ its control point. The steady term is the tilt of the
% surface, the rate term the vertical velocity of the control point swinging about the hinge.
% The sign convention follows from the hinge direction: a positive deflection is the trailing
% edge down. The columns of $\mathbf D^{Re}_{jx}$ and $\mathbf D^{Im}_{jx}$ are then multiplied
% by the same aerodynamics as the structural downwash, which is how a flap command reaches the
% modal equations.
%
% **Accelerometer.** The sensor is an inertial measurement unit, here a single-axis
% accelerometer, that reads the vertical acceleration of the point it sits on. By the rigid-bar
% kinematics of chapter 4,
%
% $$\ddot u_{z,IMU} = \ddot u_z(y) - x_{IMU}\,\ddot\psi_y(y),$$
%
% which |build_IMU| writes as a row $\boldsymbol\Phi_{zg}$ over the physical degrees of freedom
% (subscript $z$ for sensors), evaluated through the shape functions at the sensor's span
% station and multiplied by the lever arm at its chord position. The reading is along $z_s$, so
% a wing that accelerates upward reads a negative value. The sensor sits at the hinge midpoint,
% 0.8 chord from the leading edge at 85 % of the span, where both the bending and the twist of
% the flutter mode are visible.
%
% **Actuator.** A commanded deflection is not an instant deflection. Control surfaces are driven
% by servos whose bandwidth is of the order of the frequencies of flutter, so the actuator model
% is not optional! The simplest one that is not too simple is second order,
%
% $$G_{act}(s) = \frac{K\,\omega_0^2}{s^2 + 2\,d\,\omega_0\,s + \omega_0^2}, \qquad K = 1,\quad d = 0.9,\quad \omega_0 = 2\pi\cdot 32\ \mathrm{rad/s},$$
%
% from the demanded angle to the actual one, well damped ($d = 0.9$) with a 32 Hz bandwidth. It
% is realised as a two-state model whose outputs are the deflection, its rate and its
% acceleration, because the aerodynamics needs all three: the deflection drives
% $\mathbf D^{Re}_{jx}$, the rate $\mathbf D^{Im}_{jx}$, and the acceleration the apparent-mass
% term. Its rate state is scaled by $\omega_0$ (|StateScale = [1; w0]|) so that both states stay
% of comparable size.
%
% A typical actuator model must include deflection and rate limits. Both are missing from this
% linear model for now.
%
% **Source:** |build/build_DReDIm_jx.m|, |build/build_IMU.m|, |build/beamPHI.m|,
% |define/define_PT2Actuator.m|; the hinge and sensor definitions at the end of
% |define/define_Goland_Structure_Aero.m|.

rot_axis_a{1}{1} = Pa{cspanels{1}(1)}{1};
rot_axis_a{1}{2} = Pa{cspanels{1}(end)}{4};
[DRe_jx,DIm_jx] = build_DReDIm_jx(cspanels,rot_axis_a,Pa);

Structure.DRe_jx = DRe_jx;
Structure.DIm_jx = DIm_jx;

rot_axis_flap_1_s = Ps{cspanels{1}(1)}{1};
rot_axis_flap_2_s = Ps{cspanels{1}(end)}{4};

HP_flap = 0.5*(rot_axis_flap_1_s+rot_axis_flap_2_s);
x_pos_IMU=HP_flap(1);
y_pos_IMU=HP_flap(2);
[PHIzg, ~] = build_IMU(x_pos_IMU, y_pos_IMU, E, ele);

Structure.PHIzg = PHIzg;

% PT2 actuator G_act(s) = K*w0^2/(s^2 + 2*d*w0*s + w0^2), K = 1, d = 0.9 (define/define_PT2Actuator.m)
G_act      = define_PT2Actuator(num_AIL, 32);    % 32 Hz: the servo the controller was tuned with
G_act_slow = define_PT2Actuator(1, 10);          % 10 Hz: Bode comparison and experiment of section (10)
G_act_slow = G_act_slow{1};

% both from demanded to actual deflection, against the two wind-off frequencies
plot_actuator_bode({G_act{1}, G_act_slow}, {'32 Hz','10 Hz'}, OMEGA(1:2)/(2*pi), 'Name', 'actuator_bode');

%% (9) The plant G and why AFS is hard
% The previous chapters have treated the coupling of the aerodynamics to the structure, the
% vertical acceleration measurements, the integration of control surface deflections as
% localized boundary conditions on the panel downwash, and the simplest suitable actuator model.
% It is now time to stitch it all together. The resulting aeroservoelastic plant $G$ takes flap
% command inputs and outputs accelerometer readings. The plant is assembled at
% |build_G_Goland(V_inf, 1, 1)|.
%
% **Scaling.** The raw matrices of |build_ABCD_G| mix metres, radians and generalised forces;
% for numerical reasons a diagonal state and output scaling is introduced.
%
% $$\tilde{\mathbf A} = \mathbf S_x^{-1}\mathbf A\,\mathbf S_x, \qquad \tilde{\mathbf B} = \mathbf S_x^{-1}\mathbf B, \qquad \tilde{\mathbf C} = \mathbf S_o^{-1}\mathbf C\,\mathbf S_x, \qquad \tilde{\mathbf D} = \mathbf S_o^{-1}\mathbf D,$$
%
% with
% |StateScale = [1 (five modal coordinates); OMEGA (five rates); q_bar_ref*c_ref^2 (thirty lag states)]|
% and |OutScale = Kff(1,1)/Mff(1,1)|, which is $\omega_1^2$: the accelerometer output is
% measured in units of mode 1's acceleration at unit amplitude. The state scaling is a
% similarity transform for numerical reasons only. The output scaling does change the units of
% $y$, and the shipped controllers were tuned on the scaled measurement, so keep it. The
% |V_inf_ref = 100| inside the builders enters only through $\bar q_\mathrm{ref}$ in the
% lag-state scale; it is not the velocity of the model.
%
% **Names.** The scaled model gets |StateName| (|q_f1| to |q_f5|, |q_f1_dot| to |q_f5_dot|,
% |aero_lag1| to |aero_lag30|), the three actuator-side inputs |flap1|, |flap1_dot|,
% |flap1_ddot| and the output |u_z1_ddot|.
%
% **Actuator in.** |connect| wires the actuator's three outputs to those three inputs by name
% and returns the model from |flap1_d|, the demanded deflection, to |u_z1_ddot|. The result has
% 42 states, one input and one output, one |ss| per velocity, stacked into an |ss| array over
% the grid.
%
% **Why active flutter suppression is hard**
%
% The pole map shows how the dynamics change with velocity:
%
% - The dynamics depend strongly on the freestream velocity, and in a real aircraft on the
%   flight condition and the loading (fuel, payload) as well.
% - This Goland model has 42 states. Industrial models have several hundreds.
% - The 30 aerodynamic lag states have no direct physical interpretation and cannot be measured.
% - The supplied sensor models measure accelerations, so modal deflections and rates are not
%   available as direct feedback signals.
% - Each control surface excites several aeroelastic modes at once with varying gain and lag per
%   frequency and freestream velocity.
% - A feedback from distributed accelerations to control-surface deflections has no intuitive
%   connection to the modes it is supposed to damp.
% - Controller design also has to account for measurement noise, unintended interference with
%   the rigid body dynamics (primary flight controller) and avoid excessive control activity.
%
% **Source:** |build/build_G_Goland.m| (the scaling, the names, the |connect| loop),
% |build/build_ABCD_G.m|; |util/pole_plot.m|.

num_x_L         = num_modes*num_poles;      % number of aerodynamic lag states
q_bar_ref           = 0.5*rho*V_inf_ref^2;
StateScale = [ones(num_modes,1);
              OMEGA;
              repmat(q_bar_ref*c_ref^2,num_x_L,1)];

u_z_ddot_scale   = repmat(Kff(1,1)/Mff(1,1), num_IMU,1);
OutScale = u_z_ddot_scale;
InputNameNoAct = {'flap1', 'flap1_dot', 'flap1_ddot'};
OutputNameNoAct = 'u_z1_ddot';

StateNameNoAct = cell(1,2*num_modes+num_x_L);
for i = 1:num_modes
    StateNameNoAct{i} = ['q_f',num2str(i)];
    StateNameNoAct{num_modes+i} = ['q_f',num2str(i),'_dot'];
end
for i = 1:num_x_L
    StateNameNoAct{2*num_modes+i} = ['aero_lag',num2str(i)];
end

ny = num_IMU;
nu = num_AIL;
nv = length(V_inf);
G = ss(zeros(ny,nu,nv));
for i_v = 1:length(V_inf)
    V_inf_iv = V_inf(i_v);
    [A_noS,B_noS,C_noS,D_noS] = build_ABCD_G(V_inf_iv,Structure,Aero);
    G_noAct = ss(diag(1./StateScale)*A_noS*diag(StateScale), ...
                 diag(1./StateScale)*B_noS, ...
                 diag(1./OutScale)*C_noS*diag(StateScale), ...
                 diag(1./OutScale)*D_noS);

    G_noAct.InputName = InputNameNoAct;
    G_noAct.StateName = StateNameNoAct;
    G_noAct.OutputName = OutputNameNoAct;
    G_iv = connect(G_noAct,G_act{1:end},{'flap1_d'},{'u_z1_ddot'});
    G(:,:,i_v) = G_iv;
end
pole_plot(G, 'open loop', 'Velocity', V_inf, 'Limits', [-80 20 -130 130], ...
    'Title', 'Aeroelastic poles over velocity', 'Name', 'pzmap_G');

[EV_G, ~] = getEigenvalueModeshape(G,num_modes);      % the actuated plant, tracked over velocity

%% (10) Closing the loop
% Let's now do control. Feedback control to be precise.
%
% |CL = feedback(G, -Goland_Cont_nO);|
%
% Closes the loop from acceleration measurements to control surface deflection commands using
% the provided linear controller. Note the minus sign. MATLAB's |feedback(G, K)| uses negative
% feedback, $u = -K\,y$. The supplied controllers use $u = K\,y$, the convention of |lft(P, K)|.
%
% |CL| is an |ss| array over the same velocities as |G|, so |getEigenvalueModeshape| and
% |Vg_plot| apply as before. Every sampled model is stable through 190 m/s. Above 160 m/s the
% former flutter mode has at least 20 % damping; the bending mode gains damping and its
% frequency moves toward 6 Hz. The pole map shows the result at 190 m/s.
%
% The synthesis script |R00_Synthesize_Goland_Controller.m| ends with the same test in code: it
% rebuilds the plant on |linspace(20,190,42)|, closes the loop and asserts that every
% closed-loop eigenvalue at every velocity has a negative real part.
%
% **Experiment: a slower actuator and other sensor positions**
%
% The loop above is the one the controller was designed for: a 32 Hz servo and the accelerometer
% at the flap hinge. The actuator and the sensor position are two things a control engineer does
% not get to choose freely, so keep the controller and change those instead. First the servo
% bandwidth from 32 to 10 Hz, the same second-order model with $\omega_0 = 2\pi\cdot 10$: at the
% flutter frequency it now delivers half the demanded angle almost a quarter cycle late (chapter
% 8). Then, back at 32 Hz, the accelerometer moved 0.3 chord toward the leading edge, and,
% separately, to mid span at the hinge's chord position. Each variant is a new plant with the
% same controller.
%
% With the 10 Hz servo the closed loop becomes unstable at 165.2 m/s and 12.1 Hz, within half a
% meter per second of the open-loop flutter speed. The slower actuator reduces the flap response
% and adds phase lag, flutter suppression becomes infeasible.
%
% The sensor relocations preserve stability on the sampled grid but reduce damping. Moving the
% accelerometer 0.3 chord forward reduces its lever arm behind the flexural axis from 0.47 c to
% 0.17 c, so it measures less twist. Flutter-mode damping falls to about 5 % at 190 m/s. At mid
% span both mode shapes are smaller, and the damping falls to about 2 %. These cases show why
% actuator dynamics and exact sensor position are crucial when tuning an AFS controller.
%
% **Source:** |R00_Synthesize_Goland_Controller.m| (the stability assert),
% |util/getEigenvalueModeshape.m|, |util/Vg_plot.m| (|MaxDamping|), |build/build_IMU.m|.

load('Goland_Cont_nO4_lpf.mat')                  % Goland_Cont_nO: order 4 plus low-pass, from R00_Synthesize_Goland_Controller.m
CL = feedback(G,-Goland_Cont_nO);                % u = K*y, the sign convention the controller was tuned with (same as lft(P,K))
pole_plot({G(:,:,end), CL(:,:,end)}, {'open loop', 'closed loop'}, 'Limits', [-80 20 -130 130], ...
    'LegendLocation', 'southwest', 'Name', 'pzmap_ol_cl', ...
    'Title', sprintf('Aeroelastic poles at $V_\\infty = %g$ m/s', V_inf(end)));
[EV_CL, ~] = getEigenvalueModeshape(CL,num_modes);
% MaxDamping 15 hides the controller's own pole pair at 8.4 Hz with 24 % damping
[~, cr_cl] = Vg_plot({V_inf, V_inf}, {EV_G, EV_CL}, {'open loop', 'closed loop'}, ...
    'Band', [3 20], 'MaxDamping', 15, 'Name', 'vg_ol_vs_cl');

% experiment: the same controller on a slower actuator and on two other sensor positions
[PHIzg_LE,~]  = build_IMU(x_pos_IMU + 0.3*c, y_pos_IMU, E, ele);   % IMU moved 0.3 c toward the leading edge (x_s points upstream)
[PHIzg_mid,~] = build_IMU(x_pos_IMU, 0.5*s, E, ele);               % IMU at mid span, same chord position

variants = struct('label', {'32 Hz actuator, IMU at the hinge', '10 Hz actuator', ...
                            'IMU $0.3c$ toward the leading edge', 'IMU at mid span'}, ...
                  'G_act', {G_act{1}, G_act_slow, G_act{1}, G_act{1}}, ...
                  'PHIzg', {PHIzg, PHIzg, PHIzg_LE, PHIzg_mid});

EV_var = cell(1,numel(variants));
for i = 1:numel(variants)
    Structure_var = Structure;  Structure_var.PHIzg = variants(i).PHIzg;
    G_var = ss(zeros(ny,nu,nv));
    for i_v = 1:nv
        [A_v,B_v,C_v,D_v] = build_ABCD_G(V_inf(i_v),Structure_var,Aero);
        G_v = ss(diag(1./StateScale)*A_v*diag(StateScale), diag(1./StateScale)*B_v, ...
                 diag(1./OutScale)*C_v*diag(StateScale), diag(1./OutScale)*D_v);
        G_v.InputName = InputNameNoAct;  G_v.StateName = StateNameNoAct;  G_v.OutputName = OutputNameNoAct;
        G_var(:,:,i_v) = connect(G_v, variants(i).G_act, {'flap1_d'}, {'u_z1_ddot'});
    end
    EV_var{i} = getEigenvalueModeshape(feedback(G_var,-Goland_Cont_nO),num_modes);
    re_max = max(real(EV_var{i}(:)));
    verdict = {'UNSTABLE somewhere on the grid','stable over the whole grid'};
    fprintf('%-38s max Re = %+8.3f rad/s   %s\n', variants(i).label, re_max, verdict{1 + (re_max < 0)})
end
[~, cr_act] = Vg_plot({V_inf,V_inf}, EV_var(1:2), {variants(1:2).label}, ...
    'Band', [3 20], 'MaxDamping', 15, 'Name', 'vg_experiment_actuator');
[~, cr_imu] = Vg_plot({V_inf,V_inf,V_inf}, EV_var([1 3 4]), {variants([1 3 4]).label}, ...
    'Band', [3 20], 'MaxDamping', 15, 'Name', 'vg_experiment_imu');

%% (11) Virtual flight: AFS in the time domain
% Finally let's run some time simulations and visualize the controller in action.
%
% The code uses modal multisine disturbances to illustrate turbulence-like excitation; neither
% includes a physical gust model. The time simulation complements the fixed-velocity eigenvalue
% analysis. Velocity now varies with time, so the model is linear parameter-varying,
%
% $$\dot{\mathbf x} = \mathbf A\big(V_\infty(t)\big)\,\mathbf x + \mathbf B\big(V_\infty(t)\big)\,\mathbf u ,$$
%
% which the Control System Toolbox represents as an |lpvss| object: a data function that returns
% the state-space matrices at any value of the scheduling parameter.
% |build_LPV_P_Goland(100, 1, 1, 1:5)| returns the generalized plant of chapter 12 in that form.
% |lsim| integrates it along a velocity trajectory, and |lft(LPV_P, K)| closes the loop.
%
% Our first virtual flight is 60 s. The velocity ramps steadily from 150 to 185 m/s while a
% continuous multisine disturbance, twelve sinusoids between 3 and 20 Hz on each of the first
% two modes, shakes the wing like light turbulence. Below the flutter boundary the wing moves at
% the centimeter level. Beyond it the flutter mode grows exponentially. Once the tip deflection
% exceeds 2.5 % of the semispan (15 cm), the controller of chapter 10 is switched on.
%
% Damping this oscillation after the tip reaches 15 cm requires a peak flap command of about 26
% degrees. The simulation imposes no actuator deflection or rate limits, so it does not
% establish whether a particular servo could actually deliver that response.
%
% The |ss| array supports analysis and synthesis at sampled velocities; the |lpvss| supports
% simulation along a trajectory. Both use |build_ABCD_P| with the same state ordering and
% scaling, so the supplied controller connects to either form.
%
% **Source:** |build/build_LPV_P_Goland.m|, |build/build_ABCD_P.m|;
% |docs/figures/make_virtual_flight_gif.m| for the animation.

% LPV generalized plant with all 5 modes as performance outputs
% (V_inf_ref = 100 matches the controller synthesis; runs its own DLM solve)
LPV_P = build_LPV_P_Goland(100,1,1,1:5);

% vertical tip deflection at leading and trailing edge from the modal states
% (structural x points upstream: LE = max x, TE = min x)
xLE_G = max(cellfun(@(p) p{1}(1), Ps));
xTE_G = min(cellfun(@(p) p{2}(1), Ps));
PHIz_tip = build_IMU([xLE_G;xTE_G],[s;s],E,ele);
Ctip = PHIz_tip*PHIgf;
q_f_scale = 1./(sqrt(diag(Kff))/sqrt(Kff(1,1)));    % modal outputs of LPV_P are scaled

% the virtual flight: velocity ramp, continuous multisine disturbance
T = 60; t_sim = (0:0.002:T)'; Nt = length(t_sim);
V0 = 150; V1 = 185; V_t = V0 + (V1-V0)*t_sim/T;

rng(42)                     % reproducible disturbance
n_sin = 12; sigma_d = 0.08; % amplitude chosen for ~1 cm buffet below flutter
f_sin = 3 + 17*rand(n_sin,2); ph_sin = 2*pi*rand(n_sin,2);
uOL = zeros(Nt,7);          % inputs: [dist_q_f1..5_ddot, noise_u_z1_ddot, flap1_d]
uOL(:,1) = sigma_d/sqrt(n_sin)*sum(sin(2*pi*t_sim*f_sin(:,1)' + ph_sin(:,1)'),2);
uOL(:,2) = sigma_d/sqrt(n_sin)*sum(sin(2*pi*t_sim*f_sin(:,2)' + ph_sin(:,2)'),2);

% open-loop flight: one LTV simulation along the velocity trajectory
[yOL,~,xOL] = lsim(LPV_P, uOL, t_sim, [], V_t);
z_OL = (yOL(:,1:num_modes).*q_f_scale') * Ctip';    % physical tip deflection [LE, TE]

% controller ON once the deflection gets significant in the flutter regime
z_thresh = 0.025*s;         % 2.5% of semi span tip deflection (~0.15 m)
i_on = find(abs(z_OL(:,1)) > z_thresh, 1);

% closed-loop flight, restarted from the switch-on state: u = K*y via lft
CL = lft(LPV_P, Goland_Cont_nO);                    % lpvss, states = [plant; controller]
x_on = [xOL(i_on,:)'; zeros(order(Goland_Cont_nO),1)];
[yCL,~] = lsim(CL, uOL(i_on:end,1:6), t_sim(i_on:end), x_on, V_t(i_on:end));
z_CL = (yCL(:,1:num_modes).*q_f_scale') * Ctip';

zrec = [z_OL(1:i_on-1,:); z_CL];                    % stitch OL + CL histories
urec = [zeros(i_on-1,1);  yCL(:,2*num_modes+1)];    % flap1_d demand (rad)
disp(['Controller switched on at t = ',num2str(t_sim(i_on),3),' s  (V_inf = ',num2str(V_t(i_on),4),' m/s)'])
disp(['Peak tip deflection: ',num2str(max(abs(zrec(:,1))),3),' m,  peak flap command: ',num2str(max(abs(urec))*180/pi,3),' deg'])

% tip deflection and flap command in the 1.5 s around the switch-on
plot_virtual_flight(t_sim, zrec, urec, t_sim(i_on), V_t(i_on), z_thresh, 'Name', 'virtual_flight');

%% (12) How the shipped controller was designed (optional)
% Controller synthesis in this repository is fixed-structure H-Infinity optimization using
% |systune|
%
% **The channels.** |build_P_Goland| adds to $G$ a disturbance input on the acceleration of each
% selected mode, a noise input on the sensor, and as outputs the selected modal coordinates,
% their rates and the demand itself. For modes 1 and 2:
%
% $$\begin{bmatrix} \mathbf w \\ u \end{bmatrix} = \begin{bmatrix} \texttt{dist\_q\_f1\_ddot} \\ \texttt{dist\_q\_f2\_ddot} \\ \texttt{noise\_u\_z1\_ddot} \\ \texttt{flap1\_d} \end{bmatrix}, \qquad \begin{bmatrix} \mathbf z \\ y \end{bmatrix} = \begin{bmatrix} \texttt{q\_f1} \\ \texttt{q\_f2} \\ \texttt{q\_f1\_dot} \\ \texttt{q\_f2\_dot} \\ \texttt{flap1\_d} \\ \texttt{u\_z1\_ddot} \end{bmatrix} .$$
%
% The partition is by size, the convention |lft(P, K)| expects: the last input is $u$, the last
% output is $y$, and everything before them is $w$ and $z$.
%
% **The scaling.** The modal channels are scaled so that one unit of every structural channel
% carries the same energy. With $K_{ff,11}$ the modal stiffness of mode 1,
%
% - modal coordinate $q_{f,i}$: $\sqrt{K_{ff,ii}/K_{ff,11}}\;q_{f,i}$, one unit is the strain
%   energy $\tfrac12 K_{ff,11}$
% - modal rate $\dot q_{f,i}$: $\sqrt{M_{ff,ii}/K_{ff,11}}\;\dot q_{f,i}$, one unit is the
%   kinetic energy $\tfrac12 K_{ff,11}$
% - modal disturbance: $\sqrt{M_{ff,ii}/K_{ff,11}}\;w_i$, one unit is the same power into every
%   mode
% - measured acceleration, noise: physical value divided by $\omega_1^2 = K_{ff,11}/M_{ff,11}$,
%   one unit is mode 1's acceleration at unit amplitude
% - flap demand: unscaled, one unit is one radian
%
% Without it the modal stiffnesses, which span a factor 27 between mode 1 and mode 5, would make
% the optimiser weigh the modes by a numerical accident of the mode-shape normalisation. The
% chosen scaling norms ensure |systune| minimises physical energy.
%
% **The design.** |R00_Synthesize_Goland_Controller.m| tunes a fourth-order state-space block
% with a fixed first-order low-pass, on $P$ at six velocities from 90 to 190 m/s. The low-pass
% reduces controller activity above the actuator bandwidth.
%
% Three "soft goals" or objectives limit the modal response to disturbances and the flap demand
% caused by disturbances and sensor noise. The control-effort weight spans 20 to 150 rad/s,
% around the flutter frequency of 70 rad/s. A "hard goal" or constraint requests a disk margin
% of 6 dB and 45 degrees at the plant input. This is the chosen design target: the margin
% measures tolerance to simultaneous gain and phase variations (disk margin) before instability.
%
% The script seeds the random number generator, runs three random starts, asserts stability on
% the 42-point grid and saves the result. The paper scripts apply related tuning workflows to
% the eight-surface wing with blending structures.
%
% **Source:** |build/build_P_Goland.m|, |build/build_ABCD_P.m|,
% |R00_Synthesize_Goland_Controller.m|.

P = build_P_Goland(100,1,1,[1,2]);              % generalized plant at 100 m/s, modes 1 and 2 as performance outputs (runs its own DLM solve)
disp(P.InputName'), disp(P.OutputName')
fprintf('%6s %14s %14s\n', 'mode', 'q_f scale', 'q_f_dot scale')
for i = 1:num_modes
    fprintf('%6d %14.4g %14.4g\n', i, sqrt(Kff(i,i)/Kff(1,1)), sqrt(Mff(i,i)/Kff(1,1)))
end
fprintf('u_z_ddot scale Kff11/Mff11 = %.4g (= omega_1^2)\n', Kff(1,1)/Mff(1,1))

%% Where to go next
% **The papers**
%
% |research-paper-code/| holds one folder per paper, the shipped controllers and the runtimes.
% The |R*| scripts are the synthesis and sweep runs, the |fig*| scripts regenerate the figures
% from the shipped data.
%
% - **Open-Source Benchmark Model for Active Flutter Suppression** (AIAA SciTech 2026): this
%   model, its flutter mechanism, and a static gain from two accelerometers to flap 4 as the
%   first controller.
% - **Modal Blending for Active Flutter Suppression** (AIAA SciTech 2026): all eight
%   accelerometers blended into a few signals that look like the modes, so that a low-order
%   controller on flap 4 sees the flutter mode and little else.
% - **Structural Blending for Active Flutter Suppression** (IFASD 2026): blending vectors taken
%   from the structural mode shapes on both sides, eight accelerometers in and eight surfaces
%   out, with a geometric mode-isolation study and a sensor-failure test.
% - **Subspace Geometry and Performance Limits of Modal Blending for AFS** (Aerospace Systems,
%   to be released): what blending can and cannot achieve, derived from the geometry of the
%   modal subspaces.
%
% **Adding a wing**
%
% Copy |define/define_Goland_Structure_Aero.m|, change the numbers of chapter 1, adjust the flap
% panels and the hinge points to your grid. Then copy |build_G_Goland.m| and |build_P_Goland.m|
% and replace the define call and the surface count. The new wing then has entry points with the
% same signatures, so |getEigenvalueModeshape|, |Vg_plot|, |simulate_afs_switch| and
% |animate_wing| work on it unchanged.
%
% **Citing**
%
% If FlutterWing helps your research, cite the benchmark paper: J. Eichelsdoerfer, *Open-Source
% Benchmark Model for Active Flutter Suppression*, AIAA SciTech Forum 2026, AIAA 2026-1555,
% [doi:10.2514/6.2026-1555](https://doi.org/10.2514/6.2026-1555); the README carries the BibTeX.
%
% **Disclaimer:** I wrote this tutorial in its original form in german language. The tutorial
% was then translated to english and extended with the help of generative AI.
