%% FlutterWing quick start
% Builds the rectangular research wing and determines its flutter speed. Then closes the loop with the
% shipped controller and animates the result.
% The physics is explained in |docs/TUTORIAL.md| and |TUTORIAL.m|.

%% 1. Build the state space model
% |build_G_RectWing| returns the linear plant $G$. A rectangular wing with four trailing-edge flaps and
% four leading-edge slats. Colocated with the control surfaces are eight IMUs (here plain accelerometers in the vertical direction). 
% Flutter depends on airspeed, so $G$ is an |ss| array with one model per velocity.

V_inf  = linspace(20,160,32);        % freestream velocity grid (m/s)
imuIDX = 1:8;  ailIDX = 1:8;         % all 8 accelerometers, all 8 control surfaces installed.
G = build_G_RectWing(V_inf, imuIDX, ailIDX);   % 8x8x32 ss array, 56 states each


%% 2. Draw the wing
% The scene drawn by |wing_scene| comes from the model definition: 200 aerodynamic panels (10 chordwise, 20 spanwise),
% flaps 1 to 4 and slats 5 to 8 from root to tip, and an IMU at the centre of each hinge line.

[Structure, ~] = define_RectWing_Structure_Aero(5, 6);   % 5 structural modes, 6 aerodynamic lag poles
wing_scene(Structure, 'Surfaces', 'all', 'IMUs', 'all', 'Triad', 'both', 'Freestream', true, 'Name', 'rectwing_layout');

%% 3. The flutter speed
% The vibration modes frequency and a damping change with airspeed because of aerodynamic stiffness and damping forces. 
% Flutter starts where the damping of one aeroelastic mode crosses zero. |getEigenvalueModeshape|
% computes the eigenvalues of the plant at every velocity and tracks each mode. |Vg_plot| shows frequency
% (top) and damping (bottom) of the bending and torsion modes,(here hard-coded between 2 and 10 Hz).

EV = getEigenvalueModeshape(G, 5);                          % eigenvalues of every model, tracked mode by mode over velocity
[~, cr] = Vg_plot(V_inf, EV, {'open loop'}, 'Band', [2 10]); % V-g plot of the modes between 2 and 10 Hz
fprintf('Flutter onset: %.1f m/s at %.2f Hz \n', cr.V, cr.f)

%% 4. What flutter looks like
% Just above flutter onset, the coalescing bending and torsion modes are about an eighth of a cycle phase-shifted. 
% |simulate_afs_switch| simulates the (open loop) wing at 108 m/s. |animate_wing| plays the result.

G108   = build_G_RectWing(108, imuIDX, ailIDX);            % one model, 3.7 m/s past the onset
sim_ol = simulate_afs_switch(G108, [], 'GrowthCycles', 4, 'Verbose', false); % start in the flutter mode shape, no controller
fprintf('At 108 m/s the flutter mode oscillates at %.2f Hz and grows %.0f%% per cycle\n', sim_ol.f_Hz, 100*(sim_ol.growthPerCycle-1))
animate_wing(Structure, sim_ol.q_f, 'Skin', 'flat', 'Output', 'figure', 'Verbose', false);

%% 5. Close the loop with the shipped controller
% Active flutter suppression (AFS) measures the wing motion with the IMUs and moves the
% control surfaces to actively damp the oscillations.
% |data/RectWing_Cont_imu18_ail18.mat| holds a static 8 x 8 gain (scaled IMU signals in,
% surface commands out). The controller was tuned with |systune| see |R01_Synthesize_RectWing_Controller.m|. 
% |feedback(G,-K)| is equivalent to $u = K y$.

load RectWing_Cont_imu18_ail18.mat          % variable RectWing_Cont: 8 IMU signals in, 8 surface commands out
K  = RectWing_Cont;
CL = feedback(G, -K);                       % u = K*y, the sign convention the controller was tuned with
EV_cl = getEigenvalueModeshape(CL, 5);
Vg_plot({V_inf, V_inf}, {EV, EV_cl}, {'open loop', 'closed loop'}, 'Band', [2 10]);


%% 6. Simulate flutter suppression
% Time Simulation at 108 m/s. The AFS controller is switched on after 6 cycles of open-loop flutter. 
% The simulation shows the wing motion, the IMU accelerations, and the control surface deflections. 

sim = simulate_afs_switch(G108, K, 'GrowthCycles', 6, 'ControlCycles', 8, 'FramesPerCycle', 12, 'Verbose', false);
S = fw_style();                                              % colour tokens of the figure family
Readout = repmat({'$V_\infty = 108$ m/s'}, 1, sim.Nt);      % velocity readout, one string per frame
animate_wing(Structure, sim.q_f, 'Surfaces', sim.delta, 'IMUAccel', sim.acc, 'Active', sim.active, ...
    'Skin', 'foil', 'NACA', '2408', 'FlapChord', 0.25, 'SlatChord', 0.20, ...
    'TipAmplitude', 0.07, 'SurfaceAmplitude', deg2rad([33 25]), ...
    'Readout', Readout, 'ReadoutColor', S.muted, 'LabelOff', 'AFS OFF', 'LabelOn', 'AFS ON', ...
    'Output', 'figure', 'Verbose', false);

%% 7. Where to go next:
% - |docs/TUTORIAL.md| and |TUTORIAL.m| build the Goland benchmark wing.
% - |docs/CONVENTIONS.md| defines frames, signs, numbering and scaling.
% - |research-paper-code/| holds the synthesis workflows of four papers on this wing.
% - |R01_Synthesize_RectWing_Controller.m| retunes the shipped 8 x 8 gain with |systune|.
