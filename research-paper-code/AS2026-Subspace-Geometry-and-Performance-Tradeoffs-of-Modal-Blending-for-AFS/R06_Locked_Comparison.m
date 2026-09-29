%% R06: Locked comparison of the published paper
%
% Regenerates every comparison number of the published AS2026 paper: 
% Table 2 (blending methods), Table 3 (subspace corruption), 
% Table 5 (basis representations), the MIMO baseline of the
% toy oscillator, the static-gain designs and zeros of Appendix C, the
% basis-equivalence check and the V-g loci of Fig. 8.
%
% =========================================================================
% PROTOCOL
% =========================================================================
%
%   Blending directions are extracted at V_ref = 90 m/s (stable) from the
%   flutter pole pair (Re > -1, 5 < |Im| < 40 rad/s); every controller is
%   synthesized at V_synth = 130 m/s (unstable) with the same order-4
%   structured core and the same goals (Vp .2, Vn .1, Vd .5, Vu .5, Theis
%   weight 12/64, hard goal 6 dB / 45 deg at 'ud').
%
%
%   Cases (Table 2 and 5):
%     H2Pusch          rank one, unconstrained H2-optimal SISO directions
%     MVFz0            rank one, zero-aware (modal velocity feedback,
%                      r0 = 0) directions: induced zero at the origin at 90 m/s
%     R1bestRestricted rank one, best member (Delta = 115 deg, argmin of the
%                      R05 sweep) of the restricted family of R05
%     MB2orth / MB2raw / MB2pinv  rank two, flutter pole-vector bases
%     random           rank two, random directions (negative control)
%   The H2 and zero-aware directions are hard-coded (their code is
%   proprietary, as for the H2 vectors of R04 and R05); all other directions
%   are computed here.
%
% =========================================================================
% STAGES
% =========================================================================
%
%   1. Model        plants, modal data, case directions, channel zeros
%   2. Synthesis    the 7 cases
%   3. Corruption   MB2orth basis with one column replaced by a random
%                   null-space vector: U, Y and both
%   4. Toy MIMO     2x2 order-2 controller on the R04 oscillator
%   5. Statics      static modal-velocity-feedback gain at 130 m/s and a
%                   gain sweep lambda = -6:0.05:6 on the MVFz0 directions;
%                   needs diskmargin (Robust Control Toolbox), skipped
%                   without it
%   6. Equivalence  the MB2orth controller transformed exactly to the raw
%                   and pinv bases, closed-loop deviation
%   7. V-g loci     open loop and the principal cases on linspace(20,160,32)
%
% =========================================================================
% DEPENDENCIES
% =========================================================================
%
%   build_G_RectWing, build_P_RectWing, getEigenvalueModeshape,
%   h2_opt_output_siso, h2_opt_input_siso
%   Data/Quadrature_Mismatch_DeltaSweep.mat  (R05: argmin of the family)
%   Data/RHP_zeros_analysis.mat              (R04: toy oscillator)
%   Control System Toolbox; Robust Control Toolbox for stage 5 only
%
% See also: print_Paper_Values, fig06_Quadrature_Mismatch_Sweep,
%           fig08_Vg_Diagram, R04_SISO_RHP_Zeros, R05_Quadrature_Mismatch

clearvars

script_dir = fileparts(mfilename('fullpath'));
data_dir   = fullfile(script_dir, 'Data');
out_file   = fullfile(data_dir, 'Locked_Comparison_Results_test.mat');

%% Protocol settings
V_ref   = 90;     % blending extraction (stable)
V_synth = 130;    % synthesis and evaluation (unstable)
imuIDX = 1:8;
ailIDX = 1:8;
modesIDX = 1:2;

band = struct('reMin', -1, 'imAbsMin', 5, 'imAbsMax', 40);   % flutter pole pair

cfg = struct();
cfg.nO       = 4;                      % order of the structured core
cfg.seeds    = [20260305 20260716];
cfg.ladder   = [3 7 15];               % RandomStart levels
cfg.agreeTol = 0.02;                   % relative cross-seed agreement

% Open a local pool if the Parallel Computing Toolbox is installed; without it
% systune ignores 'UseParallel' and evaluates the random starts serially (slower).
if ~isempty(ver('parallel')) && isempty(gcp('nocreate')), parpool('Processes'); end
cfg.opt = @(RS) systuneOptions('RandomStart', RS, 'UseParallel', true, 'Display', 'off');

% Tuning goals (the same for every case)
w1 = 12;  w2 = 64;
s = tf('s');
invWu = ((s+0.01*w1)*(0.01*s+w2)) / ((s+w1)*(s+w2));
Vp = 0.2;  Vn = 0.1;  Vd = 0.5;  Vu = 0.5;
noise_u_z_ddot = {'noise_u_z1_ddot','noise_u_z2_ddot','noise_u_z3_ddot','noise_u_z4_ddot', ...
                  'noise_u_z5_ddot','noise_u_z6_ddot','noise_u_z7_ddot','noise_u_z8_ddot'};
flap_d_slat_d = {'flap1_d','flap2_d','flap3_d','flap4_d', ...
                 'slat1_d','slat2_d','slat3_d','slat4_d'};
ReqAttenQF      = TuningGoal.Gain({'dist_q_f1_ddot','dist_q_f2_ddot'},{'q_f1','q_f2'},Vp*(1/Vd));
ReqContEffNoise = TuningGoal.Gain(noise_u_z_ddot,flap_d_slat_d,invWu*Vu*(1/Vn));
ReqContEffDist  = TuningGoal.Gain({'dist_q_f1_ddot','dist_q_f2_ddot'},flap_d_slat_d,invWu*Vu*(1/Vd));
cfg.soft = [ReqAttenQF, ReqContEffNoise, ReqContEffDist];
cfg.hard = TuningGoal.Margins('ud',6,45);

%% 1. Model: plants, modal data, case directions
G90  = build_G_RectWing(V_ref,   imuIDX, ailIDX);
G130 = build_G_RectWing(V_synth, imuIDX, ailIDX);
P130 = build_P_RectWing(V_synth, imuIDX, ailIDX, modesIDX);
[ny, nu] = size(G90);

[V90, DD90] = eig(G90.A);
d90  = diag(DD90);
f90  = find(real(d90) > band.reMin & abs(imag(d90)) > band.imAbsMin & abs(imag(d90)) < band.imAbsMax);
Bm90 = V90 \ G90.B;
Cm90 = G90.C * V90;
W90  = inv(V90)';                     % left eigenvectors
iPos90 = f90(imag(d90(f90)) > 0);

[V130, DD130] = eig(G130.A);
d130  = diag(DD130);
f130  = find(real(d130) > band.reMin & abs(imag(d130)) > band.imAbsMin & abs(imag(d130)) < band.imAbsMax);
Bm130 = V130 \ G130.B;
Cm130 = G130.C * V130;
iPos130 = f130(imag(d130(f130)) > 0);
assert(numel(f90) == 2 && numel(f130) == 2, 'Expected one flutter pole pair at 90 and 130 m/s.');

% H2-optimal SISO directions at 90 m/s (proprietary code, hence hard-coded)
ky_H2 = [-0.099544468735065;-0.320113541267249;-0.532226133992523;-0.69969646108878; ...
          0.0764433988501878;0.180002683267535;0.213140203844075;0.176366432097578];
ku_H2 = [0.0449869237599164;0.311964052487485;0.447659257628378;0.560302488725123; ...
         0.115582157524943;0.268015357431538;0.385831844752777;0.390204097174338];

% Zero-aware (r0 = 0) directions at 90 m/s (proprietary code, hence hard-coded)
ky_MVF = [0.0989719733301003;0.31999836929846;0.534285290258859;0.704999553646188; ...
          -0.0753951487248509;-0.175671339101149;-0.204756412340957;-0.163860959003604];
ku_MVF = [0.722332644055753;0.50729834241188;0.101583820006149;-0.442345366686718; ...
          -0.0455808086072744;-0.0131737292416505;0.0407249127088668;0.104812076741021];

% Restricted family of R05: output direction fixed at the H2-separate angle,
% input direction rotated by Delta; best member = argmin of the R05 sweep
Lfam = load(fullfile(data_dir, 'Quadrature_Mismatch_DeltaSweep.mat'));
famSweep = Lfam.results;
[bestSweepMaxSoft, iFam] = min(famSweep.synth_maxSoft);
bestFamilyDeltaDeg = famSweep.Delta_synth_deg(iFam);

cvr = unitNorm(real(Cm90(:, f90(1))));
cvi = unitNorm(imag(Cm90(:, f90(1))));
bwr = unitNorm(real(Bm90(f90(1), :))');
bwi = unitNorm(imag(Bm90(f90(1), :))');
[~, ky_h2sep] = h2_opt_output_siso(G90.C, V90(:, f90(1)));
[~, ku_h2sep] = h2_opt_input_siso(G90.B, W90(:, f90(1)));
ay = [cvr, cvi] \ ky_h2sep;
theta_y = atan2(ay(2), ay(1));
ky_fam  = unitNorm(cos(theta_y)*cvr + sin(theta_y)*cvi);
th_u    = theta_y + deg2rad(bestFamilyDeltaDeg);
ku_best = unitNorm(cos(th_u)*bwr + sin(th_u)*bwi);

% Rank-two pole-vector bases (raw, orthonormalized, pseudoinverse) and the random control
tol = 1e-8;
Cmr = [real(Cm90(:, f90(1))), imag(Cm90(:, f90(1)))];
Bmr = [real(Bm90(f90(1), :)); imag(Bm90(f90(1), :))]';
rng(20260305, 'twister');
ky_rand = rand(ny, 2)*2 - 1;  ky_rand = ky_rand ./ vecnorm(ky_rand);
ku_rand = rand(nu, 2)*2 - 1;  ku_rand = ku_rand ./ vecnorm(ku_rand);

% Channel data of the rank-one cases: residues R0, R1 of the flutter mode,
% induced zero z = -r0/r1 at 90 and 130 m/s, family coordinate Delta
[R0_90,  R1_90]  = modalResidues(d90(iPos90),   Bm90(iPos90, :),   Cm90(:, iPos90));
[R0_130, R1_130] = modalResidues(d130(iPos130), Bm130(iPos130, :), Cm130(:, iPos130));
chMeta = @(ku, ky) struct( ...
    'at90',  blendedModeChannel(R0_90,  R1_90,  d90(iPos90),   ku, ky), ...
    'at130', blendedModeChannel(R0_130, R1_130, d130(iPos130), ku, ky), ...
    'familyDisplayDeltaDeg', familyDisplayDelta(ku, ky, cvr, cvi, bwr, bwi));

cases = struct('id', {}, 'KU', {}, 'KY', {}, 'channel', {});
cases(end+1) = struct('id', 'H2Pusch', 'KU', ku_H2, 'KY', ky_H2, 'channel', chMeta(ku_H2, ky_H2));
cases(end+1) = struct('id', 'MVFz0', 'KU', ku_MVF, 'KY', ky_MVF, 'channel', chMeta(ku_MVF, ky_MVF));
cases(end+1) = struct('id', 'R1bestRestricted', 'KU', ku_best, 'KY', ky_fam, 'channel', chMeta(ku_best, ky_fam));
cases(end+1) = struct('id', 'MB2orth', 'KU', orth(Bmr, tol), 'KY', orth(Cmr, tol), 'channel', []);
cases(end+1) = struct('id', 'random',  'KU', ku_rand, 'KY', ky_rand, 'channel', []);
cases(end+1) = struct('id', 'MB2raw',  'KU', Bmr, 'KY', Cmr, 'channel', []);
cases(end+1) = struct('id', 'MB2pinv', 'KU', pinv(Bmr)', 'KY', pinv(Cmr)', 'channel', []);

model = struct('flutter90', d90(iPos90), 'flutter130', d130(iPos130), ...
    'residues', struct('R0_90', R0_90, 'R1_90', R1_90, 'R0_130', R0_130, 'R1_130', R1_130), ...
    'family', struct('theta_y_deg', rad2deg(theta_y), 'bestFamilyDeltaDeg', bestFamilyDeltaDeg, ...
                     'bestSweepMaxSoft', bestSweepMaxSoft), ...
    'cases', cases);
fprintf('Model: flutter pole %.4f%+.4fi at %d m/s, %.4f%+.4fi at %d m/s\n', ...
    real(model.flutter90), imag(model.flutter90), V_ref, ...
    real(model.flutter130), imag(model.flutter130), V_synth);
fprintf('Zero-aware channel zero at %d m/s: %.2e (0 by construction)\n', V_ref, cases(2).channel.at90.z);

%% 2. Synthesis: seed ladder over the 7 cases
synth = struct();
for j = 1:numel(cases)
    synth.(cases(j).id) = seedLadder(cases(j).id, P130, cases(j).KU, cases(j).KY, cfg);
end

%% 3. Corruption of the MB2orth basis (Table 3)
QU = cases(strcmp({cases.id}, 'MB2orth')).KU;
QY = cases(strcmp({cases.id}, 'MB2orth')).KY;
NU = null(QU');
NY = null(QY');
rng(20260306, 'twister');
rU = NU * randn(size(NU, 2), 1);  rU = rU / norm(rU);
rY = NY * randn(size(NY, 2), 1);  rY = rY / norm(rY);
corruption = struct();
corruption.meta = struct('rU', rU, 'rY', rY, ...
    'thetaU_deg', rad2deg(subspace([QU(:,1), rU], QU)), ...
    'thetaY_deg', rad2deg(subspace([QY(:,1), rY], QY)));
corruption.state.corrU  = seedLadder('corrU',  P130, [QU(:,1), rU], QY, cfg);
corruption.state.corrY  = seedLadder('corrY',  P130, QU, [QY(:,1), rY], cfg);
corruption.state.corrUY = seedLadder('corrUY', P130, [QU(:,1), rU], [QY(:,1), rY], cfg);

%% 4. Toy oscillator: 2x2 MIMO baseline (SISO sweep values from R04)
Ltoy = load(fullfile(data_dir, 'RHP_zeros_analysis.mat'));
toySweep = Ltoy.results;
Gt = ss(toySweep.A, toySweep.B, toySweep.C, 0);
[nyt, nut] = size(Gt);
ReqMargT  = TuningGoal.Margins('ud', 6, 45);
ReqAttenT = TuningGoal.Gain({'ym'}, {'ud'}, 10);
for RS = cfg.ladder
    ms = nan(1, numel(cfg.seeds));
    for iS = 1:numel(cfg.seeds)
        blk = tunableSS('cont', 2, nut, nyt);
        tuneCL = feedback(Gt, AnalysisPoint('ud', nut) * blk * AnalysisPoint('ym', nyt));
        opt = systuneOptions('RandomStart', RS, 'UseParallel', false, 'Display', 'off');
        rng(cfg.seeds(iS), 'twister');
        [~, fSoft] = systune(tuneCL, ReqAttenT, ReqMargT, opt);
        ms(iS) = max(fSoft);
        fprintf('  toy MIMO seed %d RS %2d: MaxSoft %.6g\n', cfg.seeds(iS), RS, ms(iS));
    end
    if (max(ms) - min(ms)) / max(ms) <= cfg.agreeTol, break, end
end
[bestSisoSoft, iBest] = min(toySweep.sweep_softGoal);
toy = struct('MaxSoft', min(ms), 'bestSisoSoft', bestSisoSoft, ...
    'bestSisoDeltaDeg', toySweep.Delta_deg_sweep(iBest), ...
    'bestSisoZ', toySweep.sweep_z_formula(iBest), ...
    'h2DeltaDeg', toySweep.h2_Delta_deg(1), ...
    'h2Z', toySweep.sweep_z_formula(toySweep.h2_sweep_idx(1)));
toy.ratioSisoOverMimo = toy.bestSisoSoft / toy.MaxSoft;

%% 5. Static modal-velocity-feedback gains (Appendix C)
statics = struct();
if ~isempty(ver('robust'))
    p130 = d130(iPos130);
    pTarget = -real(p130) + 1i*imag(p130);    % mirrored pair: omega_n kept, z = 0
    % Zero-aware directions of the 130 m/s model (proprietary code, hence hard-coded)
    ky_st = [-0.082929882550565;-0.181586161198802;-0.194791116577673;-0.133943296498533; ...
             0.0947456409205418;0.316824138407158;0.534985876806882;0.713232469348565];
    ku_st = [-0.0176723830016931;0.27035720120657;0.505990598654695;0.68541107379208; ...
             0.0894840071099137;0.199173784544251;0.277675980963137;0.275677910339359];
    r1_st = ky_st' * R1_130 * ku_st;
    % Active-damping gain: lambda = 2*(zeta*wn - zeta_cl*wn_cl)/r1 at unchanged wn
    lamStar = 2*(-real(p130) - (-real(pTarget))) / r1_st;
    evalMin = staticEval(G130, ku_st, ky_st, lamStar, pTarget);
    statics.min = struct('lambda', lamStar, 'eval', evalMin, 'diag', bandDiagnostics(evalMin));

    % Gain sweep on the zero-aware directions of 90 m/s
    lambdaGrid = -6:0.05:6;
    sweep = staticEval(G130, ku_MVF, ky_MVF, lambdaGrid, pTarget);
    statics.tuned = struct('feasible', any([sweep.marginOK]), 'lambdaGrid', lambdaGrid, ...
        'gridMinMaxRe', min(arrayfun(@(e) max(real(e.clPoles)), sweep)));
    fprintf('Static gain lambda* = %.4f: %.2f dB / %.2f deg; sweep feasible: %d\n', ...
        lamStar, evalMin.dmGainDB, evalMin.dmPhaseDeg, statics.tuned.feasible);
else
    fprintf('Stage 5 skipped: diskmargin needs the Robust Control Toolbox.\n');
end

%% 6. Basis equivalence of the MB2orth controller
w = logspace(-1, 3, 400);
sOrth = synth.MB2orth.selected;
KUo = cases(strcmp({cases.id}, 'MB2orth')).KU;  KYo = cases(strcmp({cases.id}, 'MB2orth')).KY;
CL1 = lft(P130, sOrth.K);
equivalence = struct('checks', struct('name', {}, 'controllerGap', {}, ...
    'closedLoopGap', {}, 'principalAngleDeg', {}));
for id = {'MB2raw', 'MB2pinv'}
    KU2 = cases(strcmp({cases.id}, id{1})).KU;
    KY2 = cases(strcmp({cases.id}, id{1})).KY;
    Tu = KUo \ KU2;
    Ty = KYo \ KY2;
    K2 = KU2 * (inv(Tu) * sOrth.Kcore * inv(Ty')) * KY2'; %#ok<MINV>
    equivalence.checks(end+1) = struct('name', id{1}, ...
        'controllerGap', peakRelGap(sOrth.K, K2, w), ...
        'closedLoopGap', peakRelGap(CL1, lft(P130, K2), w), ...
        'principalAngleDeg', max(rad2deg(subspace(KUo, KU2)), rad2deg(subspace(KYo, KY2))));
    fprintf('Equivalence MB2orth -> %s: closed-loop gap %.3g\n', id{1}, equivalence.checks(end).closedLoopGap);
end

%% 7. V-g loci (Fig. 8)
Vsweep = linspace(20, 160, 32);
Gs = build_G_RectWing(Vsweep, imuIDX, ailIDX);
num_modes = 5;
vg = struct('Vsweep', Vsweep);
[vg.EV_OL, ~] = getEigenvalueModeshape(Gs, num_modes);
vg.onsetOL = onsetVelocity(vg.EV_OL, Vsweep);
for id = {'H2Pusch', 'MVFz0', 'R1bestRestricted', 'MB2orth'}
    CL = feedback(Gs, -synth.(id{1}).selected.K);
    [EV, ~] = getEigenvalueModeshape(CL, num_modes);
    vg.(id{1}) = struct('EV', EV, 'onset', onsetVelocity(EV, Vsweep), ...
        'firstInstabilityV', firstInstabilityV(CL, Vsweep));
end
fprintf('Open-loop onset %.2f m/s\n', vg.onsetOL);

%% Save
save(out_file, 'model', 'synth', 'corruption', 'toy', 'statics', ...
     'equivalence', 'vg');
fprintf('Saved %s\n', out_file);

fprintf('\n%-18s %10s %10s %4s %s\n', 'Case', 'MaxSoft', 'Hard', 'RS', 'agreement');
for id = fieldnames(synth)'
    sc = synth.(id{1});
    tag = 'agreed';
    if sc.exhausted, tag = 'OPTIMIZER-LIMITED'; end
    fprintf('%-18s %10.4f %10.4f %4d %s\n', id{1}, sc.selected.MaxSoft, ...
        sc.selected.gHard, sc.selected.randomStart, tag);
end
disp('Paper values: set archive = ''Locked_Comparison_Results_test.mat'' in print_Paper_Values.m')

%% Local functions

function st = seedLadder(id, P, KU, KY, cfg)
% Run both seeds per RandomStart level until every soft goal and the hard
% goal agree within cfg.agreeTol; report the lower MaxSoft of that level.
runs = struct('seed', {}, 'randomStart', {}, 'fSoft', {}, 'gHard', {}, ...
              'MaxSoft', {}, 'nTunableParams', {}, 'K', {}, 'Kcore', {});
st = struct('agreed', false, 'exhausted', false, 'runs', [], 'selected', []);
for RS = cfg.ladder
    lv = runs([]);
    for seed = cfg.seeds
        blk = tunableSS('Cont', cfg.nO, size(KU, 2), size(KY, 2));
        nPar = nnz(blk.A.Free) + nnz(blk.B.Free) + nnz(blk.C.Free) + nnz(blk.D.Free);
        tuneCL = lft(P, AnalysisPoint('ud', size(KU, 1)) * (KU * blk * KY') * ...
                        AnalysisPoint('ym', size(KY, 1)));
        rng(seed, 'twister');
        tic;
        [CL, fSoft, gHard] = systune(tuneCL, cfg.soft, cfg.hard, cfg.opt(RS));
        Kcore = getBlockValue(CL, 'Cont');
        run = struct('seed', seed, 'randomStart', RS, 'fSoft', fSoft(:)', ...
            'gHard', gHard, 'MaxSoft', max(fSoft), 'nTunableParams', nPar, ...
            'K', KU * Kcore * KY', 'Kcore', Kcore);
        fprintf('  %-16s seed %d RS %2d: MaxSoft %.6g Hard %.6g (%.0f s)\n', ...
            id, seed, RS, run.MaxSoft, gHard, toc);
        lv(end+1) = run; %#ok<AGROW>
    end
    runs = [runs, lv]; %#ok<AGROW>
    m = cell2mat(arrayfun(@(r) [r.fSoft, r.gHard], lv(:), 'UniformOutput', false));
    relDiff = (max(abs(m), [], 1) - min(abs(m), [], 1)) ./ max(max(abs(m), [], 1), eps);
    [~, iMin] = min([lv.MaxSoft]);
    st.selected = lv(iMin);
    if all(relDiff <= cfg.agreeTol)
        st.agreed = true;
        break
    end
    if RS == cfg.ladder(end)
        st.exhausted = true;    % optimizer-limited: reported, never hidden
        fprintf('  %-16s OPTIMIZER-LIMITED (max rel diff %.3g)\n', id, max(relDiff));
    end
end
st.runs = runs;
end

function [R0, R1] = modalResidues(p, btilde, ctilde)
% Real residue matrices of the mode pair p, conj(p): the blended channel is
% (r1*s + r0)/(s^2 + 2*zeta*wn*s + wn^2) with r0 = ky'*R0*ku, r1 = ky'*R1*ku.
Rt = ctilde * btilde;
R0 = -2 * real(conj(p) * Rt);
R1 =  2 * real(Rt);
end

function ch = blendedModeChannel(R0, R1, p, ku, ky)
% Channel coefficients, induced zero z = -r0/r1 and quadrature angle Delta.
r0 = ky' * R0 * ku;
r1 = ky' * R1 * ku;
beta = (r0 + real(p)*r1) / abs(imag(p));
ch = struct('r0', r0, 'r1', r1, 'z', -r0/r1, 'DeltaDeg', atan2d(beta, r1));
end

function dd = familyDisplayDelta(ku, ky, cvr, cvi, bwr, bwi)
% R05 family coordinate: projection onto the normalized coupling bases,
% wrapped to (-180, 180] and displayed modulo 180 deg.
ay = [cvr, cvi] \ ky;
au = [bwr, bwi] \ ku;
d = atan2(au(2), au(1)) - atan2(ay(2), ay(1));
d = mod(rad2deg(d) + 180, 360) - 180;
dd = mod(d, 180);
end

function out = staticEval(G, ku, ky, lambda, pTarget)
% Static rank-one gain K = lambda*ku*ky' on the strictly proper plant G:
% closed-loop poles of A + B*K*C, the pair nearest |pTarget|, and the
% multi-loop input disk margins of L = -K*G.
wn = abs(pTarget);
for k = numel(lambda):-1:1
    K = lambda(k) * (ku * ky');
    clPoles = eig(G.A + G.B * K * G.C);
    [~, idx] = mink(abs(abs(clPoles) - wn), 2);
    rest = clPoles(setdiff(1:numel(clPoles), idx));
    [~, MM] = diskmargin(-K * G);
    out(k) = struct('lambda', lambda(k), 'clPoles', clPoles, ...
        'maxRealOther', max(real(rest)), ...
        'dmGainDB', 20*log10(MM.GainMargin(2)), 'dmPhaseDeg', MM.PhaseMargin(2), ...
        'marginOK', 20*log10(MM.GainMargin(2)) >= 6 && MM.PhaseMargin(2) >= 45);
end
end

function d = bandDiagnostics(e)
% Stability and the lowest damping of the oscillatory poles in 5 < Im < 40.
osc = e.clPoles(abs(imag(e.clPoles)) > 5 & abs(imag(e.clPoles)) < 40 & imag(e.clPoles) > 0);
d = struct('stabilized', max(real(e.clPoles)) < 0, 'maxRe', max(real(e.clPoles)), ...
           'minBandDamping', min(-real(osc) ./ abs(osc)));
end

function g = peakRelGap(S1, S2, w)
% Peak over w of ||S1 - S2|| relative to the peak of ||S1||.
H1 = freqresp(S1, w);
H2 = freqresp(S2, w);
num = zeros(numel(w), 1);
den = zeros(numel(w), 1);
for i = 1:numel(w)
    num(i) = norm(H1(:,:,i) - H2(:,:,i));
    den(i) = norm(H1(:,:,i));
end
g = max(num) / max(den);
end

function V = onsetVelocity(EV, Vsweep)
% First velocity with a pole of Re > 0.1 and Im > 1 (Inf if none).
idx = find(any(real(EV) > 0.1 & imag(EV) > 1, 1), 1);
if isempty(idx), V = Inf; else, V = Vsweep(idx); end
end

function V = firstInstabilityV(CL, Vsweep)
% Velocity where the largest real part of the closed-loop poles crosses
% zero (linear interpolation; Inf if stable over the sweep).
maxRe = arrayfun(@(iv) max(real(pole(CL(:, :, iv)))), 1:numel(Vsweep));
ix = find(maxRe(1:end-1) <= 0 & maxRe(2:end) > 0, 1);
if maxRe(1) > 0
    V = Vsweep(1);
elseif isempty(ix)
    V = Inf;
else
    V = Vsweep(ix) + (Vsweep(ix+1) - Vsweep(ix)) * (0 - maxRe(ix)) / (maxRe(ix+1) - maxRe(ix));
end
end

function x = unitNorm(x)
x = x / norm(x);
end
