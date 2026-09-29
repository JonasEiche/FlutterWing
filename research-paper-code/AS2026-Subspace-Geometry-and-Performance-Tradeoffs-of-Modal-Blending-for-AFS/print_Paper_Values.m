%% print_Paper_Values.m
% Prints every computed number of the published AS2026 paper (Aerospace
% Systems, "Subspace Geometry and Performance Tradeoffs of Modal Blending
% for Active Flutter Suppression") from the data in Data/, grouped by the
% place where it is printed. The figures come from the figNN_* scripts.
%
% Reads (all shipped; regenerate with the R scripts named in brackets)
%   Data/Locked_Comparison_Results.mat       [R06_Locked_Comparison, writes
%                                             the *_test.mat variant]
%   Data/Modal_Obsrv_Contr_data.mat          [R01_Modal_Observability]
%   Data/Modal_Obsrv_Contr_noAct_data.mat    [R02_Modal_Observability_NoActuators]
%   Data/Geometric_Isolation.mat             [R03_Geometric_Mode_Isolation]
%   Data/Quadrature_Mismatch_DeltaSweep.mat  [R05_Quadrature_Mismatch]
% The toy-oscillator values of Sect. 4 are stored in the R06 archive, which
% takes the SISO sweep from Data/RHP_zeros_analysis.mat [R04_SISO_RHP_Zeros].
%


clearvars

archive = 'Locked_Comparison_Results.mat';

script_dir = fileparts(mfilename('fullpath'));
data_dir   = fullfile(script_dir, 'Data');
A   = load(fullfile(data_dir, archive));
M   = load(fullfile(data_dir, 'Modal_Obsrv_Contr_data.mat'));
N   = load(fullfile(data_dir, 'Modal_Obsrv_Contr_noAct_data.mat'));
iso = load(fullfile(data_dir, 'Geometric_Isolation.mat'));
iso = iso.results;

sel  = @(id) A.synth.(id).selected;
chan = @(id) A.model.cases(strcmp({A.model.cases.id}, id)).channel;
iFlut2D = 7;  iRes2D = 8;  iCrit4D = 9;     % set indices of R01
fprintf('Archive: %s\n', archive);

%% Sect. 2: model and operating points
% Open-loop onset on the V-g grid (first point with Re > 0.1 and Im > 1)
Vs = A.vg.Vsweep;
iOn = find(Vs == A.vg.onsetOL, 1);
EVon = A.vg.EV_OL(:, iOn);
pOn = EVon(real(EVon) > 0.1 & imag(EVon) > 1);
fprintf('\n=== Sect. 2 ===\n');
fprintf('Open-loop flutter onset          %.2f m/s  (paper: ~105 m/s)\n', A.vg.onsetOL);
fprintf('Coalescence frequency at onset   %.2f Hz = %.1f rad/s  (paper: 4.5 Hz, 28 rad/s)\n', ...
    abs(pOn(1))/(2*pi), abs(pOn(1)));

[~, i130] = min(abs(M.V_inf_arr - 130));
G = build_G_RectWing(M.V_inf_arr(i130), 1:8, 1:8);
p = eig(G.A);
pFlut = p(real(p) > -1  & abs(imag(p)) > 5 & abs(imag(p)) < 40 & imag(p) > 0);
pCrit = p(real(p) > -10 & abs(imag(p)) > 5 & abs(imag(p)) < 40 & imag(p) > 0);
pRes  = setdiff(pCrit, pFlut);
fprintf('Flutter pole at %.2f m/s        %.1f %+.1fj  (paper: 1.4 +/- 26.8j)\n', ...
    M.V_inf_arr(i130), real(pFlut), imag(pFlut));
fprintf('Residual pole at %.2f m/s       %.1f %+.1fj  (paper: -5.2 +/- 23.1j)\n', ...
    M.V_inf_arr(i130), real(pRes), imag(pRes));

%% Sect. 3: subspace geometry (Fig. 3, Table 4)
V = M.V_inf_arr;
mask = V > 106;
pB_flut = polyfit(log(V(mask)), log(M.SigMin_B(iFlut2D, mask)), 1);
pB_res  = polyfit(log(V(mask)), log(M.SigMin_B(iRes2D,  mask)), 1);
fprintf('\n=== Sect. 3 ===\n');
fprintf('Basis principal angle            1e%d deg  (paper: < 1e-13 deg)\n', ...
    ceil(log10(max([A.equivalence.checks.principalAngleDeg]))));
fprintf('sigma_min(U) growth exponent     flutter %.2f, residual %.2f  (paper: 1.5, 4.3)\n', ...
    pB_flut(1), pB_res(1));
fprintf('Min principal angle              full state %.2f deg, 8 sensors %.2f deg  (paper: 37.91, 0.16)\n', ...
    iso.full.minPrincipalAngleDeg, iso.real.minPrincipalAngleDeg);

fprintf('\nTable 4 (isolation)   sigma_min  projected  loss    min angle\n');
for k = {'full', 'real'}
    r = iso.(k{1});
    fprintf('  %-18s %8.3f %9.3f %6.1f%% %9.2f deg\n', r.name, ...
        r.sigminCflut, r.sigminCflutProj, 100*r.signalLoss, r.minPrincipalAngleDeg);
end
fprintf('Noise amplification (8 sensors)  %.0f  (paper: 335)\n', 1/iso.real.sigminCflutProj);

%% Sect. 4: toy oscillator (Figs. 4-5)
fprintf('\n=== Sect. 4 ===\n');
fprintf('Best SISO Delta                  %.0f deg  (paper: 18)\n', A.toy.bestSisoDeltaDeg);
fprintf('Best SISO zero                   %.1f  (paper: -1.9)\n', A.toy.bestSisoZ);
fprintf('Best SISO soft goal              %.3f  (paper: 0.137)\n', A.toy.bestSisoSoft);
fprintf('MIMO 2x2 soft goal               %.3f  (paper: 0.062)\n', A.toy.MaxSoft);
fprintf('SISO / MIMO ratio                %.1f  (paper: 2.2)\n', A.toy.ratioSisoOverMimo);
fprintf('H2 blending Delta, zero          %.1f deg, %.1f  (paper: -40.3, 12.2)\n', ...
    A.toy.h2DeltaDeg, A.toy.h2Z);

%% Sect. 5: blending comparison (Table 2, Fig. 6)
fprintf('\n=== Sect. 5 ===\n');
fprintf('Table 2                          MaxSoft   Hard\n');
rows = {'MB2orth', 'Modal MIMO (orthonormalized)'; ...
        'H2Pusch', 'H2-optimal SISO'; ...
        'MVFz0',   'Zero-aware SISO (z = 0)'; ...
        'R1bestRestricted', 'Best restricted direction'};
for i = 1:size(rows, 1)
    s = sel(rows{i, 1});
    if A.synth.(rows{i, 1}).exhausted || s.gHard > 1
        fprintf('  %-30s %7s  %6.3f\n', rows{i, 2}, '>> 1', s.gHard);
    else
        fprintf('  %-30s %7.3f  %6.3f\n', rows{i, 2}, s.MaxSoft, s.gHard);
    end
end
fprintf('  %-30s %7s  %6.4f\n', 'Random 2x2', '>> 1', sel('random').gHard);
fprintf('(paper: 0.645/0.997, 1.149/1.000, >>1/3.651, 1.128/1.000, >>1/0.9998)\n');

famSweep = load(fullfile(data_dir, 'Quadrature_Mismatch_DeltaSweep.mat'));
famSweep = famSweep.results;
okFam = famSweep.synth_hardGoal <= 1;
fprintf('Best restricted Delta            %.0f deg  (paper: 115)\n', A.model.family.bestFamilyDeltaDeg);
fprintf('H2 pair family Delta             %.0f deg  (paper: 156)\n', chan('H2Pusch').familyDisplayDeltaDeg);
fprintf('MVF z=0 pair family Delta        %.0f deg  (paper: ~65)\n', chan('MVFz0').familyDisplayDeltaDeg);
fprintf('Family worst/best MaxSoft        %.0f  (paper: 3252)\n', ...
    max(famSweep.synth_maxSoft(okFam)) / min(famSweep.synth_maxSoft(okFam)));
fprintf('H2 / MIMO MaxSoft ratio          %.2f  (paper: 1.78)\n', ...
    sel('H2Pusch').MaxSoft / sel('MB2orth').MaxSoft);
fprintf('MVF z=0 zero at 130 m/s          %+.2f rad/s  (paper: +50.66)\n', chan('MVFz0').at130.z);
fprintf('Tunable parameters               rank one %d, rank two %d  (paper: 19, 30)\n', ...
    sel('H2Pusch').nTunableParams, sel('MB2orth').nTunableParams);
fprintf('Closed-loop onset (Fig. 8)       MB2orth %s, H2Pusch %s  (paper: stable to 160 m/s)\n', ...
    onsetText(A.vg.MB2orth.onset), onsetText(A.vg.H2Pusch.onset));

%% Appendix A/B: Gram determinant, corruption (Table 3), conditioning (Fig. 7)
fprintf('\n=== Appendices A and B ===\n');
fprintf('4D output map at %.0f m/s         det %.1e, sigma_min %.3f, kappa %.0f  (paper: 7.4e6, 0.057, ~7300)\n', ...
    V(1), M.GramRaw_C(iCrit4D, 1), M.SigMin_C(iCrit4D, 1), M.Cond_C(iCrit4D, 1));

base = sel('MB2orth').MaxSoft;
cm = A.corruption.meta;
fprintf('Table 3                  thetaU  thetaY  MaxSoft  Degradation\n');
fprintf('  %-22s %4.0f %7.0f %8.3f  %s\n', 'Baseline', 0, 0, base, '---');
corrRows = {'corrU', 'U-corrupted', cm.thetaU_deg, 0; ...
            'corrY', 'Y-corrupted', 0, cm.thetaY_deg; ...
            'corrUY', 'Both-corrupted', cm.thetaU_deg, cm.thetaY_deg};
for i = 1:size(corrRows, 1)
    s = A.corruption.state.(corrRows{i, 1}).selected;
    fprintf('  %-22s %4.0f %7.0f %8.3f  %+.1f%%\n', corrRows{i, 2}, ...
        corrRows{i, 3}, corrRows{i, 4}, s.MaxSoft, 100*(s.MaxSoft - base)/base);
end
fprintf('(paper: 0.645, 0.659/+2.2%%, 1.125/+74.5%%, 1.133/+75.7%%)\n');
fprintf('Both-corrupted hard goal         %.5f  (paper: 0.99996)\n', ...
    A.corruption.state.corrUY.selected.gHard);

pN_flut = polyfit(log(N.V_inf_arr(N.V_inf_arr > 106)), log(N.SigMin_B(1, N.V_inf_arr > 106)), 1);
pN_res  = polyfit(log(N.V_inf_arr(N.V_inf_arr > 106)), log(N.SigMin_B(2, N.V_inf_arr > 106)), 1);
fprintf('kappa(C) at %.0f m/s              flutter %.0f, residual %.0f  (paper: ~15, ~26)\n', ...
    V(1), M.Cond_C(iFlut2D, 1), M.Cond_C(iRes2D, 1));
fprintf('kappa(C) at %.0f m/s             flutter %.1f, residual %.1f  (paper: 3 to 4)\n', ...
    V(end), M.Cond_C(iFlut2D, end), M.Cond_C(iRes2D, end));
fprintf('alpha_align = alpha - 2          flutter %.1f, residual %.1f  (paper: -0.5, +2.3)\n', ...
    pB_flut(1) - 2, pB_res(1) - 2);
fprintf('No-actuator growth exponents     flutter %.1f, residual %.1f  (paper: 1.6, 3.6)\n', ...
    pN_flut(1), pN_res(1));
post = V >= V(find(M.Cond_B(iFlut2D, :) == min(M.Cond_B(iFlut2D, :)), 1));
fprintf('kappa(U_flut) post-coalescence   %.1f -> %.1f  (paper: 3 -> 9)\n', ...
    min(M.Cond_B(iFlut2D, post)), M.Cond_B(iFlut2D, end));
fprintf('kappa(U_res) at %.0f m/s         %.1f  (paper: ~2)\n', V(end), M.Cond_B(iRes2D, end));

%% Appendix C: basis representations (Table 5), equivalence, static gains, zeros
fprintf('\n=== Appendix C ===\n');
fprintf('Table 5                  MaxSoft   Hard\n');
for b = {'MB2orth', 'Modal orth'; 'MB2raw', 'Modal raw'; 'MB2pinv', 'Modal pinv'}'
    fprintf('  %-22s %7.3f  %6.3f\n', b{2}, sel(b{1}).MaxSoft, sel(b{1}).gHard);
end
fprintf('(paper: 0.645/0.997, 0.649/0.995, 0.645/1.000)\n');
fprintf('Equivalent-controller deviation  1e%d  (paper: < 1e-13)\n', ...
    ceil(log10(max([A.equivalence.checks.closedLoopGap]))));

st = A.statics;
if isfield(st, 'min') && ~isempty(st.min)
    fprintf('Static MVF gain lambda*          %.2f  (paper: -1.87)\n', st.min.lambda);
    fprintf('Flutter-band damping             %.3f  (paper: 0.034)\n', st.min.diag.minBandDamping);
    fprintf('Max Re of the other poles        %.2f  (paper: -0.96)\n', st.min.eval.maxRealOther);
    fprintf('Input disk margins               %.1f dB / %.1f deg  (paper: 3.8 dB / 24.4 deg)\n', ...
        st.min.eval.dmGainDB, st.min.eval.dmPhaseDeg);
    fprintf('Tuned static gain on [%g, %g]    feasible: %d  (paper: none meets 6 dB / 45 deg)\n', ...
        st.tuned.lambdaGrid(1), st.tuned.lambdaGrid(end), st.tuned.feasible);
else
    fprintf('Static-gain designs not in this archive (needs the Robust Control Toolbox).\n');
end
fprintf('H2 channel zero                  %.1f rad/s at 90 m/s, %.1f at 130 m/s  (paper: -207.3, -325.4)\n', ...
    chan('H2Pusch').at90.z, chan('H2Pusch').at130.z);
fprintf('Best restricted channel zero     %.1f rad/s at 90 m/s, %.1f at 130 m/s  (paper: -354.5, -460.3)\n', ...
    chan('R1bestRestricted').at90.z, chan('R1bestRestricted').at130.z);

function s = onsetText(V)
    if isinf(V), s = 'none in [20,160]'; else, s = sprintf('%.2f m/s', V); end
end
