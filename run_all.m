%RUN_ALL  Reproduce the numerical results of
%   "A state-dependent dynamical model of sensory regulation during
%    autism therapy" (Sections 5.1 and 7, Discussion).
%   Runs in MATLAB (R2016b or later) and GNU Octave (>= 6).
%   Prints each computed value next to the value reported in the paper,
%   and saves Figures 5-9 as PNG files in ../figures.
%
%   Set SAVE_FIGS = false to skip the figures, and NSAMP to change the
%   number of joint samples of the sensitivity analysis.

if ~exist('SAVE_FIGS', 'var'), SAVE_FIGS = true; end
if ~exist('NSAMP', 'var'),     NSAMP = 2000;     end
figdir = fullfile(fileparts(mfilename('fullpath')), '..', 'figures');
if SAVE_FIGS && ~exist(figdir, 'dir'), mkdir(figdir); end
rep = @(txt, val, paper) fprintf('  %-58s %9.3f   (paper: %s)\n', txt, val, paper);

fprintf('\n=== Section 5.1: illustrative example (Table 7, Figure 4) ===\n');
ex = example_inputs(false);
rep('Child A, t0: SH', ex.A_slow_t0(1), '0.61');
rep('Child A, t0: SZ', ex.A_slow_t0(2), '0.10');
rep('Child B, t0: SH', ex.B_slow_t0(1), '0.26');
rep('Child B, t0: SZ', ex.B_slow_t0(2), '0.40');
rep('Child A, t0: sum of aversive contributions', ex.A_sum_t0, '0.73');
rep('Child B: registration reduction, crying', ex.reg_crying, '28%');
rep('Child B: registration reduction, door slam', ex.reg_slam, '8%');
rep('Child A, t1: SH', ex.A_slow_t1(1), '0.27');
rep('Child A, t1: SZ', ex.A_slow_t1(2), '0.07');
rep('Child B, t1: SH', ex.B_slow_t1(1), '0.22');
rep('Child B, t1: SZ', ex.B_slow_t1(2), '0.49');
rep('Child B, t1: visual aversive contribution', ex.E_B(5,1), '0.16');
rep('Wall: Child B SH', ex.B_wall_t1(1), '0.10');
rep('Wall: Child B SZ', ex.B_wall_t1(2), '0.48');
rep('Wall: Child A SH', ex.A_wall_t1(1), '0.21');
rep('t0: Asigma (both children)', ex.A_slow_t0(3), '0.75');
rep('Child B, t1: organizing share', ex.B_t1_organizing, '50%');
rep('Child B, t1: calming share', ex.B_t1_calming, '46%');
rep('Child B, t1: alerting share', ex.B_t1_alerting, '4%');
rep('Child B, t1: Asigma (slow swing)', ex.B_slow_t1(3), '0.29');
rep('Child B, t1: Asigma (fast swing)', ex.B_fast_t1(3), '0.75');

fprintf('\n=== Section 7.2: two children in the same session (Figure 5) ===\n');
s = fig_session(SAVE_FIGS);
if SAVE_FIGS, save_png(fullfile(figdir, 'fig5_session.png')); end
rep('Child A: peak arousal, auditory events', s.A.aud_Amax, '0.60');
rep('Child A: minimum capacity, auditory events', s.A.aud_Rmin, '0.79');
rep('Child A: capacity at minute 15', s.A.R_at15, '0.84');
rep('Child B: resting arousal', s.B.rest, '0.30');
rep('Child B: peak arousal, auditory events', s.B.aud_Amax, '0.50');
rep('Child B: minimum capacity, auditory events', s.B.aud_Rmin, '0.93');
rep('Child B: peak participation, auditory events', s.B.aud_Pmax, '0.69');
rep('Child A: mean participation, slow swing', s.A.slow_Pmean, '0.83');
rep('Child A: mean participation, fast swing', s.A.fast_Pmean, '0.64');
rep('Child B: mean participation, slow swing', s.B.slow_Pmean, '0.19');
rep('Child B: mean arousal, fast swing', s.B.fast_Amean, '0.55');
rep('Child B: mean participation, fast swing', s.B.fast_Pmean, '0.54');

fprintf('\n=== Section 7.3: bistability and hysteresis (Figure 6) ===\n');
b = fig_bifurcation(SAVE_FIGS);
if SAVE_FIGS, save_png(fullfile(figdir, 'fig6_bifurcation.png')); end
rep('Child A: lower end of bistable range (eta = 2.5)', b.range_A(1), '0.15');
rep('Child A: upper end of bistable range (eta = 2.5)', b.range_A(2), '0.38');
rep('Child B: lower end of bistable range (eta = 2.5)', b.range_B(1), '0.47');
rep('Child B: upper end of bistable range (eta = 2.5)', b.range_B(2), '0.74');
rep('Child A: lower end (eta = 0)', b.range_eta0_A(1), '0.06');
rep('Child A: upper end (eta = 0)', b.range_eta0_A(2), '0.35');
rep('Child B: lower end (eta = 0)', b.range_eta0_B(1), '0.11');
rep('Child B: upper end (eta = 0)', b.range_eta0_B(2), '0.64');
rep('Child A: ramp up, leaves regulated state at SH', b.hyst_up, '0.39');
rep('Child A: ramp down, returns at SH', b.hyst_down, '0.13');
rep('Child B: regulated participation peaks at SH', b.B_Ppeak_SH, '0.39');

% Proposition 6.3 (eta = 0): explicit criterion versus numerical count
P0 = params_default(); P0.eta = 0;
mism = 0; Sg = linspace(0.005, 1, 400);
for nm = {'A', 'B'}
    C = child_profile(nm{1});
    for S = Sg
        o = prop63_explicit(S, 0, P0.Ast, C, P0);
        n = count_zeros(S, C, P0, linspace(0, 1, 20001));
        mism = mism + (o.bistable ~= (n == 3));
    end
end
rep('Prop. 6.3: mismatches explicit vs numerical (800 cases)', mism, '0');

fprintf('\n=== Section 7.4: memory of the session (Figure 7) ===\n');
m = fig_memory(SAVE_FIGS);
if SAVE_FIGS, save_png(fullfile(figdir, 'fig7_memory.png')); end
rep('Every 3 min: capacity before the last stimulus', m.close.Rbefore(6), '0.08');
rep('Every 3 min: first peak of arousal', m.close.peaks(1), '0.65');
rep('Every 3 min: last peak of arousal', m.close.peaks(6), '0.87');
rep('Every 3 min: second minus first peak', m.close.peaks(2) - m.close.peaks(1), '0.02');
rep('Every 15 min: lowest capacity before a stimulus', min(m.spaced.Rbefore), '>0.8');
rep('Every 15 min: capacity after the first stimulus', m.spaced.Rafter(1), '0.86');
rep('Every 15 min: capacity after the last stimulus', m.spaced.Rafter(6), '0.61');
rep('Every 15 min: lowest peak', min(m.spaced.peaks), '0.65');
rep('Every 15 min: highest peak', max(m.spaced.peaks), '0.68');

fprintf('\n=== Section 7.5: timing of the intervention (Figures 8 and 9) ===\n');
tm = fig_timing(SAVE_FIGS);
if SAVE_FIGS, save_png(fullfile(figdir, 'fig8_timing.png')); end
rep('Regulated state: A', tm.reg(1), '0.58');
rep('Regulated state: R', tm.reg(2), '0.53');
rep('Regulated state: P', tm.reg(3), '0.50');
rep('Dysregulated state: A', tm.dys(1), '0.77');
rep('Dysregulated state: R', tm.dys(2), '0.02');
rep('Shortest burst that dysregulates (min)', tm.burst_threshold, '0.85');
rep('Shortest effective quiet break after 3-min burst (min)', tm.min_break, '4');
rep('Capacity after the minimal break', tm.R_after_break(1), '0.14');
rep('Capacity after a 5-min break', tm.R_after_break(2), '0.17');
rep('Capacity after a 10-min break', tm.R_after_break(3), '0.36');
rep('90% of R_reg, 5-min break (min after burst)', tm.t90_R(2), '23-28');
rep('90% of R_reg, 10-min break (min after burst)', tm.t90_R(3), '23-28');
rep('90% of P_reg, 5-min break (min after burst)', tm.t90_P(2), '35-40');
rep('90% of P_reg, 10-min break (min after burst)', tm.t90_P(3), '35-40');
rep('Bound of Prop. 6.11 for 90% of R_reg (min)', tm.bound_t90, '13');

w = fig_window(SAVE_FIGS);
if SAVE_FIGS, save_png(fullfile(figdir, 'fig9_window.png')); end
o = w.offsets;  M = w.meanP;
rep('Mean participation, no burst', w.no_burst, '0.50');
rep('Mean participation, burst without control', w.no_control, '0.07');
rep('Full control, best onset in [-3, 0]', max(M(1, o >= -3 & o <= 0)), '0.55');
rep('Full control starting 3 min before the burst', M(1, o == -3), '0.55');
rep('Full control, delay 1 min', M(1, o == 1), '0.44');
rep('Full control, delay 20 min', M(1, o == 20), '0.21');
rep('Half control, best onset', max(M(2,:)), '0.48');
ok = o(M(2,:) > 0.3);
rep('Half control effective from onset', min(ok), '-2.5');
rep('Half control effective up to onset', max(ok), '0.25');

fprintf('\n=== Discussion: eigenvalues ===\n');
e = eigen_discussion();
rep('Dysregulated state, slowest eigenvalue at SH = 0.2', e.dys_slowest_200, '-0.61');
rep('Dysregulated state, slowest eigenvalue at SH = 0.155', e.dys_slowest_155, '-0.14');
rep('Regulated state near upper end: eigenvalue 1', e.reg_eig_upper_end(1), '-3.4');
rep('Regulated state near upper end: eigenvalue 2', e.reg_eig_upper_end(2), '-0.12');

fprintf('\n=== Section 7.1: sensitivity of the bistable range ===\n');
sensitivity(NSAMP, 1, true);
fprintf('  (paper: bistable in all one-at-a-time cases; (0.22, 0.48) and (0.10, 0.25);\n');
fprintf('   joint: 100%% for A, 91%% for B; median widths 0.20 and 0.22)\n\n');
