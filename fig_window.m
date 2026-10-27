function out = fig_window(doplot, P, offsets)
%FIG_WINDOW  Section 7.5 / Figure 9: mean participation of Child A over
%   60 min as a function of the onset (relative to the onset of a 3-min
%   burst) of a 5-min control that multiplies the aversive input by 1-c.
if nargin < 1, doplot = true; end
if nargin < 2 || isempty(P), P = params_default(); end
if nargin < 3, offsets = -5:0.25:20; end
S0 = 0.22; Sb = 0.8; tb = 5; Lb = 3; Lc = 5; T = 60;
[C, xreg] = timing_setup(P, S0);
meanP = @(on, c, burst) mean_part(on, c, burst, S0, Sb, tb, Lb, Lc, T, xreg, C, P);
out.no_burst   = meanP(-100, 0, false);
out.no_control = meanP(-100, 0, true);
out.offsets = offsets;
cs = [1, 0.5];
for i = 1:2
    for j = 1:numel(offsets)
        out.meanP(i,j) = meanP(tb + offsets(j), cs(i), true);
    end
end
if doplot
    figure('Name', 'Figure 9', 'Position', [100 100 650 380]); hold on; box on;
    plot(offsets, out.meanP(1,:), 'b-', 'LineWidth', 1.4);
    plot(offsets, out.meanP(2,:), 'b--', 'LineWidth', 1.4);
    plot(offsets([1 end]), out.no_burst*[1 1], 'k:');
    plot(offsets([1 end]), out.no_control*[1 1], 'k:');
    xlabel('onset of the control relative to the onset of the burst (min)');
    ylabel('mean participation over 60 min'); ylim([0 0.6]);
    legend({'c = 1', 'c = 0.5'}, 'Location', 'northeast');
end
end

function m = mean_part(on, c, burst, S0, Sb, tb, Lb, Lc, T, x0, C, P)
base = @(tt) S0 + (Sb - S0)*(burst && tt >= tb && tt < tb + Lb);
ufun = @(tt) [base(tt)*(1 - c*(tt >= on && tt < on + Lc)), 0, P.Ast];
b = [0, T, tb, tb + Lb, on, on + Lc];
b = unique(b(b >= 0 & b <= T));
[t, X] = simulate_pw(ufun, b, x0, C, P, 0.02);
m = trapz(t, X(:,3)) / T;
end
