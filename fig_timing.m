function out = fig_timing(doplot, P)
%FIG_TIMING  Section 7.5 / Figure 8: Child A under an ambient aversive
%   input SH = 0.22 with a burst SH = 0.8 from minute 5, followed or not
%   by a quiet break (SH = 0).
if nargin < 1, doplot = true; end
if nargin < 2, P = params_default(); end
S0 = 0.22; Sb = 0.8; tb = 5; Tlong = 600;
[C, xreg, xdys] = timing_setup(P, S0);
out.reg = xreg';  out.dys = xdys';
dysreg = @(seg) final_A(seg, xreg, C, P) > P.Ad;

% shortest burst that dysregulates
lo = 0; hi = 3;
for it = 1:30
    m = (lo + hi)/2;
    if dysreg([m Sb; Tlong S0]), hi = m; else, lo = m; end
end
out.burst_threshold = hi;

% shortest quiet break after a 3-min burst
lo = 0; hi = 30;
for it = 1:30
    m = (lo + hi)/2;
    if dysreg([3 Sb; m 0; Tlong S0]), lo = m; else, hi = m; end
end
out.min_break = hi;

% capacity at the end of the break, and recovery times
[~, Xb] = sim_profile([3 Sb], xreg, C, P);  Rb = Xb(end,2);
out.R_end_burst = Rb;
brk = [out.min_break, 5, 10];
for i = 1:3
    [t, X] = sim_profile([3 Sb; brk(i) 0; 120 S0], xreg, C, P);
    [~, j] = min(abs(t - 3 - brk(i)));
    out.R_after_break(i) = X(j,2);
    k = find(t > 3 + brk(i) & X(:,2) >= 0.9*xreg(2), 1);
    out.t90_R(i) = t(k) - 3;                 % minutes after the burst ends
    k = find(t > 3 + brk(i) & X(:,3) >= 0.9*xreg(3), 1);
    out.t90_P(i) = t(k) - 3;
end
% bound of Proposition 6.11 (time for R to reach 0.9 Rreg)
out.bound_t90 = log((1 - Rb) / (1 - 0.9*xreg(2))) / P.mu;

if doplot
    figure('Name', 'Figure 8', 'Position', [100 100 700 600]);
    cases = {[tb S0; 0.5 Sb; 70 S0], [tb S0; 3 Sb; 70 S0], ...
             [tb S0; 3 Sb; 5 0; 70 S0], [tb S0; 3 Sb; 10 0; 70 S0]};
    sty = {'g-', 'r-', 'b-', 'm-'};
    lab = {'arousal A', 'capacity R', 'participation P'};
    for i = 1:4
        [t, X] = sim_profile(cases{i}, xreg, C, P);
        for k = 1:3
            subplot(3,1,k); hold on; box on;
            plot(t, X(:,k), sty{i}, 'LineWidth', 1.2);
        end
    end
    tt = linspace(tb + 3, 70, 200);
    subplot(3,1,2); plot(tt, 1 - (1 - Rb)*exp(-P.mu*(tt - tb - 3)), 'k:');
    plot([0 70], 0.9*xreg(2)*[1 1], 'k--');
    subplot(3,1,1); plot([0 70], [P.Ad P.Ad], 'k:');
    legend({'burst 0.5 min', 'burst 3 min, no break', ...
            'burst 3 min, break 5 min', 'burst 3 min, break 10 min'}, ...
            'Location', 'southeast');
    for k = 1:3, subplot(3,1,k); ylim([0 1]); xlim([0 70]); ylabel(lab{k}); end
    xlabel('time (min)');
end
end

function A = final_A(seg, x0, C, P)
[~, X] = sim_profile(seg, x0, C, P, 1);
A = X(end,1);
end
