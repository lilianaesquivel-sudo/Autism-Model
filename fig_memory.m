function out = fig_memory(doplot, P)
%FIG_MEMORY  Section 7.4 / Figure 7: Child A exposed to six aversive
%   stimuli of 1 min (SH = 0.8), every 3 min or every 15 min.
if nargin < 1, doplot = true; end
if nargin < 2, P = params_default(); end
C  = child_profile('A');
x0 = resting_state(C, P);
gaps = [3, 15];  tag = {'close', 'spaced'};  sty = {'b-', 'r--'};
if doplot, figure('Name', 'Figure 7', 'Position', [100 100 700 500]); end
for g = 1:2
    on = 1 + gaps(g)*(0:5);  off = on + 1;  T = off(end) + 20;
    ufun = @(tt) [0.8*any(tt >= on & tt < off), 0, P.Ast];
    [t, X] = simulate_pw(ufun, unique([0, T, on, off]), x0, C, P, 0.005);
    peaks = zeros(1,6); Rbefore = zeros(1,6); Rafter = zeros(1,6);
    for j = 1:6
        w = t >= on(j) & t <= off(j);
        peaks(j)   = max(X(w,1));
        [~, i0]    = min(abs(t - on(j)));   Rbefore(j) = X(i0,2);
        [~, i1]    = min(abs(t - off(j)));  Rafter(j)  = X(i1,2);
    end
    out.(tag{g}) = struct('peaks', peaks, 'Rbefore', Rbefore, 'Rafter', Rafter);
    if doplot
        lab = {'arousal A', 'capacity R', 'participation P'};
        for k = 1:3
            subplot(3,1,k); hold on; box on;
            plot(t, X(:,k), sty{g}, 'LineWidth', 1.2); ylim([0 1]); ylabel(lab{k});
        end
    end
end
if doplot
    subplot(3,1,1); plot(xlim, [P.Ad P.Ad], 'k:');
    legend({'every 3 min', 'every 15 min'}, 'Location', 'southeast');
    subplot(3,1,3); xlabel('time (min)');
end
end
