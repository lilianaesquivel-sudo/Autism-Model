function out = fig_bifurcation(doplot, P)
%FIG_BIFURCATION  Section 7.3 / Figure 6: equilibria versus a constant
%   aversive input (SZ = 0), bistable ranges for eta = 2.5 and eta = 0,
%   the hysteresis loop of Child A, and the participation of Child B
%   along the regulated branch.
if nargin < 1, doplot = true; end
if nargin < 2, P = params_default(); end
names = {'A', 'B'};  col = {'b', 'r'};
P0 = P; P0.eta = 0;
for c = 1:2
    C = child_profile(names{c});
    out.(['range_' names{c}])      = bistable_range(C, P);
    out.(['range_eta0_' names{c}]) = bistable_range(C, P0);
end
% equilibrium branches
S = linspace(0, 1, 401);
if doplot, figure('Name', 'Figure 6', 'Position', [100 100 900 430]); end
for c = 1:2
    C = child_profile(names{c});
    Sr = []; Ar = []; Pr = []; Ss = []; As = []; Ps = []; Su = []; Au = []; Pu = [];
    for i = 1:numel(S)
        E = equilibria_const([S(i), 0, P.Ast], C, P, 4001);
        for e = E
            if e.stable && e.A < P.Ad
                Sr(end+1) = S(i); Ar(end+1) = e.A; Pr(end+1) = e.Part; %#ok<AGROW>
            elseif e.stable
                Ss(end+1) = S(i); As(end+1) = e.A; Ps(end+1) = e.Part; %#ok<AGROW>
            else
                Su(end+1) = S(i); Au(end+1) = e.A; Pu(end+1) = e.Part; %#ok<AGROW>
            end
        end
    end
    if c == 2
        [~, im] = max(Pr);  out.B_Ppeak_SH = Sr(im);  out.B_Ppeak = Pr(im);
    end
    if doplot
        subplot(1,2,1); hold on; box on;
        plot(Sr, Ar, [col{c} '.'], Ss, As, [col{c} '.'], 'MarkerSize', 4);
        plot(Su, Au, [col{c} '.'], 'MarkerSize', 1);
        subplot(1,2,2); hold on; box on;
        plot(Sr, Pr, [col{c} '.'], Ss, Ps, [col{c} '.'], 'MarkerSize', 4);
        plot(Su, Pu, [col{c} '.'], 'MarkerSize', 1);
    end
end
% hysteresis of Child A: slow ramp 0 -> 0.8 -> 0, 400 min each way
C = child_profile('A');
Tr = 400; Smx = 0.8;
x0 = resting_state(C, P);
rhs = @(t, x) model_rhs(t, x, [Smx*min(t, 2*Tr - t)/Tr, 0, P.Ast], C, P);
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[t1, X1] = ode45(rhs, linspace(0, Tr, 40001), x0, opts);
[t2, X2] = ode45(rhs, linspace(Tr, 2*Tr, 40001), X1(end,:)', opts);
iu = find(X1(:,1) > P.Ad, 1);     out.hyst_up   = Smx*t1(iu)/Tr;
id = find(X2(:,1) < P.Ad, 1);     out.hyst_down = Smx*(2*Tr - t2(id))/Tr;
if doplot
    subplot(1,2,1);
    plot(Smx*t1/Tr, X1(:,1), 'k-', Smx*(2*Tr - t2)/Tr, X2(:,1), 'k--');
    plot([0 1], [P.Ad P.Ad], 'k:');
    xlabel('aversive input S^H'); ylabel('equilibrium arousal A^e'); title('(a)');
    ylim([0 1]);
    subplot(1,2,2);
    xlabel('aversive input S^H'); ylabel('equilibrium participation P^e'); title('(b)');
    ylim([0 1]);
end
end
