function [t, X] = simulate_pw(ufun, breaks, x0, C, P, dt)
%SIMULATE_PW  Caratheodory solution of (2.1) for piecewise-constant
%   inputs. ufun(t) returns u = [SH, SZ, Asigma]; breaks contains 0, T
%   and all switching times. The system is integrated separately on each
%   interval, the final state of one interval being the initial state of
%   the next (Section 7). dt is the output step (min).
if nargin < 6, dt = 0.01; end
breaks = unique(breaks(:))';
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
t = []; X = []; x = x0(:);
for i = 1:numel(breaks) - 1
    a = breaks(i); b = breaks(i+1);
    if b - a <= 0, continue; end
    u  = ufun((a + b) / 2);
    tt = linspace(a, b, max(3, ceil((b - a) / dt) + 1));
    [ti, Xi] = ode45(@(s, y) model_rhs(s, y, u, C, P), tt, x, opts);
    if isempty(t)
        t = ti(:); X = Xi;
    else
        t = [t; ti(2:end)]; X = [X; Xi(2:end, :)];
    end
    x = Xi(end, :)';
end
end
