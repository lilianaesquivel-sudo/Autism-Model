function [t, X] = sim_profile(seg, x0, C, P, dt)
%SIM_PROFILE  Simulate (2.1) with SZ = 0 and an aversive input that is
%   constant on consecutive segments. seg is an n x 2 array whose rows
%   are [duration, SH]. Time starts at 0.
if nargin < 5, dt = 0.01; end
b = [0, cumsum(seg(:,1))'];
ufun = @(tt) [seg(min(find(tt < b(2:end), 1), size(seg,1)), 2), 0, P.Ast];
[t, X] = simulate_pw(ufun, b, x0, C, P, dt);
end
