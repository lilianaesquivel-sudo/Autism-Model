function rng = bistable_range(C, P, Smax, nS)
%BISTABLE_RANGE  Range of constant aversive input SH (with SZ = 0) for
%   which F has three zeros (regulated, saddle, dysregulated), i.e. the
%   bistable range of Proposition 6.3 / Corollary 6.6. The ends are
%   refined by bisection. Returns [] if there is no bistable range in
%   [0, Smax]. Assumes the range is an interval (as in Prop. 6.3).
if nargin < 3, Smax = 1; end
if nargin < 4, nS = 201; end
Ag = linspace(0, 1, 2001);
S  = linspace(0, Smax, nS);
n  = count_zeros(S, C, P, Ag);
i3 = find(n >= 3);
if isempty(i3), rng = []; return; end
bis = @(lo, hi, inside_hi) bisect_end(lo, hi, inside_hi, C, P, Ag);
if i3(1) > 1
    lo = bis(S(i3(1)-1), S(i3(1)), true);
else
    lo = 0;
end
if i3(end) < nS
    hi = bis(S(i3(end)), S(i3(end)+1), false);
else
    hi = Smax;
end
rng = [lo, hi];
end

function x = bisect_end(a, b, inside_hi, C, P, Ag)
% inside_hi = true: a outside, b inside; false: a inside, b outside
for it = 1:30
    m = (a + b) / 2;
    in = count_zeros(m, C, P, Ag) >= 3;
    if in == inside_hi, b = m; else, a = m; end
end
x = (a + b) / 2;
end
