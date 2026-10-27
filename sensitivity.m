function out = sensitivity(nsamp, seed, verbose)
%SENSITIVITY  Sensitivity of the bistable range (Section 7.1).
%   (1) One-at-a-time: beta, alphaH, eta, mu, nuH, nuA multiplied by
%       0.5, 0.75, 1.5, 2.
%   (2) 1/mu = 10 and 60 min with nuH, nuA fixed (Child A).
%   (3) Joint: each of the six rates multiplied by a factor drawn
%       log-uniformly from [0.5, 2], and kappa ~ U[0.10, 0.20] (with
%       Ad = A* + kappa); nsamp samples (default 2000).
%   The joint results depend on the random samples; the percentages in
%   the paper (100% and 91%) are reproduced up to sampling variation.
if nargin < 1, nsamp = 2000; end
if nargin < 2, seed = 1; end
if nargin < 3, verbose = true; end
P0 = params_default();
rates = {'beta', 'alphaH', 'eta', 'mu', 'nuH', 'nuA'};
fac = [0.5 0.75 1.5 2];
names = {'A', 'B'};
out.oat_all_bistable = true;
for c = 1:2
    C = child_profile(names{c});
    for r = 1:numel(rates)
        for f = fac
            P = P0; P.(rates{r}) = P0.(rates{r})*f;
            rg = bistable_range(C, P, 1);
            out.oat.(names{c}).(rates{r})(fac == f, :) = pad(rg);
            if isempty(rg), out.oat_all_bistable = false; end
        end
    end
end
CA = child_profile('A');
P = P0; P.mu = 1/10;  out.range_A_mu10 = bistable_range(CA, P);
P = P0; P.mu = 1/60;  out.range_A_mu60 = bistable_range(CA, P);

% joint sampling
rand('seed', seed); %#ok<RAND>  works in MATLAB and Octave
U = rand(nsamp, numel(rates) + 1);
for c = 1:2
    C = child_profile(names{c});
    W = nan(nsamp, 2); above1 = false(nsamp, 1);
    for i = 1:nsamp
        P = P0;
        for r = 1:numel(rates)
            P.(rates{r}) = P0.(rates{r}) * 2^(2*U(i,r) - 1);
        end
        P.kappa = 0.10 + 0.10*U(i, end);
        P.Ad    = P.Ast + P.kappa;
        rg = bistable_range(C, P, 1);
        if ~isempty(rg)
            W(i,:) = rg;
        else
            above1(i) = ~isempty(bistable_range(C, P, 3, 301));
        end
    end
    ok = ~isnan(W(:,1));
    out.joint.(names{c}).fraction_bistable = mean(ok);
    out.joint.(names{c}).median_width      = median(W(ok,2) - W(ok,1));
    out.joint.(names{c}).only_above_1      = sum(above1);
    out.joint.(names{c}).neither           = sum(~ok & ~above1);
end
if verbose
    fprintf('One-at-a-time (x0.5..x2): bistable in all cases: %d\n', out.oat_all_bistable);
    fprintf('Child A, 1/mu = 10: (%.2f, %.2f);  1/mu = 60: (%.2f, %.2f)\n', ...
        out.range_A_mu10, out.range_A_mu60);
    for c = 1:2
        j = out.joint.(names{c});
        fprintf('Joint, Child %s: bistable %.1f%%, median width %.2f, only for SH>1: %d, none: %d\n', ...
            names{c}, 100*j.fraction_bistable, j.median_width, j.only_above_1, j.neither);
    end
end
end

function v = pad(rg)
if isempty(rg), v = [NaN NaN]; else, v = rg; end
end
