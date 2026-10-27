function out = eigen_discussion(P)
%EIGEN_DISCUSSION  Eigenvalues quoted in the Discussion (Remark 6.x on
%   critical slowing down): slowest eigenvalue of the dysregulated state
%   of Child A at SH = 0.2 and SH = 0.155, and eigenvalues of the
%   regulated state near the upper end of the bistable range.
if nargin < 1, P = params_default(); end
C  = child_profile('A');
rg = bistable_range(C, P);
for S = [0.2, 0.155]
    E = equilibria_const([S, 0, P.Ast], C, P);
    E = E([E.stable]);
    ev = E(end).eig;
    out.(sprintf('dys_slowest_%03d', round(1000*S))) = max(real(ev));
end
S = rg(2) - 1e-4;
E = equilibria_const([S, 0, P.Ast], C, P);
out.reg_eig_upper_end = sort(real(E(1).eig))';
out.SH_upper_end = rg(2);
out.capacity_eig = -(P.mu + P.nuH*S);
end
