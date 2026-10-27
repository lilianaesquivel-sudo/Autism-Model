function R = capacity_nullcline(A, SH, P)
%CAPACITY_NULLCLINE  Equilibrium capacity (6.4) in dimensional form.
R = P.mu ./ (P.mu + P.nuH*SH + P.nuA*max(A - P.Ad, 0));
end
