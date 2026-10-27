function out = prop63_explicit(SH, SZ, As, C, P)
%PROP63_EXPLICIT  Closed-form quantities of Proposition 6.3 (eta = 0)
%   for constant inputs. Returns Abar, Areg, p, q, m and the flags
%   bistable (Areg < Ad, q > p and 2q - p < m) and monostable_iii
%   (Abar <= Ad).
aH  = P.alphaH*SH/P.beta;  aZ = P.alphaZ*SZ/P.beta;  a = aH + aZ;
thH = P.nuH*SH/P.mu;       thA = P.nuA/P.mu;
out.Abar = (aH + aZ*As) / a;
c        = 1 + thH;
out.Areg = (a*c*out.Abar + C.Ab) / (a*c + 1);
out.p    = (1 + a*c) / (a*thA);
out.q    = sqrt(max(out.Abar - C.Ab, 0) / (a*thA));
out.m    = out.Abar - P.Ad;
out.bistable       = out.Areg < P.Ad && out.q > out.p && ...
                     2*out.q - out.p < out.m;
out.monostable_iii = out.Abar <= P.Ad;
end
