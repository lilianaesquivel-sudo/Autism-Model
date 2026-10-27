function E = equilibria_const(u, C, P, nA)
%EQUILIBRIA_CONST  Equilibria of (2.1) for constant inputs u.
%   Equilibria are the zeros of F(A) = f(A, R(A)) (Section 6.4). Returns
%   a struct array with fields A, R, Part (6.5), eig (eigenvalues of the
%   (A,R) Jacobian, 1/min) and stable (true/false).
if nargin < 4, nA = 20001; end
SH = u(1); SZ = u(2); As = u(3);
Fh = @(A) P.alphaH*SH*(1 - A) + P.alphaZ*SZ*(As - A) ...
        + P.eta*C.Zbar*(1 - gamma_band(A, P)).*(P.Ast - A) ...
        + P.beta*capacity_nullcline(A, SH, P).*(C.Ab - A);
Ag = linspace(0, 1, nA);
Fg = Fh(Ag);
idx = find(Fg(1:end-1) .* Fg(2:end) < 0 | Fg(1:end-1) == 0);
E = struct('A', {}, 'R', {}, 'Part', {}, 'eig', {}, 'stable', {});
for i = idx
    if Fg(i) == 0
        Ae = Ag(i);
    else
        Ae = fzero(Fh, [Ag(i), Ag(i+1)]);
    end
    Re = capacity_nullcline(Ae, SH, P);
    G  = gamma_band(Ae, P);
    Pe = P.lambda*Re*G / (P.lambda*Re*G + P.delta*(1 - G));
    hp = -(1 - G) - 2*(Ae - P.Ast)^2 / P.kappa^2 * G;     % h'(A)
    J  = [-P.alphaH*SH - P.alphaZ*SZ + P.eta*C.Zbar*hp - P.beta*Re, ...
           P.beta*(C.Ab - Ae);
          -P.nuA*(Ae > P.Ad)*Re, ...
          -P.mu - P.nuH*SH - P.nuA*max(Ae - P.Ad, 0)];
    ev = eig(J);
    E(end+1) = struct('A', Ae, 'R', Re, 'Part', Pe, 'eig', ev, ...
                      'stable', all(real(ev) < 0)); %#ok<AGROW>
end
end
