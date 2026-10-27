function dx = model_rhs(~, x, u, C, P)
%MODEL_RHS  Right-hand side of the dynamical model (2.1), time in min.
%   x = [A; R; Part], u = [SH, SZ, Asigma] (constant on the interval).
A = x(1); R = x(2); Pa = x(3);
SH = u(1); SZ = u(2); As = u(3);
G  = gamma_band(A, P);
dA = P.alphaH*SH*(1 - A) + P.alphaZ*SZ*(As - A) ...
   + P.eta*C.Zbar*(1 - G)*(P.Ast - A) + P.beta*R*(C.Ab - A);
dR = P.mu*(1 - R) - P.nuH*SH*R - P.nuA*max(A - P.Ad, 0)*R;
dP = P.lambda*R*G*(1 - Pa) - P.delta*(1 - G)*Pa;
dx = [dA; dR; dP];
end
