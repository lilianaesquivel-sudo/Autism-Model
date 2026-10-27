function G = gamma_band(A, P)
%GAMMA_BAND  Function Gamma of (2.2): proximity of arousal to A_*.
G = exp(-((A - P.Ast) / P.kappa).^2);
end
