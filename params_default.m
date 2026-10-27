function P = params_default()
%PARAMS_DEFAULT  Parameter values of Table 5 (rates in 1/min).
%   P = PARAMS_DEFAULT() returns a struct with the reference levels, the
%   two parameters of the input construction, and the nine rates of the
%   dynamical model (2.1).

% Reference levels
P.Ast   = 0.5;            % optimal arousal level A_*
P.kappa = 0.15;           % width of the optimal region
P.Ad    = P.Ast + P.kappa;% depletion threshold A_d = A_* + kappa

% Input construction
P.h0    = 0.2;            % baseline aversive susceptibility
P.gamma = 0.5;            % maximal attenuation by hyporesponsiveness

% Rates (1/min)
P.beta   = 5;             % self-regulation at full capacity
P.alphaH = 2.5;           % coupling of aversive input to arousal
P.alphaZ = 10;            % coupling of seeking-related input to arousal
P.eta    = 2.5;           % strength of endogenous seeking
P.mu     = 0.05;          % recovery of regulatory capacity
P.nuH    = 0.2;           % depletion by aversive input
P.nuA    = 16;            % depletion by over-arousal
P.lambda = 0.2;           % gain of participation
P.delta  = 0.3;           % loss of participation
end
