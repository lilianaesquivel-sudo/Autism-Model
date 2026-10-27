function [n, Az] = count_zeros(SH, C, P, Agrid)
%COUNT_ZEROS  Number of zeros of F (Section 6.4) for constant aversive
%   input SH (vector allowed) and no seeking-related input. Returns the
%   number of sign changes of F on Agrid for each SH, and, for scalar SH,
%   the approximate locations.
if nargin < 4, Agrid = linspace(0, 1, 4001); end
A  = Agrid(:);
SH = SH(:)';
G  = gamma_band(A, P);
R  = P.mu ./ (P.mu + P.nuH*SH + P.nuA*max(A - P.Ad, 0));     % nA x nS
F  = P.alphaH*(1 - A)*SH + P.eta*C.Zbar*(1 - G).*(P.Ast - A) ...
   + P.beta*R.*(C.Ab - A);
s  = sign(F);
ch = s(1:end-1, :) .* s(2:end, :) < 0;
n  = sum(ch, 1);
if nargout > 1 && numel(SH) == 1
    Az = A(find(ch));
end
end
