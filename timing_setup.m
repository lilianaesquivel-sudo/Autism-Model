function [C, xreg, xdys] = timing_setup(P, S0)
%TIMING_SETUP  Child A under a constant ambient aversive input S0
%   (Section 7.5): regulated and dysregulated equilibria [A; R; Part].
C = child_profile('A');
E = equilibria_const([S0, 0, P.Ast], C, P);
E = E([E.stable]);
xreg = [E(1).A;   E(1).R;   E(1).Part];
xdys = [E(end).A; E(end).R; E(end).Part];
end
