function x0 = resting_state(C, P)
%RESTING_STATE  Equilibrium without input (Section 2.3), with Part = 0.
E  = equilibria_const([0, 0, P.Ast], C, P);
x0 = [E(1).A; E(1).R; 0];
end
