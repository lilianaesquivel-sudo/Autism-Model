function [EH, EZ] = event_contributions(e, s, C, P)
%EVENT_CONTRIBUTIONS  Individualized contributions (3.7)-(3.8), (4.2)
%   of event e at intensity s for child C.
k  = e.k;
r  = 1 - P.gamma * C.L(k) * (1 - s);          % registration factor (3.7)
Ht = P.h0 + (1 - P.h0) * C.H(k);              % baseline susceptibility
EH = r * e.rhoH * Ht * s;                     % aversive contribution
EZ = e.rhoZ * C.Z(k) * s;                     % seeking contribution
end
