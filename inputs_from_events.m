function u = inputs_from_events(t, events, C, P)
%INPUTS_FROM_EVENTS  Session inputs u = [SH, SZ, Asigma] at time t,
%   obtained from all active events by the aggregation rule (5.1), the
%   modulatory shares (5.3) and the target level (5.4).
pH = 1; pZ = 1; W = 0; wp = 0; wm = 0;
for j = 1:numel(events)
    e = events(j);
    if t >= e.ton && t < e.toff
        [EH, EZ] = event_contributions(e, e.s, C, P);
        pH = pH * (1 - EH);
        pZ = pZ * (1 - EZ);
        W  = W + EZ;
        wp = wp + EZ * max(e.sigma, 0);
        wm = wm + EZ * max(-e.sigma, 0);
    end
end
SH = 1 - pH;
SZ = 1 - pZ;
if W > 0
    wp = wp / W;  wm = wm / W;
else
    wp = 0; wm = 0;
end
As = P.Ast + wp * (1 - P.Ast) - wm * P.Ast;
u = [SH, SZ, As];
end
