function [t, X, ev] = sim_session(childName, variant, P, T)
%SIM_SESSION  The 40-min session of Section 7.2 (Figure 5) for one child,
%   starting at the resting state with Part(0) = 0.
if nargin < 2, variant = 'slow'; end
if nargin < 3 || isempty(P), P = params_default(); end
if nargin < 4, T = 40; end
C  = child_profile(childName);
ev = session_events(variant);
x0 = resting_state(C, P);
ufun = @(tt) inputs_from_events(tt, ev, C, P);
[t, X] = simulate_pw(ufun, event_breaks(ev, T), x0, C, P, 0.005);
end
