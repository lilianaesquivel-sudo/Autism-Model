function out = example_inputs(verbose)
%EXAMPLE_INPUTS  Illustrative example of Section 5.1 (Table 7, Figure 4):
%   individualized contributions and total inputs at t0 (three auditory
%   events) and t1 (swing), for Child A and Child B.
if nargin < 1, verbose = true; end
P = params_default();
t0 = 6.01;  t1 = 20;               % instants inside the session of Sec. 7.2
names = {'A', 'B'};
for c = 1:2
    C = child_profile(names{c});
    for v = {'slow', 'fast', 'wall'}
        ev = session_events(v{1});
        for tt = [t0, t1]
            u = inputs_from_events(tt, ev, C, P);
            key = sprintf('%s_%s_t%d', names{c}, v{1}, tt == t1);
            out.(key) = u;
        end
    end
    % individual contributions at t0 and t1 (slow swing)
    ev = session_events('slow');
    E = zeros(numel(ev), 2);
    for j = 1:numel(ev)
        [E(j,1), E(j,2)] = event_contributions(ev(j), ev(j).s, C, P);
    end
    out.(['E_' names{c}]) = E;
end
% registration: reduction of the aversive contribution for Child B
CB = child_profile('B');
ev = session_events('slow');
out.reg_crying = P.gamma*CB.L(1)*(1 - ev(2).s);
out.reg_slam   = P.gamma*CB.L(1)*(1 - ev(3).s);
EA = out.E_A;
out.A_sum_t0 = sum(EA(1:3, 1));
% shares of the seeking-related input of Child B at t1
EB = out.E_B;  W = EB(4,2) + EB(5,2);
out.B_t1_calming    = EB(4,2)*0.5/W;
out.B_t1_alerting   = EB(5,2)*0.5/W;
out.B_t1_organizing = 1 - out.B_t1_calming - out.B_t1_alerting;

if verbose
    ev = session_events('slow');
    fprintf('\nIndividualized contributions E^H, E^Z (Table 7)\n');
    fprintf('%-18s  %6s %6s   %6s %6s\n', 'event', 'A:E^H', 'A:E^Z', 'B:E^H', 'B:E^Z');
    for j = 1:numel(ev)
        fprintf('%-18s  %6.3f %6.3f   %6.3f %6.3f\n', ev(j).name, ...
            out.E_A(j,1), out.E_A(j,2), out.E_B(j,1), out.E_B(j,2));
    end
    fprintf('\nTotals [SH SZ Asigma]\n');
    f = fieldnames(out);
    for i = 1:numel(f)
        if numel(out.(f{i})) == 3 && ~strncmp(f{i}, 'E_', 2)
            fprintf('  %-14s %6.3f %6.3f %6.3f\n', f{i}, out.(f{i}));
        end
    end
end
end
