function out = fig_session(doplot, P)
%FIG_SESSION  Section 7.2 / Figure 5: two children in the same session,
%   with the swing delivered slowly (calming) or fast (alerting).
if nargin < 1, doplot = true; end
if nargin < 2, P = params_default(); end
names = {'A', 'B'};
if doplot, fh = figure('Name', 'Figure 5', 'Position', [100 100 900 650]); end
for c = 1:2
    C = child_profile(names{c});
    [t, X]  = sim_session(names{c}, 'slow', P);
    [t2, X2] = sim_session(names{c}, 'fast', P);
    aud = t >= 2 & t < 10;  sw = t >= 15 & t < 25;  sw2 = t2 >= 15 & t2 < 25;
    s.rest      = X(1,1);
    s.aud_Amax  = max(X(aud,1));
    s.aud_Rmin  = min(X(aud,2));
    s.aud_Pmax  = max(X(aud,3));
    [~, i15]    = min(abs(t - 15));
    s.R_at15    = X(i15,2);
    s.slow_Pmean = trapz(t(sw), X(sw,3)) / 10;
    s.slow_Amean = trapz(t(sw), X(sw,1)) / 10;
    s.fast_Pmean = trapz(t2(sw2), X2(sw2,3)) / 10;
    s.fast_Amean = trapz(t2(sw2), X2(sw2,1)) / 10;
    out.(names{c}) = s;
    if doplot
        lab = {'arousal A', 'capacity R', 'participation P'};
        for k = 1:3
            subplot(3, 2, 2*(k-1) + c); hold on; box on;
            plot(t, X(:,k), 'b-', 'LineWidth', 1.3);
            plot(t2, X2(:,k), 'r--', 'LineWidth', 1.3);
            if k == 1
                plot([0 40], [P.Ad P.Ad], 'k:');
                plot([0 40], [P.Ast P.Ast], 'k:');
                title(['Child ' names{c}]);
            end
            ylim([0 1]); xlim([0 40]); ylabel(lab{k});
            if k == 3, xlabel('time (min)'); end
            if k == 1 && c == 1, legend({'slow swing (calming)', 'fast swing (alerting)'}, 'Location', 'southwest'); end
        end
    end
end
end
