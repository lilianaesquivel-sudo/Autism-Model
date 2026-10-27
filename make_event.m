function e = make_event(name, k, ton, toff, s, rhoH, rhoZ, sigma)
%MAKE_EVENT  Event (incidental or therapeutic component) of constant
%   intensity s on [ton, toff) in domain k, with signature
%   (rhoH, rhoZ, sigma), as in (3.2), (3.5) and (4.1).
e = struct('name', name, 'k', k, 'ton', ton, 'toff', toff, 's', s, ...
           'rhoH', rhoH, 'rhoZ', rhoZ, 'sigma', sigma);
end
