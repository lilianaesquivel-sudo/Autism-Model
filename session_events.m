function ev = session_events(variant)
%SESSION_EVENTS  Events of the illustrative example (Table 6, Appendix A)
%   arranged in the 40-min session of Section 7.2.
%   variant = 'slow' (default): swing delivered slowly (vestibular sigma=-1/2)
%           = 'fast'          : swing delivered fast   (vestibular sigma=+1/2)
%           = 'wall'          : slow swing facing a plain wall (visual s=0.1)
%   Domains: 1 auditory, 2 visual, 3 tactile, 4 gustatory/olfactory,
%            5 vestibular/proprioceptive.
if nargin < 1, variant = 'slow'; end
sv = -0.5; svis = 0.4;
switch lower(variant)
    case 'fast', sv = 0.5;
    case 'wall', svis = 0.1;
end
%                name                k  ton  toff  s    rhoH   rhoZ sigma
ev(1) = make_event('musical toy',     1, 2,   8,    0.5, 0.125, 1,   0.5);
ev(2) = make_event('crying',          1, 4,   10,   0.3, 0.75,  0,   0.5);
ev(3) = make_event('door slam',       1, 6,   6.05, 0.8, 0.75,  0,   1);
ev(4) = make_event('swing vest/prop', 5, 15,  25,   0.7, 0.375, 0.75, sv);
ev(5) = make_event('swing visual',    2, 15,  25,   svis, 0.5,  0.5, 0.5);
end
