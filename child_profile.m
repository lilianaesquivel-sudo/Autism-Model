function C = child_profile(name)
%CHILD_PROFILE  Sensory profiles of Table 4 and resting levels of Table 5.
%   C = CHILD_PROFILE('A') or CHILD_PROFILE('B') returns a struct with the
%   domain scores H, L, Z (domains: auditory, visual, tactile,
%   gustatory/olfactory, vestibular/proprioceptive), the overall seeking
%   score Zbar (3.3), and the resting level Ab.
switch upper(name)
    case 'A'
        C.H  = [0.8 0.3 0.6 0.4 0.7];
        C.L  = [0.1 0.1 0.1 0.1 0.1];
        C.Z  = [0.2 0.1 0.2 0.1 0.1];
        C.Ab = 0.5;
    case 'B'
        C.H  = [0.2 0.8 0.1 0.2 0.1];
        C.L  = [0.8 0.1 0.5 0.4 0.6];
        C.Z  = [0.8 0.2 0.6 0.5 0.9];
        C.Ab = 0.25;
    otherwise
        error('child_profile: unknown child %s', name);
end
C.name = upper(name);
C.Zbar = mean(C.Z);
end
