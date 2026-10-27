function b = event_breaks(ev, T)
%EVENT_BREAKS  Switching times of a list of events on [0, T].
b = [0, T, [ev.ton], [ev.toff]];
b = unique(b(b >= 0 & b <= T));
end
