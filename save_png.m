function save_png(fname)
%SAVE_PNG  Save the current figure as PNG (MATLAB and Octave).
try
    print(gcf, fname, '-dpng', '-r150');
catch
    try
        print(gcf, fname, '-dpngcairo');     % Octave + gnuplot without Ghostscript
    catch err
        warning('Could not save %s: %s', fname, err.message);
    end
end
end
