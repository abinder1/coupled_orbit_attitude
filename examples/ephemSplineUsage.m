% Copyright (c) 2026 Andrew Binder
function [query_times, values] = ephemSplineUsage()
    % Fit synthetic samples with ephem_integrator and evaluate derivatives.
    % values(:, :, 1:4) contains position, velocity, acceleration, and jerk.
    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(root);
    setupProjectPaths();

    t = linspace(0, 4, 101).';
    position = [sin(t), cos(t)];
    velocity = [cos(t), -sin(t)];
    opts = struct('order', 6, 'window_size', 40, 'stride', 20, ...
                  'velocity_weight', 1);
    sp = emframe.interp.fitWindowedSpline(t, position, velocity, opts);
    query_times = linspace(t(1), t(end), 201).';
    values = emframe.interp.evaluateSpline(sp, query_times, 3);
end
