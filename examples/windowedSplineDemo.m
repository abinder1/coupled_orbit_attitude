% Copyright (c) 2026 Andrew Binder

% windowedSplineDemo
%   Demonstrate a windowed quintic B-spline least-squares fit to position
%   and velocity. Each data column is fitted independently with the same
%   time basis. The fit is approximate, not exact osculatory interpolation.
%
%   Simple interior knots give C^4 continuity, including at window
%   boundaries. The batch fit uses all node times to construct the knots,
%   then freezes coefficients window by window with lookahead. It is not
%   an online append API; see coa.interp.fitWindowedSpline and
%   coa.interp.evaluateSpline.
%
%   The synthetic truth is used only for diagnostics: jerk comparisons,
%   component-wise error histories, RMS and peak errors, knot jumps, and
%   fit/evaluation time.
%   In real data, higher-derivative error requires an independent reference.
%   No toolboxes required.

clear;
clc;
close('all');
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'toolbox'));

%% Configuration and synthetic data
opts.order = 6;              % Quintic, C^4 at simple interior knots
opts.window_size = 40;       % Window size in nodes
opts.stride = 20;            % Window step in nodes
opts.velocity_weight = 1.0;  % Weight before time-step scaling
n = 500;
n_queries = 20000;

% Replace these with a strictly increasing time vector and N-by-D arrays.
% Velocity must use position units per unit of time.
rng(1);
t = cumsum(0.05 + 0.03 * rand(n, 1));
sample_truth = truth(t);
position = sample_truth(:, :, 1);
velocity = sample_truth(:, :, 2);

%% Fit and evaluate
timer = tic;
sp = coa.interp.fitWindowedSpline(t, position, velocity, opts);
fit_seconds = toc(timer);
tq = linspace(t(1), t(end), n_queries).';
timer = tic;
estimated = coa.interp.evaluateSpline(sp, tq, 3);
eval_seconds = toc(timer);
reference = truth(tq);

% Evaluate both polynomial pieces at every interior knot, not only the
% window boundaries. Exclude endpoints, where only one piece exists.
knots = t(2:end - 1);
left = coa.interp.evaluateSpline(sp, knots, 3, 'left');
right = coa.interp.evaluateSpline(sp, knots, 3, 'right');
errors = estimated - reference;
jumps = right - left;
names = {'Position', 'Velocity', 'Acceleration', 'Jerk'};
component_names = arrayfun(@(i) sprintf('Component %d', i), ...
                           1:size(position, 2), 'UniformOutput', false);

%% Report per-component performance
fprintf('Fit: %.3f s, %d windows; evaluation: %.3f s (%.3g us/query).\n', ...
        fit_seconds, numel(sp.window_starts), eval_seconds, ...
        1e6 * eval_seconds / n_queries);
fprintf('Errors against analytic truth on %d evaluation times:\n', n_queries);
for m = 1:numel(names)
    rms_error = sqrt(mean(errors(:, :, m).^2, 1));
    peak_error = max(abs(errors(:, :, m)), [], 1);
    max_jump = max(abs(jumps(:, :, m)), [], 1);
    for component = 1:size(position, 2)
        fprintf('  %-12s component %d: RMS %.3e, max %.3e, knot jump %.3e\n', ...
                names{m}, component, rms_error(component), ...
                peak_error(component), max_jump(component));
    end
end

%% Jerk comparison and window boundaries
figure('Name', 'Spline jerk comparison');
tiledlayout(size(position, 2), 1);
window_times = t(sp.window_starts(2:end));
for component = 1:size(position, 2)
    nexttile;
    plot(tq, reference(:, component, 4), 'k-', ...
         tq, estimated(:, component, 4), 'r--');
    hold on;
    for boundary = 1:numel(window_times)
        xline(window_times(boundary), ':', 'Color', [0.6, 0.6, 0.6], ...
              'HandleVisibility', 'off');
    end
    hold off;
    title(component_names{component});
    xlabel('Time');
    ylabel('Jerk');
    legend('Truth', 'Spline', 'Location', 'best');
    grid on;
end

%% Error histories
figure('Name', 'Spline error histories');
tiledlayout(2, 2);
for m = 1:numel(names)
    nexttile;
    semilogy(tq, max(abs(errors(:, :, m)), realmin));
    title(names{m});
    xlabel('Time');
    ylabel('Absolute error');
    grid on;
    legend(component_names, 'Location', 'best');
end

%% RMS and peak error by component
figure('Name', 'Spline error summary');
tiledlayout(2, 2);
for m = 1:numel(names)
    nexttile;
    rms_error = sqrt(mean(errors(:, :, m).^2, 1));
    peak_error = max(abs(errors(:, :, m)), [], 1);
    bar(max([rms_error(:), peak_error(:)], realmin));
    set(gca, 'YScale', 'log', 'XTick', 1:numel(component_names), ...
        'XTickLabel', component_names);
    title(names{m});
    ylabel('Absolute error');
    legend('RMS', 'Maximum', 'Location', 'best');
    grid on;
end

%% Continuity at every interior knot
figure('Name', 'Spline knot continuity');
tiledlayout(2, 2);
for m = 1:numel(names)
    nexttile;
    semilogy(knots, max(abs(jumps(:, :, m)), realmin));
    title(names{m});
    xlabel('Time');
    ylabel('Absolute left-right jump');
    grid on;
    legend(component_names, 'Location', 'best');
end

%% Analytic trajectory and its first three derivatives
function x = truth(t)
    % Synthetic three-component trajectory used to validate the fit.
    t = t(:);
    x = zeros(numel(t), 3, 4);
    x(:, 1, 1) = cos(t) + 0.3 * sin(3 * t);
    x(:, 1, 2) = -sin(t) + 0.9 * cos(3 * t);
    x(:, 1, 3) = -cos(t) - 2.7 * sin(3 * t);
    x(:, 1, 4) = sin(t) - 8.1 * cos(3 * t);
    x(:, 2, 1) = sin(1.3 * t);
    x(:, 2, 2) = 1.3 * cos(1.3 * t);
    x(:, 2, 3) = -1.69 * sin(1.3 * t);
    x(:, 2, 4) = -2.197 * cos(1.3 * t);
    x(:, 3, 1) = 0.2 * cos(0.7 * t);
    x(:, 3, 2) = -0.14 * sin(0.7 * t);
    x(:, 3, 3) = -0.098 * cos(0.7 * t);
    x(:, 3, 4) = 0.0686 * sin(0.7 * t);
end
