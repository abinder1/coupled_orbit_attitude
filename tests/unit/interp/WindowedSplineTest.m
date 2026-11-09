% Copyright (c) 2026 Andrew Binder

classdef WindowedSplineTest < matlab.unittest.TestCase
    % Tests for the windowed position and velocity spline.

    methods (TestMethodSetup)

        function addToolboxPath(tc)
            root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
            old_path = path;
            tc.addTeardown(@() path(old_path));
            addpath(fullfile(root, 'toolbox'));
        end

    end

    methods (Test)

        function reproducesMultidimensionalPolynomialAndDerivatives(tc)
            t = linspace(0, 4, 51).';
            position = [1 + 2 * t - t.^2 + 0.1 * t.^3, ...
                        -2 + 0.5 * t.^2 - 0.02 * t.^4, 3 + 0 * t];
            velocity = [2 - 2 * t + 0.3 * t.^2, t - 0.08 * t.^3, 0 * t];
            opts = splineOptions();
            sp = coa.interp.fitWindowedSpline(t, position, velocity, opts);
            tq = linspace(0, 4, 201).';
            actual = coa.interp.evaluateSpline(sp, tq, 3);
            expected = zeros(numel(tq), 3, 4);
            expected(:, :, 1) = [1 + 2 * tq - tq.^2 + 0.1 * tq.^3, ...
                                 -2 + 0.5 * tq.^2 - 0.02 * tq.^4, 3 + 0 * tq];
            expected(:, :, 2) = [2 - 2 * tq + 0.3 * tq.^2, ...
                                 tq - 0.08 * tq.^3, 0 * tq];
            expected(:, :, 3) = [-2 + 0.6 * tq, 1 - 0.24 * tq.^2, 0 * tq];
            expected(:, :, 4) = [0.6 + 0 * tq, -0.48 * tq, 0 * tq];
            tc.verifySize(actual, [numel(tq), 3, 4]);
            tc.verifyEqual(actual, expected, 'AbsTol', 1e-8);
        end

        function componentsAreIndependent(tc)
            t = cumsum(0.1 + (1:60).' / 1000);
            position = [sin(t), 5 * cos(0.7 * t), t.^2, zeros(size(t))];
            velocity = [cos(t), -3.5 * sin(0.7 * t), 2 * t, zeros(size(t))];
            opts = splineOptions();
            multi = coa.interp.fitWindowedSpline(t, position, velocity, opts);
            tq = linspace(t(1), t(end), 101).';
            multi_values = coa.interp.evaluateSpline(multi, tq, 3);
            for component = 1:size(position, 2)
                scalar_position = position(:, component);
                scalar_velocity = velocity(:, component);
                scalar_spline = coa.interp.fitWindowedSpline(t, scalar_position, ...
                                                             scalar_velocity, opts);
                single_values = coa.interp.evaluateSpline(scalar_spline, tq, 3);
                tc.verifyEqual(multi_values(:, component, :), single_values, ...
                               'AbsTol', 1e-9);
            end
            tc.verifyEqual(multi_values(:, 4, :), zeros(101, 1, 4));
        end

        function continuousAtEveryKnotAndWindowBoundary(tc)
            t = cumsum(0.08 + (1:60).' / 1000);
            opts = splineOptions();
            sp = coa.interp.fitWindowedSpline(t, [sin(t), cos(t)], ...
                                              [cos(t), -sin(t)], opts);
            tc.verifyGreaterThan(numel(sp.window_starts), 1);
            left = coa.interp.evaluateSpline(sp, t(2:end - 1), 4, 'left');
            right = coa.interp.evaluateSpline(sp, t(2:end - 1), 4, 'right');
            tc.verifyEqual(left(:, :, 1:4), right(:, :, 1:4), 'AbsTol', 1e-8);
            tc.verifyEqual(left(:, :, 5), right(:, :, 5), 'AbsTol', 1e-6);
        end

        function approximatesOffGridSinusoidalDerivatives(tc)
            t = linspace(0, 4, 161).';
            opts = splineOptions();
            sp = coa.interp.fitWindowedSpline(t, [sin(t), cos(2 * t)], ...
                                              [cos(t), -2 * sin(2 * t)], opts);
            tq = linspace(0, 4, 321).';
            actual = coa.interp.evaluateSpline(sp, tq, 3);
            expected = zeros(numel(tq), 2, 4);
            expected(:, :, 1) = [sin(tq), cos(2 * tq)];
            expected(:, :, 2) = [cos(tq), -2 * sin(2 * tq)];
            expected(:, :, 3) = [-sin(tq), -4 * cos(2 * tq)];
            expected(:, :, 4) = [-cos(tq), 8 * sin(2 * tq)];
            tolerances = [1e-7, 1e-6, 1e-4, 1e-2];
            for m = 1:4
                tc.verifyLessThan(max(abs(actual(:, :, m) - expected(:, :, m)), ...
                                      [], 'all'), tolerances(m));
            end
        end

        function supportsSingleWindowAndShiftedTimes(tc)
            t = 2450000 + (0:5).' / 10;
            tau = t - t(1);
            opts = splineOptions();
            sp = coa.interp.fitWindowedSpline(t, [tau, tau.^2], ...
                                              [ones(size(tau)), 2 * tau], opts);
            tc.verifyEqual(numel(sp.window_starts), 1);
            values = coa.interp.evaluateSpline(sp, t([1, end]), 2);
            tc.verifyEqual(values(:, :, 1), [tau([1, end]), tau([1, end]).^2], ...
                           'AbsTol', 1e-10);
            tc.verifyEqual(values(:, :, 2), [1, 0; 1, 2 * tau(end)], ...
                           'AbsTol', 1e-9);
            tc.verifyEqual(values(:, :, 3), [0, 2; 0, 2], 'AbsTol', 1e-8);
        end

        function rejectsInvalidInputs(tc)
            t = (0:9).';
            opts = splineOptions();
            short_velocity = t(1:end - 1);
            short_fit = @() coa.interp.fitWindowedSpline(t, t, short_velocity, opts);
            tc.verifyError(short_fit, ...
                           'coa:interp:fitWindowedSpline:SizeMismatch');
            repeated_times = [t(1:9); t(9)];
            repeated_fit = @() coa.interp.fitWindowedSpline(repeated_times, t, t, opts);
            tc.verifyError(repeated_fit, ...
                           'coa:interp:fitWindowedSpline:InvalidTimes');
            invalid_opts = struct('weight', 1);
            tc.verifyError(@() coa.interp.fitWindowedSpline(t, t, t, invalid_opts), ...
                           'coa:interp:fitWindowedSpline:InvalidOptions');
            sp = coa.interp.fitWindowedSpline(t, t, ones(size(t)), opts);
            tc.verifySize(coa.interp.evaluateSpline(sp, [], 3), [0, 1, 4]);
            tc.verifyError(@() coa.interp.evaluateSpline(sp, -1, 3), ...
                           'coa:interp:evaluateSpline:OutOfRange');
            tc.verifyError(@() coa.interp.evaluateSpline(sp, 1, 3, 'bad'), ...
                           'coa:interp:evaluateSpline:InvalidSide');
            tc.verifyError(@() coa.interp.evaluateSpline(sp, ones(2), 3), ...
                           'coa:interp:evaluateSpline:InvalidQuerySize');
        end

        function defaultsAndPartialOptions(tc)
            t = (0:0.1:5).';
            position = [sin(t), cos(t)];
            velocity = [cos(t), -sin(t)];
            sp = coa.interp.fitWindowedSpline(t, position, velocity);
            defaults = struct('order', 6, 'window_size', 40, 'stride', 20, ...
                              'velocity_weight', 1);
            explicit = coa.interp.fitWindowedSpline(t, position, velocity, defaults);
            tc.verifyEqual(sp, explicit);
            partial = coa.interp.fitWindowedSpline(t, position, velocity, ...
                                                   struct('velocity_weight', 2));
            defaults.velocity_weight = 2;
            explicit = coa.interp.fitWindowedSpline(t, position, velocity, defaults);
            tc.verifyEqual(partial, explicit);
            tc.verifyEqual(coa.interp.evaluateSpline(sp, t), ...
                           coa.interp.evaluateSpline(sp, t, 0));
        end

    end
end

function opts = splineOptions()
    opts = struct('order', 6, 'window_size', 12, 'stride', 5, ...
                  'velocity_weight', 1);
end
