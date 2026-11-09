% Copyright (c) 2026 Andrew Binder

classdef EphemIntegratorIntegrationTest < matlab.unittest.TestCase
    methods (TestMethodSetup)

        function addProjectRoot(tc)
            root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
            old_path = path;
            tc.addTeardown(@path, old_path);
            addpath(root);
        end

    end

    methods (Test)

        function fitsAndEvaluatesThroughPinnedDependency(tc)
            root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
            previous_entries = strsplit(path, pathsep);
            setupProjectPaths(root);
            dependency = fullfile(root, 'external', 'ephem_integrator');
            tc.verifyEqual(which('emframe.interp.fitWindowedSpline'), ...
                           fullfile(dependency, 'toolbox', '+emframe', '+interp', ...
                                    'fitWindowedSpline.m'));
            tc.verifyEqual(which('emframe.interp.evaluateSpline'), ...
                           fullfile(dependency, 'toolbox', '+emframe', '+interp', ...
                                    'evaluateSpline.m'));
            added_entries = setdiff(strsplit(path, pathsep), previous_entries);
            tc.verifyFalse(any(startsWith(added_entries, ...
                                          fullfile(dependency, 'tests'))));
            tc.verifyFalse(any(startsWith(added_entries, ...
                                          fullfile(dependency, 'examples'))));

            t = linspace(0, 4, 51).';
            position = [1 + 2 * t - t.^2, t.^3];
            velocity = [2 - 2 * t, 3 * t.^2];
            opts = struct('order', 6, 'window_size', 12, 'stride', 5, ...
                          'velocity_weight', 1);
            sp = emframe.interp.fitWindowedSpline(t, position, velocity, opts);
            tq = linspace(0, 4, 101).';
            actual = emframe.interp.evaluateSpline(sp, tq, 3);
            expected = zeros(numel(tq), 2, 4);
            expected(:, :, 1) = [1 + 2 * tq - tq.^2, tq.^3];
            expected(:, :, 2) = [2 - 2 * tq, 3 * tq.^2];
            expected(:, :, 3) = [-2 + 0 * tq, 6 * tq];
            expected(:, :, 4) = [0 * tq, 6 + 0 * tq];
            tc.verifySize(actual, [numel(tq), 2, 4]);
            tc.verifyEqual(actual, expected, 'AbsTol', 1e-8);
            tc.verifyError(@() emframe.interp.evaluateSpline(sp, -1), ...
                           'emframe:interp:evaluateSpline:OutOfRange');
        end

        function missingSubmoduleReportsHowToInitialize(tc)
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            old_path = path;
            tc.verifyError(@() setupProjectPaths(folder.Folder), ...
                           'coa:setupProjectPaths:MissingDependency');
            tc.verifyEqual(path, old_path);
        end

        function runnableExampleUsesDependency(tc)
            root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
            addpath(fullfile(root, 'examples'));
            [t, values] = ephemSplineUsage();
            tc.verifySize(values, [201, 2, 4]);
            tc.verifyEqual(values(:, :, 1), [sin(t), cos(t)], 'AbsTol', 1e-6);
            tc.verifyEqual(values(:, :, 2), [cos(t), -sin(t)], 'AbsTol', 1e-5);
            tc.verifyEqual(values(:, :, 3), [-sin(t), -cos(t)], 'AbsTol', 1e-4);
            tc.verifyEqual(values(:, :, 4), [-cos(t), sin(t)], 'AbsTol', 1e-2);
        end

    end
end
