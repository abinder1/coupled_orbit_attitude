% Copyright (c) 2026 Andrew Binder

classdef TetherFilterDemoEphemerisTruthTest < matlab.unittest.TestCase
    % Tests for the demo's time-dependent ephemeris truth dynamics.
    methods (TestMethodSetup)

        function addDependencyPaths(tc)
            root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
            original_path = path;
            addpath(fullfile(root, 'toolbox'));
            addpath(fullfile(root, 'examples'));
            addpath(fullfile(root, 'external', 'ephem_integrator', 'toolbox'));
            addpath(fullfile(root, 'external', 'ephem_integrator', 'tests', 'support'));
            tc.addTeardown(@path, original_path);
        end

    end

    methods (Test)

        function zeroTetherForceMatchesPublicDerivative(tc)
            [parameters, config, state] = inputs();
            truth_rhs = @tetherFilterDemoSupport.ephemerisTruthDynamics;
            stiffness = 100;
            config.spring_damper.stiffness = stiffness;
            config.spring_damper.damping = 0;
            config.tether_length_rate = 0;

            for tau = [0, 0.12, 0.3]
                free = freeDerivatives(tau, state, parameters);
                separation = state(4:6) - state(1:3);
                relative_velocity = state(13:15) - state(10:12);
                distance = norm(separation);
                direction = separation / distance;
                length_rate = dot(direction, relative_velocity);
                natural_length_acceleration = ...
                    dot(direction, free(13:15) - free(10:12)) + ...
                    (dot(relative_velocity, relative_velocity) - length_rate^2) / ...
                    distance;
                config.spring_damper.initial_length = ...
                    distance + natural_length_acceleration / stiffness;

                actual = truth_rhs(tau, state, config, parameters);
                tc.verifyEqual(actual, free, 'AbsTol', 1e-11);
            end
        end

        function tetherForceCancelsAndPayloadIsFree(tc)
            [parameters, config, state] = inputs();
            truth_rhs = @tetherFilterDemoSupport.ephemerisTruthDynamics;
            eta = config.mass_2_fraction;

            for use_spring = [false, true]
                if ~use_spring
                    config.spring_damper = struct();
                else
                    config.spring_damper = struct('initial_length', 0.08, ...
                                                  'stiffness', 3, 'damping', 4);
                end
                for tau = [0, 0.12, 0.3]
                    free = freeDerivatives(tau, state, parameters);
                    actual = truth_rhs(tau, state, config, parameters);
                    correction_1 = actual(10:12) - free(10:12);
                    correction_2 = actual(13:15) - free(13:15);

                    tc.verifyEqual(actual(1:9), state(10:18), 'AbsTol', 0);
                    tc.verifyEqual(actual(16:18), free(16:18), 'AbsTol', 1e-12);
                    tc.verifyGreaterThan(norm(correction_1), 1e-8);
                    tc.verifyEqual((1 - eta) * correction_1 + eta * correction_2, ...
                                   zeros(3, 1), 'AbsTol', 1e-12);
                    tc.verifyEqual(cross(state(4:6) - state(1:3), correction_1), ...
                                   zeros(3, 1), 'AbsTol', 1e-12);
                end
            end
        end

        function payloadAccelerationAgreesInInertialAxes(tc)
            [parameters, config, state] = inputs();
            truth_rhs = @tetherFilterDemoSupport.ephemerisTruthDynamics;
            frame = @emframe.orbit.earthMoonFrameKinematics;
            inertial_rhs = @emframe.embInertialDerivative;
            for tau = [0, 0.12, 0.3]
                derivative = truth_rhs(tau, state, config, parameters);
                [basis, omega, omega_dot] = frame(parameters.ephemeris.moon(tau));
                position = state(7:9);
                velocity = state(16:18);
                inertial_state = [basis * position; ...
                                  basis * (velocity + cross(omega, position))];
                inertial_derivative = inertial_rhs(tau, inertial_state, parameters);
                transport = 2 * cross(omega, velocity) + ...
                            cross(omega_dot, position) + ...
                            cross(omega, cross(omega, position));
                inertial_acceleration = basis * (derivative(16:18) + transport);

                tc.verifyEqual(inertial_acceleration, inertial_derivative(4:6), ...
                               'AbsTol', 1e-11);
            end
        end

        function usesSolverTimeAndRejectsOutOfWindow(tc)
            [parameters, config, state] = inputs();
            truth_rhs = @tetherFilterDemoSupport.ephemerisTruthDynamics;
            first = truth_rhs(0, state, config, parameters);
            later = truth_rhs(0.3, state, config, parameters);
            tc.verifyGreaterThan(norm(later - first), 1e-8);

            rhs = @() truth_rhs(0.5, state, config, parameters);
            tc.verifyError(rhs, 'emframe:validation:PreparedOutOfCoverage');
        end

        function scenarioSelectsEphemerisOnlyForTruth(tc)
            import tetherFilterDemoSupport.*
            [parameters, config, state] = inputs();
            config.truth_model = "ephemeris";
            config.ephemeris_parameters = parameters;
            config.mu = 0.012150585;
            noise = struct('tension_psd', 1e-6, 'transverse_psd', 1e-8, ...
                           'body_3_acceleration_psd', 1e-8);
            truth_rhs = Scenario.truthDynamicsFunction(config);
            filter_rhs = Scenario.filterDynamicsFunction(config, noise);
            cr3bp_rhs = Scenario.cr3bpDynamicsFunction(config, noise);
            expected_truth = ephemerisTruthDynamics(0.12, state, config, parameters);
            packed = [state; zeros(324, 1)];

            tc.verifyEqual(truth_rhs(0.12, state), expected_truth, 'AbsTol', 0);
            tc.verifyEqual(filter_rhs(0.12, packed), cr3bp_rhs(0.12, packed), ...
                           'AbsTol', 0);
        end

        function rendezvousUsesAbsoluteEphemerisTime(tc)
            import tetherFilterDemoSupport.*
            [parameters, config, state] = inputs();
            config.truth_model = "ephemeris";
            config.ephemeris_parameters = parameters;
            config.duration = 0.03;
            tether_state = state([1:6, 10:15]);
            initial = Scenario.rendezvousInitialState(tether_state, config);
            truth_rhs = Scenario.truthDynamicsFunction(config);
            truth = ode113(truth_rhs, [0, config.duration], initial, ...
                           odeset('RelTol', 1e-10, 'AbsTol', 1e-12));
            final = truth.y(:, end);

            tc.verifyEqual(initial([1:6, 10:15]), tether_state, 'AbsTol', 0);
            tc.verifyEqual(final(7:9), final(4:6), 'AbsTol', 1e-8);
            tc.verifyEqual(final(16:18), final(13:15), 'AbsTol', 1e-8);
        end

        function cr3bpRendezvousMatchesPreviousTimeConvention(tc)
            import tetherFilterDemoSupport.*
            [~, config, state] = inputs();
            config.truth_model = "cr3bp";
            config.mu = 0.012150585;
            config.duration = 0.03;
            tether_state = state([1:6, 10:15]);
            initial = Scenario.rendezvousInitialState(tether_state, config);
            truth_rhs = Scenario.truthDynamicsFunction(config);
            truth = ode113(truth_rhs, [0, config.duration], initial, ...
                           odeset('RelTol', 1e-10, 'AbsTol', 1e-12));
            final = truth.y(:, end);
            catch_state = [final(4:6); final(13:15)];
            old_rhs = @(~, x) [x(4:6); coa.orbit.cr3bp(x, config.mu)];
            old_solution = ode113(old_rhs, [0, -config.duration], catch_state, ...
                                  odeset('AbsTol', 1e-10, 'RelTol', 1e-10));

            tc.verifyEqual(initial([7:9, 16:18]), old_solution.y(:, end), ...
                           'AbsTol', 1e-9);
        end

        function preparesScaledDefaultForceFromDeConstants(tc)
            import tetherFilterDemoSupport.*
            fixture = refinementPointwiseInputs();
            data = addDeConstants(fixture);
            config = struct('ephemeris_start_et', 0, 'duration_s', 86400, ...
                            'length_scale', 384400, 'time_scale', 375190.258892627);
            parameters = Ephemeris.prepareParameters(config, data);
            force = emframe.internal.makeForceModel(fixture.constants);
            scales = struct('length', config.length_scale, 'time', config.time_scale);
            reference = emframe.prepareScaledDynamics(data, [0, 86400], force, scales);
            tau = 0.1;
            state = [0.6; 0.1; 0.12; 0.02; 0.24; 0.01];
            actual = emframe.embRotatingDerivative(tau, state, parameters);
            expected = emframe.embRotatingDerivative(tau, state, reference);

            tc.verifyEqual(actual, expected, 'RelTol', 1e-12, 'AbsTol', 1e-12);
            tc.verifyError(@() parameters.ephemeris.moon(-0.01), ...
                           'emframe:validation:PreparedOutOfCoverage');
            tc.verifyError(@() parameters.ephemeris.moon(86400 / config.time_scale + ...
                                                         1e-6), ...
                           'emframe:validation:PreparedOutOfCoverage');

            data.constantNames = data.constantNames(2:end);
            data.constantValues = data.constantValues(2:end);
            tc.verifyError(@() Ephemeris.prepareParameters(config, data), ...
                           'tetherFilterDemoSupport:Ephemeris:MissingConstant');
        end

        function initialTetherSpinUsesActualFrameRate(tc)
            import tetherFilterDemoSupport.*
            [parameters, config, ~] = inputs();
            config.truth_model = "ephemeris";
            config.ephemeris_parameters = parameters;
            config.length_scale = 384400;
            config.time_scale = 375190.258892627;
            config.tether_length_km = 200;
            config.spin_rate = 2 * pi / (55 * 60);
            config.tether_direction = [1; 2; 3] / sqrt(14);
            config.tether_length_rate = 0.001;
            config.halo_state = [0.6; 0.1; 0.12; 0.02; 0.24; 0.01];
            config.duration = 0.003;
            config.initial_sigma = struct('com_position', 100, 'com_velocity', 1, ...
                                          'relative_position', 0.01, ...
                                          'relative_velocity', 0.1, ...
                                          'payload_position', 100, ...
                                          'payload_velocity', 1);
            initial = Scenario.makeInitialEstimate(config);
            [~, frame_rate] = ...
                emframe.orbit.earthMoonFrameKinematics(parameters.ephemeris.moon(0));
            separation = initial(4:6) - initial(1:3);
            rotating_relative_velocity = initial(13:15) - initial(10:12);
            inertial_relative_velocity = rotating_relative_velocity + ...
                                         cross(frame_rate, separation);
            spin = [0; 0; config.spin_rate * config.time_scale];
            expected = cross(spin, separation) + ...
                       config.tether_length_rate * separation / norm(separation);
            eta = config.mass_2_fraction;

            tc.verifyEqual(inertial_relative_velocity, expected, 'AbsTol', 1e-12);
            tc.verifyEqual((1 - eta) * initial(1:3) + eta * initial(4:6), ...
                           config.halo_state(1:3), 'AbsTol', 1e-14);
            tc.verifyEqual((1 - eta) * initial(10:12) + eta * initial(13:15), ...
                           config.halo_state(4:6), 'AbsTol', 1e-14);
            tc.verifyGreaterThan(abs(frame_rate(1)), 1e-8);

            config.truth_model = "cr3bp";
            config.mu = 0.012150585;
            initial = Scenario.makeInitialEstimate(config);
            initializer = @coa.orbit.tetherInitialState;
            legacy = initializer(config.halo_state, ...
                                 config.tether_length_km / config.length_scale, ...
                                 eta, config.spin_rate * config.time_scale, ...
                                 config.tether_length_rate, config.tether_direction);
            tc.verifyEqual(initial([1:6, 10:15]), legacy, 'AbsTol', 0);
        end

        function rejectsMissingOrInvalidLocalConfig(tc)
            import tetherFilterDemoSupport.*
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            config_path = fullfile(folder.Folder, 'config.json');
            tc.verifyError(@() Ephemeris.loadParameters(struct(), config_path), ...
                           'tetherFilterDemoSupport:Ephemeris:MissingConfig');
            writelines("{}", config_path);
            tc.verifyError(@() Ephemeris.loadParameters(struct(), config_path), ...
                           'tetherFilterDemoSupport:Ephemeris:InvalidConfig');
            writelines('{"EMFRAME_DE440_HEADER":"","EMFRAME_DE440_ASCII":"x"}', ...
                       config_path);
            tc.verifyError(@() Ephemeris.loadParameters(struct(), config_path), ...
                           'tetherFilterDemoSupport:Ephemeris:InvalidConfig');
        end

    end
end

function data = addDeConstants(fixture)
    data = fixture.ephemeris;
    gm = fixture.constants.gm;
    combined_gm = [gm(1:3), gm(4) + gm(5), gm(6:11)];
    au = 149597870.7;
    data.constantNames = {'GMS', 'GM1', 'GM2', 'GMB', 'GM4', 'GM5', ...
                          'GM6', 'GM7', 'GM8', 'GM9', 'AU', 'EMRAT'};
    data.constantValues = [combined_gm * 86400^2 / au^3, au, gm(4) / gm(5)];
end

function [parameters, config, state] = inputs()
    fixture = refinementPointwiseInputs();
    force = emframe.internal.makeForceModel(fixture.constants, [10, 399, 301], true);
    scales = struct('length', 384400, 'time', 375190.258892627);
    parameters = emframe.prepareScaledDynamics(fixture.ephemeris, ...
                                               [0, 172800], force, scales);
    config.mass_2_fraction = 0.3;
    config.tether_length_rate = 0;
    config.spring_damper = struct('initial_length', 0.08, ...
                                  'stiffness', 3, 'damping', 4);
    state = [0.6; 0.1; 0.12; 0.67; 0.04; 0.17; 0.8; 0.2; 0.1; ...
             0.02; 0.24; 0.01; 0.03; 0.2; 0.02; -0.02; 0.15; 0.1];
end

function derivative = freeDerivatives(tau, state, parameters)
    derivative = zeros(18, 1);
    for body = 1:3
        position_indices = (3 * body - 2):(3 * body);
        velocity_indices = position_indices + 9;
        body_state = [state(position_indices); state(velocity_indices)];
        body_derivative = emframe.embRotatingDerivative(tau, body_state, parameters);
        derivative(position_indices) = body_derivative(1:3);
        derivative(velocity_indices) = body_derivative(4:6);
    end
end
