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

    end
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
