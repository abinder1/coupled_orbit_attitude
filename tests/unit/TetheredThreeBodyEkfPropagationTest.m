% Copyright (c) 2026 Andrew Binder

classdef TetheredThreeBodyEkfPropagationTest < matlab.unittest.TestCase
    % Tests for coa.orbit.tetheredThreeBodyEkfPropagation

    properties (TestParameter)
        mass_2_fraction = {0.1, 0.5, 0.9}
        tether_length_rate = {0, 0.05, -0.2}
        state = {
                 [0.3; 0.2; 0.1; 0.35; 0.1; 0.15; -0.5; 0.4; 0.2; ...
                  0.01; -0.03; 0.04; -0.02; 0.05; 0.0; 0.1; 0.0; -0.1], ...
                 [-0.8; 0.5; -0.2; -0.7; 0.45; -0.25; 0.6; -0.3; 0.1; ...
                  0.1; 0.2; -0.3; 0.0; 0.1; 0.2; -0.05; 0.15; 0.0]
                }
        bad_packed_state = {zeros(17, 1), zeros(19, 1), zeros(341, 1), ...
                            zeros(343, 1), zeros(18, 2), zeros(1, 18), ...
                            struct("x", zeros(18, 1))}
        bad_scalar = {[0.1; 0.2], [], struct("a", 1), {0.1}}
        missing_field = {"tension_psd", "transverse_psd", "body_3_acceleration_psd"}
    end

    methods (Test)

        % (1) Input validation

        function rejectsInvalidPackedState(tc, bad_packed_state)
            [~, ~, fn, noise] = defaultInputs();
            tc.verifyError(@() ekf(0, bad_packed_state, 0.3, 0.1, fn, noise), ...
                           ?MException);
        end

        function rejectsInvalidScalars(tc, bad_scalar)
            [x, ~, fn] = defaultInputs();
            tc.verifyError(@() ekf(bad_scalar, x, 0.3, 0.1, fn), ?MException);
            tc.verifyError(@() ekf(0, x, bad_scalar, 0.1, fn), ?MException);
            tc.verifyError(@() ekf(0, x, 0.3, bad_scalar, fn), ?MException);
        end

        function rejectsNonFunctionHandle(tc)
            x = defaultInputs();
            tc.verifyError(@() ekf(0, x, 0.3, 0.1, 1), ?MException);
        end

        function rejectsTooFewInputs(tc)
            x = defaultInputs();
            tc.verifyError(@() ekf(0, x, 0.3, 0.1), ?MException);
        end

        function rejectsNonStructProcessNoise(tc)
            [x, P, fn] = defaultInputs();
            tc.verifyError(@() ekf(0, [x; P(:)], 0.3, 0.1, fn, 1), ?MException);
        end

        function covarianceModeRequiresEveryPsdField(tc, missing_field)
            [x, P, fn, noise] = defaultInputs();
            noise = rmfield(noise, missing_field);
            tc.verifyError(@() ekf(0, [x; P(:)], 0.3, 0.1, fn, noise), ...
                           ?MException);
        end

        function stateOnlyModeIgnoresProcessNoise(tc)
            [x, ~, fn] = defaultInputs();
            without_noise = ekf(0, x, 0.3, 0.1, fn);
            with_bad_noise = ekf(0, x, 0.3, 0.1, fn, struct("unrelated", 1));
            tc.verifyEqual(with_bad_noise, without_noise);
        end

        % (2) Output sizes and mode consistency

        function outputSizes(tc, state, mass_2_fraction, tether_length_rate)
            [~, P, fn, noise] = defaultInputs();
            d_state = ekf(0, state, mass_2_fraction, tether_length_rate, fn);
            d_packed = ekf(0, [state; P(:)], mass_2_fraction, ...
                           tether_length_rate, fn, noise);
            tc.verifySize(d_state, [18, 1]);
            tc.verifySize(d_packed, [342, 1]);
        end

        function stateDerivativeIndependentOfMode(tc, state, mass_2_fraction, ...
                                                  tether_length_rate)
            [~, P, fn, noise] = defaultInputs();
            d_state = ekf(0, state, mass_2_fraction, tether_length_rate, fn);
            d_packed = ekf(0, [state; P(:)], mass_2_fraction, ...
                           tether_length_rate, fn, noise);
            tc.verifyEqual(d_packed(1:18), d_state);
        end

        function ignoresNormalizedTime(tc, state)
            [~, P, fn, noise] = defaultInputs();
            packed = [state; P(:)];
            d0 = ekf(0, packed, 0.3, 0.1, fn, noise);
            d1 = ekf(7.5, packed, 0.3, 0.1, fn, noise);
            tc.verifyEqual(d1, d0);
        end

        function worksAsOdeFunction(tc, state)
            % The signature must be usable directly by the ODE suite.
            [~, ~, fn] = defaultInputs();
            rhs = @(t, y) ekf(t, y, 0.3, 0, fn);
            [~, y] = ode45(rhs, [0, 1e-3], state);
            tc.verifySize(y, [size(y, 1), 18]);
            tc.verifyTrue(all(isfinite(y), "all"));
        end

        % (3) State derivative

        function stateDerivativeMatchesBodyModels(tc, state, mass_2_fraction, ...
                                                  tether_length_rate)
            [~, ~, fn] = defaultInputs();
            d = ekf(0, state, mass_2_fraction, tether_length_rate, fn);

            [x1, x2, x3] = bodyStates(state);
            tether_acceleration = coa.orbit.tetherAccel(x1, x2, ...
                                                        mass_2_fraction, ...
                                                        tether_length_rate, fn);

            tc.verifyEqual(d(1:9), state(10:18));
            tc.verifyEqual(d(10:15), tether_acceleration);
            tc.verifyEqual(d(16:18), fn(x3));
        end

        function body3IsDecoupledFromTether(tc, state)
            % Moving bodies 1 and 2 must not change body 3's derivative.
            [~, ~, fn] = defaultInputs();
            d = ekf(0, state, 0.3, 0.1, fn);
            moved = state;
            moved([1:6, 10:15]) = moved([1:6, 10:15]) + 0.05;
            d_moved = ekf(0, moved, 0.3, 0.1, fn);
            tc.verifyEqual(d_moved([7:9, 16:18]), d([7:9, 16:18]));
        end

        % (4) Linearized dynamics F, via the zero-noise covariance derivative

        function dynamicsJacobianMatchesCentralDifference(tc, state, ...
                                                          mass_2_fraction, ...
                                                          tether_length_rate)
            [~, ~, fn, ~, no_noise] = defaultInputs();
            jacobian_fd = recoverDynamicsJacobian(state, mass_2_fraction, ...
                                                  tether_length_rate, fn);

            % With P = I and no process noise, P_dot = F + F'.
            packed = [state; reshape(eye(18), [], 1)];
            d = ekf(0, packed, mass_2_fraction, tether_length_rate, fn, no_noise);
            covariance_derivative = reshape(d(19:end), [18, 18]);

            tc.verifyEqual(covariance_derivative, jacobian_fd + jacobian_fd', ...
                           "AbsTol", 1e-6, "RelTol", 1e-6);
        end

        function dynamicsJacobianHasKinematicStructure(tc, state)
            % Top-right 9x9 block of F is I, top-left is 0: positions are
            % driven only by velocities.
            [~, ~, fn] = defaultInputs();
            jacobian = recoverDynamicsJacobian(state, 0.3, 0.1, fn);
            tc.verifyEqual(jacobian(1:9, 1:9), zeros(9), "AbsTol", 1e-8);
            tc.verifyEqual(jacobian(1:9, 10:18), eye(9), "AbsTol", 1e-8);
        end

        function dynamicsJacobianDecouplesBody3(tc, state)
            [~, ~, fn] = defaultInputs();
            jacobian = recoverDynamicsJacobian(state, 0.3, 0.1, fn);
            body_3 = [7:9, 16:18];
            tethered = [1:6, 10:15];
            tc.verifyEqual(jacobian(body_3, tethered), zeros(6, 12), "AbsTol", 1e-8);
            tc.verifyEqual(jacobian(tethered, body_3), zeros(12, 6), "AbsTol", 1e-8);
        end

        % (5) Covariance propagation

        function covarianceDerivativeMatchesLyapunovEquation(tc, state, ...
                                                             mass_2_fraction, ...
                                                             tether_length_rate)
            [~, ~, fn, noise] = defaultInputs();
            covariance = randomCovariance();

            d = ekf(0, [state; covariance(:)], mass_2_fraction, ...
                    tether_length_rate, fn, noise);
            covariance_derivative = reshape(d(19:end), [18, 18]);

            jacobian = recoverDynamicsJacobian(state, mass_2_fraction, ...
                                               tether_length_rate, fn);
            expected = jacobian * covariance + covariance * jacobian' + ...
                       independentNoiseCovariance(state, mass_2_fraction, noise);

            % Tolerance reflects the finite-difference error in F.
            tc.verifyEqual(covariance_derivative, expected, ...
                           "AbsTol", 1e-8, "RelTol", 1e-6);
        end

        function covarianceDerivativeIsSymmetric(tc, state)
            [~, ~, fn, noise] = defaultInputs();
            covariance = randomCovariance();
            d = ekf(0, [state; covariance(:)], 0.3, 0.1, fn, noise);
            covariance_derivative = reshape(d(19:end), [18, 18]);
            tc.verifyEqual(covariance_derivative, covariance_derivative', ...
                           "AbsTol", 1e-12);
        end

        function zeroCovarianceGivesProcessNoiseOnly(tc, state, mass_2_fraction)
            [~, ~, fn, noise] = defaultInputs();
            d = ekf(0, [state; zeros(324, 1)], mass_2_fraction, 0.1, fn, noise);
            covariance_derivative = reshape(d(19:end), [18, 18]);
            expected = independentNoiseCovariance(state, mass_2_fraction, noise);
            tc.verifyEqual(covariance_derivative, expected, "AbsTol", 1e-12);
        end

        function processNoiseIsPositiveSemidefinite(tc, state)
            [~, ~, fn, noise] = defaultInputs();
            d = ekf(0, [state; zeros(324, 1)], 0.3, 0.1, fn, noise);
            noise_covariance = reshape(d(19:end), [18, 18]);
            tc.verifyGreaterThanOrEqual(min(eig(noise_covariance)), -1e-14);
        end

        function eachPsdScalesItsOwnContribution(tc, state)
            % Doubling one PSD must change only that source's contribution.
            [~, ~, fn, noise] = defaultInputs();
            packed = [state; zeros(324, 1)];
            base = noiseTerm(packed, fn, noise);

            for field = ["tension_psd", "transverse_psd", "body_3_acceleration_psd"]
                doubled = noise;
                doubled.(field) = 2 * noise.(field);
                zeroed = noise;
                zeroed.(field) = 0;
                only_this = base - noiseTerm(packed, fn, zeroed);
                tc.verifyEqual(noiseTerm(packed, fn, doubled), base + only_this, ...
                               "AbsTol", 1e-14);
            end
        end

        function transverseNoiseHasNoRadialComponent(tc, state)
            % P_perp * e = 0, so transverse noise must not drive the
            % relative acceleration along the tether (Eq. 25 of the .tex).
            [~, ~, fn] = defaultInputs();
            only_transverse = struct("tension_psd", 0, "transverse_psd", 1, ...
                                     "body_3_acceleration_psd", 0);
            noise_covariance = noiseTerm([state; zeros(324, 1)], fn, only_transverse);

            tether_direction = state(4:6) - state(1:3);
            tether_direction = tether_direction / norm(tether_direction);
            relative_map = [-eye(3), eye(3)];
            velocity_rows = 10:15;
            relative_noise = relative_map * ...
                noise_covariance(velocity_rows, velocity_rows) * relative_map';
            tc.verifyEqual(tether_direction' * relative_noise * tether_direction, ...
                           0, "AbsTol", 1e-14);
        end

        function tensionNoiseIsAlongTether(tc, state)
            [~, ~, fn] = defaultInputs();
            only_tension = struct("tension_psd", 1, "transverse_psd", 0, ...
                                  "body_3_acceleration_psd", 0);
            noise_covariance = noiseTerm([state; zeros(324, 1)], fn, only_tension);

            tether_direction = state(4:6) - state(1:3);
            tether_direction = tether_direction / norm(tether_direction);
            transverse_basis = null(tether_direction');
            body_1_noise = noise_covariance(10:12, 10:12);
            tc.verifyEqual(transverse_basis' * body_1_noise * transverse_basis, ...
                           zeros(2), "AbsTol", 1e-14);
        end

        function body3NoiseOnlyEntersBody3Velocity(tc, state)
            [~, ~, fn] = defaultInputs();
            only_body_3 = struct("tension_psd", 0, "transverse_psd", 0, ...
                                 "body_3_acceleration_psd", 2);
            noise_covariance = noiseTerm([state; zeros(324, 1)], fn, only_body_3);

            expected = zeros(18);
            expected(16:18, 16:18) = 2 * eye(3);
            tc.verifyEqual(noise_covariance, expected, "AbsTol", 1e-14);
        end

    end
end

function d = ekf(varargin)
    % Short alias that keeps the test call sites within the line-length limit.
    d = coa.orbit.tetheredThreeBodyEkfPropagation(varargin{:});
end

function [x, covariance, fn, noise, no_noise] = defaultInputs()
    x = [0.3; 0.2; 0.1; 0.35; 0.1; 0.15; -0.5; 0.4; 0.2; ...
         0.01; -0.03; 0.04; -0.02; 0.05; 0.0; 0.1; 0.0; -0.1];
    covariance = 1e-4 * eye(18);
    fn = @(state) coa.orbit.cr3bp(state, 0.0121505856);
    noise = struct("tension_psd", 1e-6, "transverse_psd", 2e-7, ...
                   "body_3_acceleration_psd", 3e-7);
    no_noise = struct("tension_psd", 0, "transverse_psd", 0, ...
                      "body_3_acceleration_psd", 0);
end

function [x1, x2, x3] = bodyStates(state)
    x1 = state([1:3, 10:12]);
    x2 = state([4:6, 13:15]);
    x3 = state([7:9, 16:18]);
end

function covariance = randomCovariance()
    % Deterministic, symmetric positive definite.
    generator = RandStream("mt19937ar", "Seed", 7);
    factor = randn(generator, 18);
    covariance = 1e-3 * (factor * factor' + 18 * eye(18));
end

function jacobian = recoverDynamicsJacobian(state, mass_2_fraction, ...
                                            tether_length_rate, fn)
    % F by central differences of the state-only derivative, independent of
    % the analytic Jacobians used in the covariance path.
    f = @(x) ekf(0, x, mass_2_fraction, tether_length_rate, fn);
    jacobian = centralDifferenceJacobian(f, state);
end

function noise_covariance = independentNoiseCovariance(state, mass_2_fraction, noise)
    % G*Q*G' built directly from Eqs. 20 and 26 of ekf_equations.tex.
    separation = state(4:6) - state(1:3);
    e = separation / norm(separation);
    projector = eye(3) - e * e';
    mass_1_fraction = 1 - mass_2_fraction;

    tension_shape = zeros(18, 1);
    tension_shape(10:12) = e / mass_1_fraction;
    tension_shape(13:15) = -e / mass_2_fraction;

    transverse_shape = zeros(18, 3);
    transverse_shape(10:12, :) = projector / mass_1_fraction;
    transverse_shape(13:15, :) = -projector / mass_2_fraction;

    body_3_shape = zeros(18, 3);
    body_3_shape(16:18, :) = eye(3);

    noise_covariance = ...
        noise.tension_psd * (tension_shape * tension_shape') + ...
        noise.transverse_psd * (transverse_shape * transverse_shape') + ...
        noise.body_3_acceleration_psd * (body_3_shape * body_3_shape');
end

function noise_covariance = noiseTerm(packed, fn, noise)
    d = ekf(0, packed, 0.3, 0.1, fn, noise);
    noise_covariance = reshape(d(19:end), [18, 18]);
end

function jacobian_fd = centralDifferenceJacobian(f, x)
    n = numel(x);
    jacobian_fd = zeros(numel(f(x)), n);
    for k = 1:n
        h = 1e-6 * max(1, abs(x(k)));
        e = zeros(n, 1);
        e(k) = h;
        jacobian_fd(:, k) = (f(x + e) - f(x - e)) / (2 * h);
    end
end
