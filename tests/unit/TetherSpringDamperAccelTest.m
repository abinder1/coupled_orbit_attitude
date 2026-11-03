% Copyright (c) 2026 Andrew Binder

classdef TetherSpringDamperAccelTest < matlab.unittest.TestCase
    % Tests for coa.orbit.tetherSpringDamperAccel

    properties (TestParameter)
        eta = {0.1, 0.5, 0.9}
        bad_gain = {-1, [1; 2], [], "a", NaN}
    end

    methods (Test)

        function rejectsInvalidInputs(tc, bad_gain)
            [x1, x2, fn] = defaultInputs();
            accel = @coa.orbit.tetherSpringDamperAccel;
            tc.verifyError(@() accel(0, x1, x2, 0.4, 1, 0, bad_gain, 1, fn), ...
                           ?MException);
            tc.verifyError(@() accel(0, x1, x2, 0.4, 1, 0, 1, bad_gain, fn), ...
                           ?MException);
        end

        function rejectsInvalidEtaAndLength(tc)
            [x1, x2, fn] = defaultInputs();
            accel = @coa.orbit.tetherSpringDamperAccel;
            tc.verifyError(@() accel(0, x1, x2, 1.5, 1, 0, 1, 1, fn), ?MException);
            tc.verifyError(@() accel(0, x1, x2, 0.4, -1, 0, 1, 1, fn), ?MException);
        end

        function matchesRigidConstraintWhenOnReference(tc, eta)
            % Zero gains, correct length, rate and spin: identical to tetherAccel.
            [x1, x2] = spinningPair(2, 0.7, 0.4);
            fn = @(x) -x(1:3);
            length = 2;
            expected = coa.orbit.tetherAccel(x1, x2, eta, 0, fn);
            actual = coa.orbit.tetherSpringDamperAccel(0, x1, x2, eta, length, ...
                                                       0, 0, 0, fn);
            tc.verifyEqual(actual, expected, 'AbsTol', 1e-12);
        end

        function conservesMassWeightedAcceleration(tc, eta)
            [x1, x2, fn] = defaultInputs();
            x2(1:3) = x2(1:3) * 1.1;
            a = coa.orbit.tetherSpringDamperAccel(0.3, x1, x2, eta, 0.5, 0.01, ...
                                                  3, 2, fn);
            free = [fn(x1); fn(x2)];
            tc.verifyEqual((1 - eta) * a(1:3) + eta * a(4:6), ...
                           (1 - eta) * free(1:3) + eta * free(4:6), ...
                           'AbsTol', 1e-12);
        end

        function stretchedTetherPullsMassesTogether(tc)
            x1 = [0; 0; 0; 0; 0; 0];
            x2 = [1.2; 0; 0; 0; 0; 0];
            a = coa.orbit.tetherSpringDamperAccel(0, x1, x2, 0.5, 1, 0, 4, 0, ...
                                                  @freeSpace);
            tc.verifyGreaterThan(a(1), 0);
            tc.verifyLessThan(a(4), 0);
            tc.verifyEqual(a(4) - a(1), -4 * 0.2, 'AbsTol', 1e-12);
        end

        function returnsTensionAlongLineOfSight(tc)
            x1 = zeros(6, 1);
            x2 = [1.2; 0; 0; 0; 0; 0];
            fn = @freeSpaceWithJacobians;
            [a, ~, ~, tension] = coa.orbit.tetherSpringDamperAccel(0, x1, x2, 0.5, ...
                                                                   1, 0, 4, 0, fn);
            tc.verifyEqual(tension, 4 * 0.2, 'AbsTol', 1e-12);
            tc.verifyEqual(a(4) - a(1), -tension, 'AbsTol', 1e-12);
        end

        function dampsRadialVelocityError(tc)
            x1 = zeros(6, 1);
            x2 = [1; 0; 0; 0.3; 0; 0];
            a = coa.orbit.tetherSpringDamperAccel(0, x1, x2, 0.5, 1, 0, 0, 5, ...
                                                  @freeSpace);
            tc.verifyEqual(a(4) - a(1), -5 * 0.3, 'AbsTol', 1e-12);
        end

        function convergesToReferenceLengthWhileSpinning(tc)
            % Free space, 5 percent stretched and spinning, critically damped.
            length_0 = 2;
            [x1, x2] = spinningPair(1.05 * length_0, 0.7, 0.5);
            rhs = @(t, x) derivative(t, x, length_0, 0, 9, 6);
            [~, y] = ode45(rhs, [0, 12], [x1; x2], integratorOptions());
            final_length = norm(y(end, 7:9) - y(end, 1:3));
            tc.verifyEqual(final_length, length_0, 'AbsTol', 1e-6);
        end

        function tracksRampingReferenceLength(tc)
            % Starting on the reference ramp, the length error stays near zero.
            length_0 = 2;
            length_dot = 0.1;
            [x1, x2] = spinningPair(length_0, 0.7, 0.5);
            x1(4) = x1(4) - 0.5 * length_dot;
            x2(4) = x2(4) + 0.5 * length_dot;
            rhs = @(t, x) derivative(t, x, length_0, length_dot, 9, 6);
            [t, y] = ode45(rhs, [0, 5], [x1; x2], integratorOptions());
            lengths = vecnorm(y(:, 7:9) - y(:, 1:3), 2, 2);
            tc.verifyEqual(lengths, length_0 + length_dot * t, 'AbsTol', 1e-6);
        end

        function jacobiansMatchFiniteDifferences(tc, eta)
            % Off-constraint state, nonzero gains, and the full CR3BP field.
            [x1, x2, fn] = jacobianInputs();
            args = {eta, 0.1, 1e-3, 3, 4, fn};
            [~, da_dp, da_dv] = coa.orbit.tetherSpringDamperAccel(0.2, x1, x2, ...
                                                                  args{:});
            [numeric_dp, numeric_dv] = finiteDifferenceJacobians(x1, x2, args);
            scale = max(abs([numeric_dp(:); numeric_dv(:)]));
            tc.verifyEqual(da_dp, numeric_dp, 'AbsTol', 1e-5 * scale);
            tc.verifyEqual(da_dv, numeric_dv, 'AbsTol', 1e-5 * scale);
        end

        function jacobiansMatchFiniteDifferencesWithZeroGains(tc)
            [x1, x2, fn] = jacobianInputs();
            args = {0.4, 0.1, 1e-3, 0, 0, fn};
            [~, da_dp, da_dv] = coa.orbit.tetherSpringDamperAccel(0.2, x1, x2, ...
                                                                  args{:});
            [numeric_dp, numeric_dv] = finiteDifferenceJacobians(x1, x2, args);
            scale = max(abs([numeric_dp(:); numeric_dv(:)]));
            tc.verifyEqual(da_dp, numeric_dp, 'AbsTol', 1e-5 * scale);
            tc.verifyEqual(da_dv, numeric_dv, 'AbsTol', 1e-5 * scale);
        end

    end
end

function [x1, x2, fn] = jacobianInputs()
    x1 = [0.6; 0.1; 0.12; 0.02; 0.24; 0.01];
    x2 = [0.67; 0.04; 0.17; 0.03; 0.2; 0.02];
    fn = @(x) coa.orbit.cr3bp(x, 0.0121);
end

function [da_dp, da_dv] = finiteDifferenceJacobians(x1, x2, args)
    % Central differences of the stacked acceleration with respect to
    % [p1; p2] and [v1; v2].
    x = [x1; x2];
    step = 1e-6;
    numeric = zeros(6, 12);
    for k = 1:12
        shift = zeros(12, 1);
        shift(k) = step;
        numeric(:, k) = (stackedAccel(x + shift, args) - ...
                         stackedAccel(x - shift, args)) / (2 * step);
    end
    % x is ordered [p1; v1; p2; v2].
    da_dp = numeric(:, [1:3, 7:9]);
    da_dv = numeric(:, [4:6, 10:12]);
end

function a = stackedAccel(x, args)
    a = coa.orbit.tetherSpringDamperAccel(0.2, x(1:6), x(7:12), args{:});
end

function a = freeSpace(~)
    a = zeros(3, 1);
end

function [a, da_dp, da_dv] = freeSpaceWithJacobians(~)
    a = zeros(3, 1);
    da_dp = zeros(3);
    da_dv = zeros(3);
end

function [x1, x2] = spinningPair(length, spin, eta)
    % Masses on the x axis about their center of mass, spinning about z.
    x1 = [-eta * length; 0; 0; 0; -eta * length * spin; 0];
    x2 = [(1 - eta) * length; 0; 0; 0; (1 - eta) * length * spin; 0];
end

function dxdt = derivative(t, x, length_0, length_dot, stiffness, damping)
    % Free-space dynamics of the packed state [x1; x2].
    a = coa.orbit.tetherSpringDamperAccel(t, x(1:6), x(7:12), 0.5, length_0, ...
                                          length_dot, stiffness, damping, ...
                                          @freeSpace);
    dxdt = [x(4:6); a(1:3); x(10:12); a(4:6)];
end

function options = integratorOptions()
    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
end

function [x1, x2, fn] = defaultInputs()
    x1 = [0.3; 0.2; 0.1; 0.01; -0.03; 0.04];
    x2 = [0.35; 0.1; 0.15; -0.02; 0.05; 0.0];
    fn = @(x) -x(1:3);
end
