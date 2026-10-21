% Copyright (c) 2026 Andrew Binder

classdef TetherAccelTest < matlab.unittest.TestCase
    % Tests for coa.orbit.tetherAccel

    properties (TestParameter)
        eta = {0.1, 0.5, 0.9}
        length_dot = {0, 0.05, -0.2}
        states = {
                  {[0.3; 0.2; 0.1; 0.01; -0.03; 0.04], ...
                   [0.35; 0.1; 0.15; -0.02; 0.05; 0.0]}, ...
                  {[-0.8; 0.5; -0.2; 0.1; 0.2; -0.3], ...
                   [-0.7; 0.45; -0.25; 0.0; 0.1; 0.2]}
                 }
        bad_state = {zeros(5, 1), zeros(7, 1), zeros(6, 2), struct("x", zeros(6, 1))}
        bad_scalar = {[0.1; 0.2], [], struct("a", 1), {0.1}}
    end

    methods (Test)

        % (1) Input validation

        function rejectsInvalidStates(tc, bad_state)
            [x1, x2, fn] = defaultInputs();
            tc.verifyError(@() coa.orbit.tetherAccel(bad_state, x2, 0.3, 0.1, fn), ...
                           ?MException);
            tc.verifyError(@() coa.orbit.tetherAccel(x1, bad_state, 0.3, 0.1, fn), ...
                           ?MException);
        end

        function rejectsInvalidEta(tc, bad_scalar)
            [x1, x2, fn] = defaultInputs();
            tc.verifyError(@() coa.orbit.tetherAccel(x1, x2, bad_scalar, 0.1, fn), ...
                           ?MException);
        end

        function rejectsInvalidLengthDot(tc, bad_scalar)
            [x1, x2, fn] = defaultInputs();
            tc.verifyError(@() coa.orbit.tetherAccel(x1, x2, 0.3, bad_scalar, fn), ...
                           ?MException);
        end

        function rejectsNonFunctionHandle(tc)
            [x1, x2] = defaultInputs();
            tc.verifyError(@() coa.orbit.tetherAccel(x1, x2, 0.3, 0.1, 1), ?MException);
        end

        function rejectsTooFewInputs(tc)
            [x1, x2] = defaultInputs();
            tc.verifyError(@() coa.orbit.tetherAccel(x1, x2, 0.3, 0.1), ?MException);
        end

        % (2) Output handling via nargout

        function outputSizes(tc, states, eta, length_dot)
            [x1, x2, fn] = unpackInputs(states);
            a = coa.orbit.tetherAccel(x1, x2, eta, length_dot, fn);
            tc.verifySize(a, [6, 1]);

            [a, da_dp, da_dv] = coa.orbit.tetherAccel(x1, x2, eta, length_dot, fn);
            tc.verifySize(a, [6, 1]);
            tc.verifySize(da_dp, [6, 6]);
            tc.verifySize(da_dv, [6, 6]);
        end

        function outputsConsistentAcrossNargout(tc, states, eta, length_dot)
            [x1, x2, fn] = unpackInputs(states);
            a1 = coa.orbit.tetherAccel(x1, x2, eta, length_dot, fn);
            [a2, da_dp2] = coa.orbit.tetherAccel(x1, x2, eta, length_dot, fn);
            [a3, da_dp3, ~] = coa.orbit.tetherAccel(x1, x2, eta, length_dot, fn);
            tc.verifyEqual(a2, a1);
            tc.verifyEqual(a3, a1);
            tc.verifyEqual(da_dp3, da_dp2);
        end

        function rejectsMoreThanThreeOutputs(tc)
            [x1, x2, fn] = defaultInputs();
            tc.verifyError(@() requestFour(x1, x2, fn), ?MException);
        end

        % (3) Physics of the constrained accelerations

        function matchesLinearSystemSolve(tc, states, eta, length_dot)
            % Solve the 7x7 system of Sec. 2.3 directly and compare.
            [x1, x2, fn] = unpackInputs(states);
            a = coa.orbit.tetherAccel(x1, x2, eta, length_dot, fn);

            d = x2(1:3) - x1(1:3);
            d_dot = x2(4:6) - x1(4:6);
            a1_free = fn(x1);
            a2_free = fn(x2);
            A = [(1 - eta) * eye(3), zeros(3), 2 * d; ...
                 zeros(3), eta * eye(3), -2 * d; ...
                 -2 * d', 2 * d', 0];
            gap = 2 * length_dot^2 - 2 * (d_dot' * d_dot);
            rhs = [(1 - eta) * a1_free; eta * a2_free; gap];
            sol = A \ rhs;
            tc.verifyEqual(a, sol(1:6), "AbsTol", 1e-12, "RelTol", 1e-12);
        end

        function satisfiesAccelerationLevelConstraint(tc, states, eta, length_dot)
            [x1, x2, fn] = unpackInputs(states);
            a = coa.orbit.tetherAccel(x1, x2, eta, length_dot, fn);
            d = x2(1:3) - x1(1:3);
            d_dot = x2(4:6) - x1(4:6);
            tc.verifyEqual(d' * (a(4:6) - a(1:3)), length_dot^2 - d_dot' * d_dot, ...
                           "AbsTol", 1e-12);
        end

        function constraintForceIsInternal(tc, states, eta, length_dot)
            % The constraint force acts along the tether with equal and
            % opposite impulse, so the mass-weighted acceleration is unchanged.
            [x1, x2, fn] = unpackInputs(states);
            a = coa.orbit.tetherAccel(x1, x2, eta, length_dot, fn);
            free = (1 - eta) * fn(x1) + eta * fn(x2);
            tc.verifyEqual((1 - eta) * a(1:3) + eta * a(4:6), free, "AbsTol", 1e-12);
        end

        % (4) Partials vs. central finite differences of the accelerations

        function partialsMatchCentralDifference(tc, states, eta, length_dot)
            [x1, x2, fn] = unpackInputs(states);
            [~, da_dp, da_dv] = coa.orbit.tetherAccel(x1, x2, eta, length_dot, fn);

            z = [x1(1:3); x2(1:3); x1(4:6); x2(4:6)];
            accel = @(z) coa.orbit.tetherAccel([z(1:3); z(7:9)], [z(4:6); z(10:12)], ...
                                               eta, length_dot, fn);
            jacobian_fd = centralDifferenceJacobian(accel, z);

            tc.verifyEqual(da_dp, jacobian_fd(:, 1:6), "AbsTol", 1e-6, "RelTol", 1e-6);
            tc.verifyEqual(da_dv, jacobian_fd(:, 7:12), "AbsTol", 1e-6, "RelTol", 1e-6);
        end

        function partialsHoldForGenericAccelerations(tc, eta, length_dot)
            % The Jacobian derivation is valid for any a_i(x_i), not just CR3BP.
            fn = @genericAccel;
            x1 = [0.3; 0.2; 0.1; 0.01; -0.03; 0.04];
            x2 = [0.9; -0.1; 0.4; -0.02; 0.05; 0.1];
            [~, da_dp, da_dv] = coa.orbit.tetherAccel(x1, x2, eta, length_dot, fn);

            z = [x1(1:3); x2(1:3); x1(4:6); x2(4:6)];
            accel = @(z) coa.orbit.tetherAccel([z(1:3); z(7:9)], [z(4:6); z(10:12)], ...
                                               eta, length_dot, fn);
            jacobian_fd = centralDifferenceJacobian(accel, z);

            tc.verifyEqual(da_dp, jacobian_fd(:, 1:6), "AbsTol", 1e-6, "RelTol", 1e-6);
            tc.verifyEqual(da_dv, jacobian_fd(:, 7:12), "AbsTol", 1e-6, "RelTol", 1e-6);
        end

    end
end

function [x1, x2, fn] = defaultInputs()
    x1 = [0.3; 0.2; 0.1; 0.01; -0.03; 0.04];
    x2 = [0.35; 0.1; 0.15; -0.02; 0.05; 0.0];
    fn = @(x) coa.orbit.cr3bp(x, 0.0121505856);
end

function [x1, x2, fn] = unpackInputs(states)
    x1 = states{1};
    x2 = states{2};
    fn = @(x) coa.orbit.cr3bp(x, 0.0121505856);
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

function [a, da_dp, da_dv] = genericAccel(x)
    % Arbitrary smooth acceleration with closed-form partials.
    p = x(1:3);
    v = x(4:6);
    a = [p(1) * v(2)^2 + sin(p(2)); ...
         p(3) * p(1) + v(1)^3; ...
         exp(p(1)) * v(3) + p(2)^2 * v(1)];
    if nargout > 1
        da_dp = [v(2)^2, cos(p(2)), 0; ...
                 p(3), 0, p(1); ...
                 exp(p(1)) * v(3), 2 * p(2) * v(1), 0];
        da_dv = [0, 2 * p(1) * v(2), 0; ...
                 3 * v(1)^2, 0, 0; ...
                 p(2)^2, 0, exp(p(1))];
    end
end

function requestFour(x1, x2, fn)
    [~, ~, ~, ~] = coa.orbit.tetherAccel(x1, x2, 0.3, 0.1, fn);
end
