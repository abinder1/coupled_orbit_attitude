% Copyright (c) 2026 Andrew Binder

classdef Cr3bpTest < matlab.unittest.TestCase
    % Tests for coa.orbit.cr3bp

    properties (TestParameter)
        mu = {0.0121505856, 0.3, 1e-6}
        state = {
                 [0.3; 0.2; 0.1; 0.01; -0.03; 0.04], ...
                 [-0.8; 0.5; -0.2; 0.1; 0.2; -0.3], ...
                 [1.1; -0.4; 0.0; 0.0; 0.0; 0.0]}
        bad_state = {
                     zeros(5, 1), ...
                     zeros(7, 1), ...
                     zeros(6, 2), ...
                     struct("x", zeros(6, 1)), ...
                     {zeros(6, 1)}
                    }
        bad_mu = {[0.1; 0.2], [0.1, 0.2], [], struct("mu", 0.1), {0.1}}
    end

    methods (Test)

        % (1) Input validation

        function rejectsInvalidState(tc, bad_state)
            tc.verifyError(@() coa.orbit.cr3bp(bad_state, 0.01), ?MException);
        end

        function rejectsInvalidMu(tc, bad_mu)
            x = [0.3; 0.2; 0.1; 0.01; -0.03; 0.04];
            tc.verifyError(@() coa.orbit.cr3bp(x, bad_mu), ?MException);
        end

        function rejectsTooFewInputs(tc)
            tc.verifyError(@() coa.orbit.cr3bp([0.3; 0.2; 0.1; 0; 0; 0]), ?MException);
        end

        function rejectsTooManyInputs(tc)
            tc.verifyError(@() coa.orbit.cr3bp(zeros(6, 1), 0.01, 1), ?MException);
        end

        % (2) Output handling via nargout

        function singleOutputIsAcceleration(tc, state, mu)
            a = coa.orbit.cr3bp(state, mu);
            tc.verifySize(a, [3, 1]);
        end

        function derivativeOutputsAreThreeByThree(tc, state, mu)
            [a, da_dp, da_dv] = coa.orbit.cr3bp(state, mu);
            tc.verifySize(a, [3, 1]);
            tc.verifySize(da_dp, [3, 3]);
            tc.verifySize(da_dv, [3, 3]);
        end

        function twoOutputsMatchThreeOutputs(tc, state, mu)
            [a2, da_dp2] = coa.orbit.cr3bp(state, mu);
            [a3, da_dp3, ~] = coa.orbit.cr3bp(state, mu);
            tc.verifyEqual(a2, a3);
            tc.verifyEqual(da_dp2, da_dp3);
        end

        function accelerationIndependentOfNargout(tc, state, mu)
            a1 = coa.orbit.cr3bp(state, mu);
            [a2, ~, ~] = coa.orbit.cr3bp(state, mu);
            tc.verifyEqual(a2, a1);
        end

        function rejectsMoreThanThreeOutputs(tc)
            % Requesting three outputs is valid; four is not.
            tc.verifyError(@() requestFour(), ?MException);
        end

        function velocityPartialIsCoriolisTerm(tc, state, mu)
            [~, ~, da_dv] = coa.orbit.cr3bp(state, mu);
            % Coriolis term: -2*skew([0; 0; 1])
            tc.verifyEqual(da_dv, -2 * coa.util.skew([0; 0; 1]));
        end

        % (3) Partials vs. central finite differences of the acceleration

        function positionPartialMatchesCentralDifference(tc, state, mu)
            [~, da_dp] = coa.orbit.cr3bp(state, mu);
            accel = @(r) coa.orbit.cr3bp([r; state(4:6)], mu);
            da_dp_fd = centralDifferenceJacobian(accel, state(1:3));
            tc.verifyEqual(da_dp, da_dp_fd, "AbsTol", 1e-6, "RelTol", 1e-6);
        end

        function velocityPartialMatchesCentralDifference(tc, state, mu)
            [~, ~, da_dv] = coa.orbit.cr3bp(state, mu);
            accel = @(v) coa.orbit.cr3bp([state(1:3); v], mu);
            da_dv_fd = centralDifferenceJacobian(accel, state(4:6));
            tc.verifyEqual(da_dv, da_dv_fd, "AbsTol", 1e-6, "RelTol", 1e-6);
        end

    end
end

function jacobian_fd = centralDifferenceJacobian(f, x)
    n = numel(x);
    f0 = f(x);
    jacobian_fd = zeros(numel(f0), n);
    for k = 1:n
        h = 1e-6 * max(1, abs(x(k)));
        e = zeros(n, 1);
        e(k) = h;
        jacobian_fd(:, k) = (f(x + e) - f(x - e)) / (2 * h);
    end
end

function requestFour()
    [~, ~, ~, ~] = coa.orbit.cr3bp(zeros(6, 1) + 0.5, 0.01);
end
