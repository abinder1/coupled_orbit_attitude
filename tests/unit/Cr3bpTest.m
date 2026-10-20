% Copyright (c) 2026 Andrew Binder

classdef Cr3bpTest < matlab.unittest.TestCase
    % Tests for coa.+orbit.cr3bp

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

        function singleOutputIsStateDerivative(tc, state, mu)
            xdot = coa.orbit.cr3bp(state, mu);
            tc.verifySize(xdot, [6, 1]);
            tc.verifyEqual(xdot(1:3), state(4:6));
        end

        function secondOutputIsJacobian(tc, state, mu)
            [xdot, J] = coa.orbit.cr3bp(state, mu);
            tc.verifySize(xdot, [6, 1]);
            tc.verifySize(J, [6, 6]);
        end

        function stateDerivativeIndependentOfNargout(tc, state, mu)
            xdot1 = coa.orbit.cr3bp(state, mu);
            [xdot2, ~] = coa.orbit.cr3bp(state, mu);
            tc.verifyEqual(xdot2, xdot1);
        end

        function rejectsMoreThanTwoOutputs(tc)
            % Requesting two outputs is valid; three is not.
            tc.verifyError(@() requestThree(), ?MException);
        end

        function jacobianHasKinematicBlocks(tc, state, mu)
            [~, J] = coa.orbit.cr3bp(state, mu);
            tc.verifyEqual(J(1:3, 1:3), zeros(3));
            tc.verifyEqual(J(1:3, 4:6), eye(3));
            % Coriolis block: -2*skew([0; 0; 1])
            tc.verifyEqual(J(4:6, 4:6), -2 * coa.util.skew([0; 0; 1]));
        end

        % (3) Jacobian vs. central finite differences of the state derivative

        function jacobianMatchesCentralDifference(tc, state, mu)
            [~, J] = coa.orbit.cr3bp(state, mu);
            jacobian_fd = centralDifferenceJacobian(@(z) coa.orbit.cr3bp(z, mu), state);
            tc.verifyEqual(J, jacobian_fd, "AbsTol", 1e-6, "RelTol", 1e-6);
        end

        function accelerationJacobianMatchesCentralDifference(tc, state, mu)
            [~, J] = coa.orbit.cr3bp(state, mu);
            accel = @(z) selectRows(coa.orbit.cr3bp(z, mu), 4:6);
            jacobian_fd = centralDifferenceJacobian(accel, state);
            tc.verifyEqual(J(4:6, :), jacobian_fd, "AbsTol", 1e-6, "RelTol", 1e-6);
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

function y = selectRows(x, idx)
    y = x(idx);
end

function requestThree()
    [~, ~, ~] = coa.orbit.cr3bp(zeros(6, 1) + 0.5, 0.01);
end
