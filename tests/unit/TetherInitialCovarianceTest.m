% Copyright (c) 2026 Andrew Binder

classdef TetherInitialCovarianceTest < matlab.unittest.TestCase
    % Tests for coa.filter.tetherInitialCovariance.

    properties (Constant)
        eta = 0.3
        sigma = struct("com_position", 5, "com_velocity", 0.5, ...
                       "relative_position", 0.01, "relative_velocity", 0.002, ...
                       "payload_position", 2, "payload_velocity", 0.2)
    end

    methods (Test)

        function covarianceEqualsFactorProductAndIsSymmetric(tc)
            [p, a] = coa.filter.tetherInitialCovariance(tc.eta, tc.sigma);

            tc.verifySize(p, [18, 18]);
            tc.verifyEqual(p, a * a', "AbsTol", 1e-14);
            tc.verifyEqual(p, p', "AbsTol", 0);
            tc.verifyGreaterThan(min(eig(p)), 0);
        end

        function relativeAndComStatesHaveRequestedSigmas(tc)
            p = coa.filter.tetherInitialCovariance(tc.eta, tc.sigma);
            e = tc.eta;
            % Map the body errors back to c = (1 - eta) p1 + eta p2 and
            % r = p2 - p1 (positions), and likewise for velocities.
            eye_3 = eye(3);
            body_to_com = [(1 - e) * eye_3, e * eye_3];
            body_to_rel = [-eye_3, eye_3];
            position_block = p(1:6, 1:6);
            velocity_block = p(10:15, 10:15);

            tc.verifyEqual(body_to_com * position_block * body_to_com', ...
                           tc.sigma.com_position^2 * eye_3, "AbsTol", 1e-12);
            tc.verifyEqual(body_to_rel * position_block * body_to_rel', ...
                           tc.sigma.relative_position^2 * eye_3, "AbsTol", 1e-12);
            tc.verifyEqual(body_to_com * velocity_block * body_to_com', ...
                           tc.sigma.com_velocity^2 * eye_3, "AbsTol", 1e-12);
            tc.verifyEqual(body_to_rel * velocity_block * body_to_rel', ...
                           tc.sigma.relative_velocity^2 * eye_3, "AbsTol", 1e-12);
            % COM and relative errors are uncorrelated.
            tc.verifyEqual(body_to_com * position_block * body_to_rel', ...
                           zeros(3), "AbsTol", 1e-12);
        end

        function comErrorLeavesSeparationUnchanged(tc)
            sigma = tc.sigma;
            sigma.relative_position = 0;
            sigma.relative_velocity = 0;
            [~, a] = coa.filter.tetherInitialCovariance(tc.eta, sigma);

            error_state = a * randn(18, 1);
            tc.verifyEqual(error_state(4:6) - error_state(1:3), zeros(3, 1), ...
                           "AbsTol", 1e-12);
            tc.verifyEqual(error_state(13:15) - error_state(10:12), ...
                           zeros(3, 1), "AbsTol", 1e-12);
        end

        function payloadIsIndependentOfThePair(tc)
            p = coa.filter.tetherInitialCovariance(tc.eta, tc.sigma);

            tc.verifyEqual(p(7:9, 1:6), zeros(3, 6), "AbsTol", 0);
            tc.verifyEqual(p(7:9, 7:9), tc.sigma.payload_position^2 * eye(3), ...
                           "AbsTol", 1e-14);
            tc.verifyEqual(p(16:18, 16:18), ...
                           tc.sigma.payload_velocity^2 * eye(3), "AbsTol", 1e-14);
            tc.verifyEqual(p(1:9, 10:18), zeros(9), "AbsTol", 0);
        end

        function equalMassesWithEqualSigmasGiveExpectedBodyVariance(tc)
            sigma = struct("com_position", 3, "com_velocity", 1, ...
                           "relative_position", 4, "relative_velocity", 2, ...
                           "payload_position", 1, "payload_velocity", 1);
            p = coa.filter.tetherInitialCovariance(0.5, sigma);

            % Var(p1) = sigma_c^2 + eta^2 * sigma_r^2
            tc.verifyEqual(p(1, 1), 3^2 + 0.25 * 4^2, "AbsTol", 1e-12);
        end

        function invalidInputsErrorWithIdentifiers(tc)
            bad = tc.sigma;
            bad.com_position = -1;
            tc.verifyError(@() coa.filter.tetherInitialCovariance(tc.eta, bad), ...
                           "coa:filter:tetherInitialCovariance:InvalidSigma");

            missing = rmfield(tc.sigma, "payload_velocity");
            args = {tc.eta, missing};
            tc.verifyError(@() coa.filter.tetherInitialCovariance(args{:}), ...
                           "coa:filter:tetherInitialCovariance:MissingField");

            args = {1.2, tc.sigma};
            tc.verifyError(@() coa.filter.tetherInitialCovariance(args{:}), ...
                           "coa:filter:tetherInitialCovariance:" + ...
                           "InvalidMassFraction");
        end

    end
end
