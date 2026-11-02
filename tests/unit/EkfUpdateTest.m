% Copyright (c) 2026 Andrew Binder

classdef EkfUpdateTest < matlab.unittest.TestCase
    % Tests for coa.filter.ekfUpdate.

    properties (Constant)
        state = [0.3; 0.2; 0.1; 0.35; 0.1; 0.15; -0.5; 0.4; 0.2; ...
                 0.1; -0.2; 0.3; 0.05; 0.4; -0.1; 0.2; 0.1; -0.3]
    end

    methods (Test)

        function scalarLinearCaseMatchesClosedForm(tc)
            prior_variance = 4;
            noise_variance = 1;
            model = @(x) affineModel(x, @(y) y, 1, noise_variance);

            [x_post, p_post, info] = coa.filter.ekfUpdate(2, prior_variance, ...
                                                          3, model);

            gain = prior_variance / (prior_variance + noise_variance);
            tc.verifyEqual(x_post, 2 + gain * 1, "AbsTol", 1e-14);
            tc.verifyEqual(p_post, (1 - gain) * prior_variance, "AbsTol", 1e-14);
            tc.verifyEqual(info.innovation, 1);
            tc.verifyEqual(info.innovation_covariance, 5);
            tc.verifyEqual(info.gain, gain, "AbsTol", 1e-14);
            tc.verifyEqual(info.residual, 3 - x_post, "AbsTol", 1e-14);
            % For the optimal gain, R - H*P+*H' = R*(1 - H*K) variance.
            tc.verifyEqual(info.residual_covariance, ...
                           noise_variance - p_post, "AbsTol", 1e-14);
        end

        function matchesInlineJosephUpdateForMeasurementModels(tc)
            mask = 2^8 - 1;
            covs = coa.meas.measurementCovariances(1, 1);
            prior_covariance = diag(1e-4 * (1:18));
            truth = tc.state + 1e-3 * (1:18)';
            prior = tc.state;

            model = @(x) coa.meas.measurementModels(x, mask, covs);
            observation = model(truth);
            angle_rows = (1:10)';

            args = {prior, prior_covariance, observation, model, angle_rows};
            [x_post, p_post, info] = coa.filter.ekfUpdate(args{:});

            [h, H, R] = model(prior);
            nu = observation - h;
            nu(angle_rows) = mod(nu(angle_rows) + pi, 2 * pi) - pi;
            S = H * prior_covariance * H' + R;
            K = (S \ (H * prior_covariance))';
            A = eye(18) - K * H;
            expected_covariance = A * prior_covariance * A' + K * R * K';
            expected_covariance = (expected_covariance + expected_covariance') / 2;

            tc.verifyEqual(x_post, prior + K * nu, "RelTol", 1e-10, "AbsTol", 1e-12);
            tc.verifyEqual(p_post, expected_covariance, ...
                           "RelTol", 1e-10, "AbsTol", 1e-14);
            tc.verifyEqual(info.innovation, nu, "AbsTol", 1e-12);
            tc.verifyEqual(info.gain, K, "RelTol", 1e-10, "AbsTol", 1e-12);
        end

        function posteriorCovarianceIsSymmetricPositiveDefinite(tc)
            model = @(x) affineModel(x, @(y) y(1:2) + y(2:3), [1, 1, 0; 0, 1, 1], ...
                                     0.1 * eye(2));
            prior_covariance = [2, 0.5, 0; 0.5, 1, 0.2; 0, 0.2, 3];

            [~, p_post] = coa.filter.ekfUpdate([1; 2; 3], prior_covariance, ...
                                               [3; 5], model);

            tc.verifyEqual(p_post, p_post');
            tc.verifyGreaterThan(min(eig(p_post)), 0);
            tc.verifyLessThan(trace(p_post), trace(prior_covariance));
        end

        function anglesAreWrappedAcrossBranchCut(tc)
            model = @(x) affineModel(x, @(y) y, 1, 1e-4);
            prior_angle = pi - 0.01;
            observation = -pi + 0.01;

            args = {prior_angle, 1e-2, observation, model};
            [x_wrapped, ~, wrapped] = coa.filter.ekfUpdate(args{:}, 1);
            [~, ~, unwrapped] = coa.filter.ekfUpdate(args{:});

            tc.verifyEqual(wrapped.innovation, 0.02, "AbsTol", 1e-12);
            tc.verifyLessThan(abs(wrapped.residual), pi);
            tc.verifyEqual(unwrapped.innovation, observation - prior_angle, ...
                           "AbsTol", 1e-12);
            tc.verifyLessThan(abs(x_wrapped - prior_angle), 0.05);
        end

        function emptyMeasurementReturnsPrior(tc)
            model = @(x) affineModel(x, @(~) zeros(0, 1), zeros(0, 2), zeros(0, 0));
            prior_covariance = [1, 0.1; 0.1, 2];

            args = {[1; 2], prior_covariance, zeros(0, 1), model};
            [x_post, p_post, info] = coa.filter.ekfUpdate(args{:});

            tc.verifyEqual(x_post, [1; 2]);
            tc.verifyEqual(p_post, prior_covariance);
            tc.verifyEmpty(info.innovation);
        end

        function singularInnovationCovarianceErrors(tc)
            model = @(x) affineModel(x, @(y) y(1), [1, 0], 0);
            tc.verifyError(@() coa.filter.ekfUpdate([1; 2], zeros(2), 1, model), ...
                           'coa:filter:ekfUpdate:SingularInnovationCovariance');
        end

        function inconsistentSizesError(tc)
            model = @(x) affineModel(x, @(y) y(1), [1, 0], 1);

            tc.verifyError(@() coa.filter.ekfUpdate([1; 2], eye(3), 1, model), ...
                           'coa:filter:ekfUpdate:SizeMismatch');
            tc.verifyError(@() coa.filter.ekfUpdate([1; 2], eye(2), [1; 2], model), ...
                           'coa:filter:ekfUpdate:SizeMismatch');
            tc.verifyError(@() coa.filter.ekfUpdate([1; 2], eye(2), 1, model, 2), ...
                           'coa:filter:ekfUpdate:SizeMismatch');
        end

        function nonfiniteInputsError(tc)
            model = @(x) affineModel(x, @(y) y(1), [1, 0], 1);

            tc.verifyError(@() coa.filter.ekfUpdate([1; 2], eye(2), NaN, model), ...
                           'MATLAB:validators:mustBeFinite');
        end

    end
end

function [h, jacobian, noise_covariance] = affineModel(x, h_fn, jacobian, ...
                                                       noise_covariance)
    % Test measurement model returning every output regardless of nargout.
    h = h_fn(x);
end
