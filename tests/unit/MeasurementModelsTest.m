% Copyright (c) 2026 Andrew Binder

classdef MeasurementModelsTest < matlab.unittest.TestCase
    % Tests for coa.meas measurement models and primitives.

    properties (TestParameter)
        state = {
                 [0.3; 0.2; 0.1; 0.35; 0.1; 0.15; -0.5; 0.4; 0.2; ...
                  0.1; -0.2; 0.3; 0.05; 0.4; -0.1; 0.2; 0.1; -0.3]
                }
    end

    methods (Test)

        function primitiveJacobiansMatchFiniteDifferences(tc)
            rho = [0.3; -0.4; 0.7];
            rho_2 = [-0.2; 0.8; 0.5];
            rho_dot = [0.1; 0.2; -0.3];

            [~, angles_jacobian] = coa.meas.angles(rho);
            angles_fd = centralDifferenceJacobian(@(x) coa.meas.angles(x), rho);
            tc.verifyEqual(angles_jacobian, angles_fd, ...
                           "AbsTol", 1e-7, "RelTol", 1e-7);

            [~, size_jacobian_1, size_jacobian_2] = ...
                coa.meas.angularSize(rho, rho_2);
            size_of_1 = @(x) coa.meas.angularSize(x, rho_2);
            size_of_2 = @(x) coa.meas.angularSize(rho, x);
            size_fd_1 = centralDifferenceJacobian(size_of_1, rho);
            size_fd_2 = centralDifferenceJacobian(size_of_2, rho_2);
            tc.verifyEqual(size_jacobian_1, size_fd_1, ...
                           "AbsTol", 1e-7, "RelTol", 1e-7);
            tc.verifyEqual(size_jacobian_2, size_fd_2, ...
                           "AbsTol", 1e-7, "RelTol", 1e-7);

            [~, range_jacobian] = coa.meas.range(rho);
            range_fd = centralDifferenceJacobian(@(x) coa.meas.range(x), rho);
            tc.verifyEqual(range_jacobian, range_fd, ...
                           "AbsTol", 1e-7, "RelTol", 1e-7);

            [~, rate_jacobian, rate_velocity_jacobian] = ...
                coa.meas.rangeRate(rho, rho_dot);
            rate_of_position = @(x) coa.meas.rangeRate(x, rho_dot);
            rate_of_velocity = @(x) coa.meas.rangeRate(rho, x);
            rate_fd = centralDifferenceJacobian(rate_of_position, rho);
            rate_velocity_fd = centralDifferenceJacobian(rate_of_velocity, rho_dot);
            tc.verifyEqual(rate_jacobian, rate_fd, ...
                           "AbsTol", 1e-7, "RelTol", 1e-7);
            tc.verifyEqual(rate_velocity_jacobian, rate_velocity_fd, ...
                           "AbsTol", 1e-7, "RelTol", 1e-7);

            [~, pnt_jacobian] = coa.meas.pnt(rho);
            tc.verifyEqual(pnt_jacobian, eye(3));
        end

        function primitiveValueOnlyCallsAllowJacobianSingularities(tc)
            tc.verifyEqual(coa.meas.range(zeros(3, 1)), 0);
            tc.verifyError(@() requestRangeJacobianAtOrigin(), ...
                           'coa:meas:range:Singular');

            tc.verifyEqual(coa.meas.angularSize([1; 0; 0], [2; 0; 0]), 0);
            tc.verifyEqual(coa.meas.angularSize([1; 0; 0], [-2; 0; 0]), pi);
            tc.verifyError(@() requestAngularSizeJacobianAtAlignment(), ...
                           'coa:meas:angularSize:Singular');
        end

        function allMeasurementJacobiansMatchFiniteDifferences(tc, state)
            mask = 2^8 - 1;
            covariances = coa.meas.measurementCovariances(1, 1);
            model = @coa.meas.measurementModels;
            [h, dh_dstate, noise_cov] = model(state, mask, covariances);

            tc.verifySize(h, [30, 1]);
            tc.verifySize(dh_dstate, [30, 18]);
            tc.verifySize(noise_cov, [30, 30]);

            model_of_state = @(x) coa.meas.measurementModels(x, mask);
            numeric_jacobian = centralDifferenceJacobian(model_of_state, state);
            tc.verifyEqual(dh_dstate, numeric_jacobian, ...
                           "AbsTol", 2e-6, "RelTol", 2e-6);
        end

        function nominalMeasurementCovariancesUseConfiguredSensorAccuracies(tc)
            length_scale = 384400;
            speed_scale = 2;
            covariances = coa.meas.measurementCovariances(length_scale, speed_scale);
            arcsec = pi / 180 / 3600;
            meter = 1e-3 / length_scale;
            meter_per_sec = 1e-3 / speed_scale;

            tc.verifyEqual(covariances.ground_angles, (2 * arcsec)^2 * eye(2));
            tc.verifyEqual(covariances.camera_angles, (5 * arcsec)^2 * eye(2));
            tc.verifyEqual(covariances.star_tracker_angles, (10 * arcsec)^2 * eye(2));
            tc.verifyEqual(covariances.camera_angular_size, 2 * (5 * arcsec)^2);
            tc.verifyEqual(covariances.dsn_range, (10 * meter)^2);
            tc.verifyEqual(covariances.dsn_range_rate, (1e-3 * meter_per_sec)^2);
            tc.verifyEqual(covariances.crosslink_range, (1 * meter)^2);
            tc.verifyEqual(covariances.crosslink_range_rate, (1e-3 * meter_per_sec)^2);
            tc.verifyEqual(covariances.gnss_position, (10 * meter)^2 * eye(3));
        end

        function groupMetadataMatchesAllMeasurementMasks(tc, state)
            sizes = [6, 4, 1, 3, 3, 2, 2, 9];
            names = ["Ground optical angles", "Camera / star-tracker angles", ...
                     "Camera angular size", "DSN range", "DSN range rate", ...
                     "Crosslink range", "Crosslink range rate", "GNSS position"];
            covariances = coa.meas.measurementCovariances(1, 1);
            single_h = cell(1, 8);
            single_jacobian = cell(1, 8);
            single_covariance = cell(1, 8);
            for group = 1:8
                [single_h{group}, single_jacobian{group}, ...
                 single_covariance{group}] = ...
                    coa.meas.measurementModels(state, 2^(group - 1), covariances);
            end

            for mask = 0:255
                [groups, angle_rows] = coa.meas.measurementGroups(mask);
                [h, jacobian, covariance] = ...
                    coa.meas.measurementModels(state, mask, covariances);
                active = find(bitget(mask, 1:8));
                tc.verifyEqual(reshape([groups.index], 1, []), active);
                tc.verifyEqual(angle_rows, ...
                               (1:(6 * bitget(mask, 1) + 4 * bitget(mask, 2)))');
                offset = 0;
                for k = 1:numel(active)
                    group = active(k);
                    rows = offset + (1:sizes(group));
                    tc.verifyEqual(groups(k).name, names(group));
                    tc.verifyEqual(groups(k).rows, rows);
                    tc.verifyEqual(h(rows), single_h{group});
                    tc.verifyEqual(jacobian(rows, :), single_jacobian{group});
                    tc.verifyEqual(covariance(rows, rows), single_covariance{group});
                    offset = offset + sizes(group);
                end
                tc.verifySize(h, [offset, 1]);
                tc.verifySize(jacobian, [offset, 18]);
                tc.verifySize(covariance, [offset, offset]);
            end
        end

        function invalidMeasurementGroupMasksAreRejected(tc, state)
            tc.verifyError(@() coa.meas.measurementGroups(256), ...
                           'coa:meas:measurementGroups:UnknownGroup');
            tc.verifyError(@() coa.meas.measurementGroups(-1), ?MException);
            tc.verifyError(@() coa.meas.measurementGroups(0.5), ?MException);
            tc.verifyError(@() coa.meas.measurementModels(state, 256), ...
                           'coa:meas:measurementModels:UnknownGroup');
        end

        function measurementGroupsOnlyValidateRequestedCovariances(tc, state)
            [h, dh_dstate] = coa.meas.measurementModels(state, 2^7);
            tc.verifySize(h, [9, 1]);
            tc.verifySize(dh_dstate, [9, 18]);

            gnss_cov = struct("gnss_position", eye(3));
            model = @coa.meas.measurementModels;
            [h, dh_dstate, noise_cov] = model(state, 2^7, gnss_cov);
            tc.verifySize(h, [9, 1]);
            tc.verifySize(dh_dstate, [9, 18]);
            tc.verifyEqual(noise_cov, eye(9));

            [h_empty, dh_empty, cov_empty] = ...
                coa.meas.measurementModels(state, 0);
            tc.verifyEmpty(h_empty);
            tc.verifySize(dh_empty, [0, 18]);
            tc.verifySize(cov_empty, [0, 0]);
        end

        function measurementValueOnlyCallDoesNotRequireJacobian(tc, state)
            state(1:3) = 0;
            range_measurements = coa.meas.measurementModels(state, 2^3);
            tc.verifyEqual(range_measurements(1), 0);
            tc.verifyError(@() requestRangeMeasurementJacobian(state), ...
                           'coa:meas:range:Singular');
        end

        function measurementCovarianceMustBePositiveDefinite(tc, state)
            mask = 2^7;
            bad_cov_1 = struct("gnss_position", [1, 1, 0; 0, 1, 0; 0, 0, 1]);
            bad_cov_2 = struct("gnss_position", diag([1, 1, 0]));
            tc.verifyError(@() requestMeasurementCovariance(state, mask, struct()), ...
                           'coa:meas:measurementModels:MissingCovariance');
            tc.verifyError(@() requestMeasurementCovariance(state, mask, bad_cov_1), ...
                           'coa:meas:measurementModels:InvalidCovariance');
            tc.verifyError(@() requestMeasurementCovariance(state, mask, bad_cov_2), ...
                           'coa:meas:measurementModels:InvalidCovariance');
        end

        function invalidInputsAreRejected(tc, state)
            tc.verifyError(@() coa.meas.angles([1; 2; NaN]), ?MException);
            tc.verifyError(@() coa.meas.measurementModels(zeros(17, 1), 0), ...
                           ?MException);
            tc.verifyError(@() coa.meas.measurementModels(state, 1.5), ...
                           ?MException);
            tc.verifyError(@() coa.meas.measurementCovariances(0, 1), ...
                           ?MException);
        end

    end
end

function jacobian = centralDifferenceJacobian(f, x)
    output = f(x);
    jacobian = zeros(numel(output), numel(x));
    for k = 1:numel(x)
        step = 1e-6 * max(1, abs(x(k)));
        perturbation = zeros(size(x));
        perturbation(k) = step;
        jacobian(:, k) = (f(x + perturbation) - f(x - perturbation)) / ...
                         (2 * step);
    end
end

function requestRangeJacobianAtOrigin()
    [~, ~] = coa.meas.range(zeros(3, 1));
end

function requestAngularSizeJacobianAtAlignment()
    [~, ~, ~] = coa.meas.angularSize([1; 0; 0], [2; 0; 0]);
end

function requestMeasurementCovariance(state, mask, meas_covs)
    [~, ~, ~] = coa.meas.measurementModels(state, mask, meas_covs);
end

function requestRangeMeasurementJacobian(state)
    [~, ~] = coa.meas.measurementModels(state, 2^3);
end
