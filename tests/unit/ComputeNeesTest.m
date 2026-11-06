% Copyright (c) 2026 Andrew Binder

classdef ComputeNeesTest < matlab.unittest.TestCase
    % Tests for coa.filter.computeNees.

    methods (Test)

        function usesEachRunAndUpdatesCovariance(tc)
            errors = zeros(2, 2, 2);
            errors(:, :, 1) = [2, 4; 2, 4];
            errors(:, :, 2) = [3, 6; 3, 6];
            covariances = zeros(2, 2, 2, 2);
            covariances(:, :, 1, 1) = diag([4, 9]);
            covariances(:, :, 2, 1) = diag([4, 9]);
            covariances(:, :, 1, 2) = diag([1, 9]);
            covariances(:, :, 2, 2) = diag([16, 36]);
            blocks = struct('rows', {1:2, 1, 2});

            actual = coa.filter.computeNees(errors, covariances, blocks);

            tc.verifySize(actual, [2, 2, 3]);
            tc.verifyEqual(actual(:, :, 1), [2, 8; 5, 2], "AbsTol", 1e-14);
            tc.verifyEqual(actual(:, :, 2), [1, 4; 4, 1], "AbsTol", 1e-14);
            tc.verifyEqual(actual(:, :, 3), [1, 4; 1, 1], "AbsTol", 1e-14);
        end

        function handlesSingletonDimensionsAndCorrelations(tc)
            errors = reshape([2, 3], [1, 1, 2]);
            covariance = [4, 1; 1, 9];
            blocks = struct('rows', {[1, 2], [2, 1], 2});
            actual = coa.filter.computeNees(errors, covariance, blocks);

            tc.verifyEqual(reshape(actual, 1, []), [60 / 35, 60 / 35, 1], ...
                           "AbsTol", 1e-14);
            tc.verifyEqual(coa.filter.computeNees(2, 4, struct('rows', 1)), 1);
        end

        function invalidShapesAndBlocksAreRejected(tc)
            errors = zeros(2, 3, 4);
            covariance = repmat(eye(4), [1, 1, 3, 2]);
            valid = struct('rows', 1:4);
            tc.verifyError(@() coa.filter.computeNees(errors, eye(4), valid), ...
                           'coa:filter:computeNees:SizeMismatch');
            tc.verifyError(@() coa.filter.computeNees(errors, covariance, struct()), ...
                           'coa:filter:computeNees:InvalidBlocks');
            for rows = {[], 0, 5, 1.5, NaN, [1, 2; 3, 4]}
                block = struct('rows', rows{1});
                calculate = @() coa.filter.computeNees(errors, covariance, block);
                tc.verifyError(calculate, ...
                               'coa:filter:computeNees:InvalidBlocks');
            end
        end

    end
end
