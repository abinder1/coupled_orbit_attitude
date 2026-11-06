% Copyright (c) 2026 Andrew Binder

classdef MeanChiSquareBoundsTest < matlab.unittest.TestCase
    % Tests for coa.filter.meanChiSquareBounds.

    methods (Test)

        function matchesClosedFormForTwoTotalDegreesOfFreedom(tc)
            [lower, upper] = coa.filter.meanChiSquareBounds(2, 1, 0.9);
            tc.verifyEqual(lower, -log(0.95), "AbsTol", 1e-14);
            tc.verifyEqual(upper, -log(0.05), "AbsTol", 1e-14);
        end

        function preservesDegreesOfFreedomShape(tc)
            dof = [1, 2; 6, 18];
            [lower, upper] = coa.filter.meanChiSquareBounds(100, dof, 0.95);
            tc.verifySize(lower, size(dof));
            tc.verifySize(upper, size(dof));
            tc.verifyLessThan(lower, dof);
            tc.verifyGreaterThan(upper, dof);
            for k = 1:numel(dof)
                [scalar_lower, scalar_upper] = ...
                    coa.filter.meanChiSquareBounds(100, dof(k), 0.95);
                tc.verifyEqual(lower(k), scalar_lower);
                tc.verifyEqual(upper(k), scalar_upper);
            end
        end

        function invalidInputsAreRejected(tc)
            for count = {0, 1.5, Inf}
                calculate = @() coa.filter.meanChiSquareBounds(count{1}, 2, 0.95);
                tc.verifyError(calculate, ?MException);
            end
            for dof = {[], 0, [1, -2], 1.5, NaN}
                tc.verifyError(@() coa.filter.meanChiSquareBounds(1, dof{1}, 0.95), ...
                               ?MException);
            end
            for confidence = {0, 1, -0.1, NaN}
                calculate = @() coa.filter.meanChiSquareBounds(1, 2, confidence{1});
                tc.verifyError(calculate, ?MException);
            end
        end

    end
end
