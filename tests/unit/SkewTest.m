classdef SkewTest < matlab.unittest.TestCase
    % Tests for coa.util.skew (eq:skew)

    properties (TestParameter)
        v = {[1;2;3], [0;0;0], [-4.2;0.1;7]}
    end

    methods (Test)
        function matchesCrossProduct(tc, v)
            w = [0.3; -1.1; 2.5];
            tc.verifyEqual(coa.util.skew(v)*w, cross(v, w), "AbsTol", 1e-12);
        end

        function isSkewSymmetric(tc, v)
            S = coa.util.skew(v);
            tc.verifyEqual(S.', -S);
        end
    end
end
