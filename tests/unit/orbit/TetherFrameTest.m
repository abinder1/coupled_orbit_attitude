% Copyright (c) 2026 Andrew Binder

classdef TetherFrameTest < matlab.unittest.TestCase
    % Tests for coa.orbit.tetherFrame.

    methods (Test)

        function coordinateAxesMatchExpectedDirection(tc)
            tc.verifyEqual(coa.orbit.tetherFrame([2; 0; 0], [1; 3; 0]), eye(3));
            tc.verifyEqual(coa.orbit.tetherFrame([2; 0; 0], [1; -3; 0]), ...
                           diag([1, -1, -1]));
        end

        function obliqueFrameIsRightHandedAndOrthonormal(tc)
            position = [1; 2; 2];
            velocity = [-2; 1; 0.5];
            frame = coa.orbit.tetherFrame(position, velocity);
            normal = cross(position, velocity);
            tc.verifyEqual(frame(1, :)', position / norm(position), "AbsTol", 1e-15);
            tc.verifyEqual(frame(3, :)', normal / norm(normal), "AbsTol", 1e-15);
            tc.verifyEqual(frame * frame', eye(3), "AbsTol", 1e-14);
            tc.verifyEqual(det(frame), 1, "AbsTol", 1e-14);
        end

        function nearParallelVectorsAreNotRejected(tc)
            frame = coa.orbit.tetherFrame([1; 0; 0], [1; 1e-14; 0]);
            tc.verifyEqual(frame, eye(3));
        end

        function singularAndInvalidInputsAreRejected(tc)
            tc.verifyError(@() coa.orbit.tetherFrame(zeros(3, 1), [1; 0; 0]), ...
                           'coa:orbit:tetherFrame:Singular');
            tc.verifyError(@() coa.orbit.tetherFrame([1; 2; 3], [2; 4; 6]), ...
                           'coa:orbit:tetherFrame:Singular');
            tc.verifyError(@() coa.orbit.tetherFrame([1; 0; 0], zeros(3, 1)), ...
                           'coa:orbit:tetherFrame:Singular');
            tc.verifyError(@() coa.orbit.tetherFrame([1; NaN; 0], [0; 1; 0]), ...
                           ?MException);
        end

    end
end
