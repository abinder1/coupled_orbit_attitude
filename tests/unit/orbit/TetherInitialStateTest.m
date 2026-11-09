% Copyright (c) 2026 Andrew Binder

classdef TetherInitialStateTest < matlab.unittest.TestCase
    % Tests for coa.orbit.tetherInitialState.

    methods (Test)

        function matchesPreExtractionReference(tc)
            com = [0.2; -0.1; 0.3; 0.04; 0.2; -0.03];
            state = coa.orbit.tetherInitialState(com, 0.002, 0.27, 1.3, ...
                                                 -0.003, [1; 2; 2] / 3);
            expected = [0.20018000000000002; -0.099640000000000006; ...
                        0.30036000000000002; 0.19951333333333335; ...
                        -0.10097333333333335; 0.29902666666666666; ...
                        0.039621999999999997; 0.19951400000000002; ...
                        -0.030539999999999998; 0.041022000000000003; ...
                        0.20131400000000002; -0.028539999999999999];
            tc.verifyEqual(state, expected, "AbsTol", 1e-15);
        end

        function preservesCenterOfMassAndRelativeMotion(tc)
            com = [0.2; -0.1; 0.3; 0.04; 0.2; -0.03];
            direction = [1; 2; 2] / 3;
            length_nd = 0.002;
            spin_rate = 1.3;
            for eta = [0, 0.27, 0.5, 1]
                for length_dot = [-0.003, 0, 0.005]
                    state = coa.orbit.tetherInitialState(com, length_nd, eta, ...
                                                         spin_rate, length_dot, ...
                                                         direction);
                    p1 = state(1:3);
                    p2 = state(4:6);
                    v1 = state(7:9);
                    v2 = state(10:12);
                    tc.verifyEqual((1 - eta) * p1 + eta * p2, com(1:3), ...
                                   "AbsTol", 1e-15);
                    tc.verifyEqual((1 - eta) * v1 + eta * v2, com(4:6), ...
                                   "AbsTol", 1e-15);
                    tc.verifyEqual(p1 - p2, length_nd * direction, ...
                                   "AbsTol", 1e-15);
                    expected_velocity = cross((spin_rate - 1) * [0; 0; 1], ...
                                              length_nd * direction) + ...
                                        length_dot * direction;
                    tc.verifyEqual(v1 - v2, expected_velocity, "AbsTol", 1e-15);
                end
            end
        end

        function directionIsNotNormalized(tc)
            state = coa.orbit.tetherInitialState(zeros(6, 1), 2, 0.5, 1, ...
                                                 0, [2; 0; 0]);
            tc.verifyEqual(state(1:3) - state(4:6), [4; 0; 0]);
        end

        function invalidInputsAreRejected(tc)
            args = {zeros(6, 1), 1, 0.5, 1, 0, [1; 0; 0]};
            invalid = args;
            invalid{2} = 0;
            tc.verifyError(@() coa.orbit.tetherInitialState(invalid{:}), ?MException);
            invalid = args;
            invalid{3} = 1.1;
            tc.verifyError(@() coa.orbit.tetherInitialState(invalid{:}), ?MException);
            invalid = args;
            invalid{6} = [1; 0; NaN];
            tc.verifyError(@() coa.orbit.tetherInitialState(invalid{:}), ?MException);
        end

    end
end
