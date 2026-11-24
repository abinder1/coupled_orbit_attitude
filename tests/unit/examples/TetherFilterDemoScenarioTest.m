% Copyright (c) 2026 Andrew Binder

classdef TetherFilterDemoScenarioTest < matlab.unittest.TestCase
    % Tests for the filter demo's separate truth and filter dynamics paths.
    methods (Test)

        function cr3bpPathsMatchLegacyRhs(tc)
            import tetherFilterDemoSupport.*
            [config, state, process_noise] = inputs();
            truth_rhs = Scenario.truthDynamicsFunction(config);
            filter_rhs = Scenario.filterDynamicsFunction(config, process_noise);
            legacy_truth_rhs = legacyRhs(config, struct());
            legacy_filter_rhs = legacyRhs(config, process_noise);

            tc.verifyEqual(truth_rhs(0.2, state), legacy_truth_rhs(0.2, state), ...
                           "AbsTol", 0);
            packed_state = [state; zeros(324, 1)];
            tc.verifyEqual(filter_rhs(0.2, packed_state), ...
                           legacy_filter_rhs(0.2, packed_state), "AbsTol", 0);
        end

        function shortIntegrationsMatchLegacyRhs(tc)
            import tetherFilterDemoSupport.*
            [config, state, process_noise] = inputs();
            options = odeset("RelTol", 1e-10, "AbsTol", 1e-12);
            legacy_truth_rhs = legacyRhs(config, struct());
            truth_rhs = Scenario.truthDynamicsFunction(config);
            [expected_time, expected_state] = ...
                ode113(legacy_truth_rhs, [0, 0.01], state, options);
            [actual_time, actual_state] = ...
                ode113(truth_rhs, [0, 0.01], state, options);
            tc.verifyEqual(actual_time, expected_time, "AbsTol", 0);
            tc.verifyEqual(actual_state, expected_state, "AbsTol", 0);

            packed_state = [state; zeros(324, 1)];
            options = odeset("RelTol", 1e-9, ...
                             "AbsTol", [1e-10 * ones(18, 1); 1e-24 * ones(324, 1)]);
            legacy_filter_rhs = legacyRhs(config, process_noise);
            filter_rhs = Scenario.filterDynamicsFunction(config, process_noise);
            [expected_time, expected_state] = ...
                ode113(legacy_filter_rhs, [0, 0.01], packed_state, options);
            [actual_time, actual_state] = ...
                ode113(filter_rhs, [0, 0.01], packed_state, options);
            tc.verifyEqual(actual_time, expected_time, "AbsTol", 0);
            tc.verifyEqual(actual_state, expected_state, "AbsTol", 0);
        end

    end
end

function [config, state, process_noise] = inputs()
    config.truth_model = "cr3bp";
    config.mu = 0.012150585;
    config.mass_2_fraction = 0.5;
    config.tether_length_rate = 0.001;
    config.spring_damper = struct("initial_length", 0.1, ...
                                  "stiffness", 3, "damping", 4);
    state = [0.6; 0.1; 0.12; 0.67; 0.04; 0.17; 0.8; 0.2; 0.1; ...
             0.02; 0.24; 0.01; 0.03; 0.2; 0.02; -0.02; 0.15; 0.1];
    process_noise = struct("tension_psd", 1e-6, ...
                           "transverse_psd", 1e-8, ...
                           "body_3_acceleration_psd", 1e-8);
end

function rhs = legacyRhs(config, process_noise)
    accel_fn = @(~, x) coa.orbit.cr3bp(x, config.mu);
    propagate = @coa.orbit.tetheredThreeBodyEkfPropagation;
    args = {config.mass_2_fraction, config.tether_length_rate, ...
            accel_fn, process_noise, config.spring_damper};
    rhs = @(t, x) propagate(t, x, args{:});
end
