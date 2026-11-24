% Copyright (c) 2026 Andrew Binder

classdef Scenario
    % Initial conditions and truth dynamics for the tether scenario.
    methods (Static)

        function [state, covariance, factor] = makeInitialEstimate(config)
            % Initialize the tethered bodies, rendezvous payload, and covariance.
            import tetherFilterDemoSupport.*
            accel_fn = @(x) coa.orbit.cr3bp(x, config.mu);

            tether_length = config.tether_length_km / config.length_scale;
            spin_rate = config.spin_rate * config.time_scale;
            tether_rate = config.tether_length_rate;
            tether_state = coa.orbit.tetherInitialState(config.halo_state, ...
                                                        tether_length, ...
                                                        config.mass_2_fraction, ...
                                                        spin_rate, ...
                                                        tether_rate, ...
                                                        config.tether_direction);
            state = Scenario.rendezvousInitialState(tether_state, config, accel_fn);

            sigma = config.initial_sigma;
            for field = ["com_position", "relative_position", "payload_position"]
                sigma.(field) = sigma.(field) / config.length_scale;
            end
            for field = ["com_velocity", "relative_velocity", "payload_velocity"]
                sigma.(field) = sigma.(field) / 1000 * ...
                                config.time_scale / config.length_scale;
            end
            [covariance, factor] = ...
                coa.filter.tetherInitialCovariance(config.mass_2_fraction, sigma);
        end

        function state = rendezvousInitialState(tether_state, config, accel_fn)
            % Back-propagate the catch mechanism to initialize a rendezvous payload.
            import tetherFilterDemoSupport.*
            duration = config.duration;
            tether_rhs = Scenario.truthDynamicsFunction(config);

            % Body 3 is a placeholder here; it does not affect the tether dynamics.
            placeholder_state = [tether_state(1:6); tether_state(4:6); ...
                                 tether_state(7:12); tether_state(10:12)];
            tether_solution = ode113(tether_rhs, [0, duration], placeholder_state, ...
                                     odeset('RelTol', 1e-10, 'AbsTol', 1e-12));
            final_tether_state = tether_solution.y(:, end);

            catch_state = [final_tether_state(4:6); final_tether_state(13:15)];
            payload_state = Scenario.propagateSingleBody(catch_state, -duration, ...
                                                         accel_fn);

            % Reorder [p1; p2; v1; v2; p3; v3] as [p1; p2; p3; v1; v2; v3].
            combined_state = [tether_state; payload_state];
            state = [combined_state(1:6); combined_state(13:15); ...
                     combined_state(7:12); combined_state(16:18)];
        end

        function scenario = propagateScenario(initial_state, config)
            % Propagate the nominal initial condition; every trial shares this truth.
            import tetherFilterDemoSupport.*
            scenario.truth_solution = ...
                ode113(Scenario.truthDynamicsFunction(config), ...
                       [0, config.duration], initial_state, ...
                       odeset('RelTol', 1e-10, 'AbsTol', 1e-12));
        end

        function rhs = truthDynamicsFunction(config)
            % Select the truth dynamics model for the filter scenario.
            import tetherFilterDemoSupport.*
            switch config.truth_model
                case "cr3bp"
                    rhs = Scenario.cr3bpDynamicsFunction(config, struct());
                otherwise
                    error('tetherFilterDemoSupport:Scenario:UnknownTruthModel', ...
                          'Unknown truth model: %s.', config.truth_model);
            end
        end

        function rhs = filterDynamicsFunction(config, process_noise)
            % The filter dynamics remain the CR3BP model with covariance propagation.
            import tetherFilterDemoSupport.*
            rhs = Scenario.cr3bpDynamicsFunction(config, process_noise);
        end

        function state = perturbInitialState(initial_state, initial_factor)
            % Draw a zero-mean initial error with the covariance used by the filter.
            state = initial_state + initial_factor * randn(18, 1);
        end

        function state = propagateSingleBody(state_0, duration, accel_fn)
            % Coast a single body [p; v]; negative duration propagates backward.
            rhs = @(~, x) [x(4:6); accel_fn(x)];
            solution = ode113(rhs, [0, duration], state_0, ...
                              odeset('AbsTol', 1e-10, 'RelTol', 1e-10));
            state = solution.y(:, end);
        end

    end

    methods (Static)

        function rhs = cr3bpDynamicsFunction(config, process_noise)
            % Construct the CR3BP tether RHS for truth or filter propagation.
            accel_fn = @(x) coa.orbit.cr3bp(x, config.mu);
            args = {config.mass_2_fraction, config.tether_length_rate, accel_fn, ...
                    process_noise, config.spring_damper};
            rhs = @(t, x) coa.orbit.tetheredThreeBodyEkfPropagation(t, x, args{:});
        end

    end
end
