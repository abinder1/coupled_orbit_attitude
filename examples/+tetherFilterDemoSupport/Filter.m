% Copyright (c) 2026 Andrew Binder

classdef Filter
    % Reference EKF, measurement updates, and prediction history.
    methods (Static)

        function filter_run = runFilter(scenario, initial_state, ...
                                        initial_covariance, config)
            % Propagate to each measurement, update, then continue from the posterior.
            import tetherFilterDemoSupport.*
            filter_rhs = Scenario.dynamicsFunction(config, config.process_noise);
            abs_tolerance = [1e-10 * ones(18, 1); 1e-24 * ones(324, 1)];
            filter_options = odeset('RelTol', 1e-9, 'AbsTol', abs_tolerance);

            length_scale = config.length_scale;
            velocity_scale = length_scale / config.time_scale;
            covariance_fn = @coa.meas.measurementCovariances;
            measurement_covariances = covariance_fn(length_scale, velocity_scale);
            packed_state = [initial_state; initial_covariance(:)];
            previous_time = 0;
            filter_run.history = Filter.emptyHistory();

            for k = 1:config.measurement_count
                measurement_time = config.measurement_times(k);
                filter_solution = ode113(filter_rhs, [previous_time, ...
                                                      measurement_time], ...
                                         packed_state, filter_options);
                sample_interval = linspace(previous_time, measurement_time, ...
                                           config.plot_samples_per_interval);
                filter_samples = deval(filter_solution, sample_interval);
                truth_samples = deval(scenario.truth_solution, sample_interval);
                first_sample = 1 + (k > 1);
                samples = first_sample:numel(sample_interval);
                filter_run.history = ...
                    Filter.appendPredictionHistory(filter_run.history, ...
                                                   sample_interval(samples), ...
                                                   filter_samples(:, samples), ...
                                                   truth_samples(:, samples), config);

                prior_state = filter_solution.y(1:18, end);
                prior_covariance = ...
                    Filter.unpackCovariance(filter_solution.y(19:end, end));
                truth_state = deval(scenario.truth_solution, measurement_time);
                update = Filter.measurementUpdate(prior_state, prior_covariance, ...
                                                  truth_state, ...
                                                  measurement_covariances, config);
                filter_run.updates(k, 1) = update;

                packed_state = [update.posterior_state; ...
                                update.posterior_covariance(:)];
                filter_run.history = ...
                    Filter.appendPredictionHistory(filter_run.history, ...
                                                   measurement_time, packed_state, ...
                                                   truth_state, config);
                previous_time = measurement_time;
            end

            % Predict (no update) from the last measurement to the end of the scenario.
            final_packed = packed_state;
            if previous_time < config.duration
                filter_solution = ode113(filter_rhs, [previous_time, ...
                                                      config.duration], ...
                                         packed_state, filter_options);
                sample_interval = linspace(previous_time, config.duration, ...
                                           config.plot_samples_per_interval);
                filter_samples = deval(filter_solution, sample_interval);
                truth_samples = deval(scenario.truth_solution, sample_interval);
                samples = 2:numel(sample_interval);
                filter_run.history = ...
                    Filter.appendPredictionHistory(filter_run.history, ...
                                                   sample_interval(samples), ...
                                                   filter_samples(:, samples), ...
                                                   truth_samples(:, samples), config);
                final_packed = filter_solution.y(:, end);
            end
            filter_run.final_state = final_packed(1:18);
            filter_run.final_covariance = Filter.unpackCovariance(final_packed(19:end));
        end

        function update = measurementUpdate(prior_state, prior_covariance, ...
                                            truth_state, ...
                                            measurement_covariances, config)
            % Simulate one noisy observation from truth, then apply the shared update.
            mask = config.measurement_mask;
            model = @(state) coa.meas.measurementModels(state, mask, ...
                                                        measurement_covariances);
            [true_measurement, ~, noise_covariance] = model(truth_state);
            noise_chol = chol(noise_covariance, 'lower');
            observation = true_measurement + noise_chol * randn(size(true_measurement));
            [~, angle_rows] = coa.meas.measurementGroups(mask);

            [posterior_state, posterior_covariance, info] = ...
                coa.filter.ekfUpdate(prior_state, prior_covariance, observation, ...
                                     model, angle_rows, config.ekf_iterations);

            update = info;
            update.truth_state = truth_state;
            update.prior_state = prior_state;
            update.posterior_state = posterior_state;
            update.posterior_covariance = posterior_covariance;
            update.noise_chol = noise_chol;
        end

        function history = emptyHistory()
            % Allocate arrays; update times appear twice for prior and posterior.
            history.time = zeros(1, 0);
            history.position_error = zeros(9, 0);
            history.velocity_error = zeros(9, 0);
            history.position_sigma = zeros(9, 0);
            history.velocity_sigma = zeros(9, 0);
        end

        function history = appendPredictionHistory(history, time, packed, truth, config)
            % Append signed Cartesian errors and componentwise standard deviations.
            state_error = packed(1:18, :) - truth;
            variance = packed(18 + 1 + 19 * (0:17), :);
            position_scale = config.length_scale;
            velocity_scale = config.length_scale / config.time_scale * 1000;
            history.time = [history.time, time];
            position_error = state_error(config.position_rows, :) * position_scale;
            velocity_error = state_error(config.velocity_rows, :) * velocity_scale;
            history.position_error = [history.position_error, position_error];
            history.velocity_error = [history.velocity_error, velocity_error];
            position_sigma = sqrt(max(variance(config.position_rows, :), 0));
            velocity_sigma = sqrt(max(variance(config.velocity_rows, :), 0));
            history.position_sigma = [history.position_sigma, ...
                                      position_sigma * position_scale];
            history.velocity_sigma = [history.velocity_sigma, ...
                                      velocity_sigma * velocity_scale];
        end

        function covariance = unpackCovariance(packed_covariance)
            % Reshape a column-stacked covariance and remove integration asymmetry.
            covariance = reshape(packed_covariance, [18, 18]);
            covariance = (covariance + covariance') / 2;
        end

    end
end
