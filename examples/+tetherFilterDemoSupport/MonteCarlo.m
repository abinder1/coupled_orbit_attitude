% Copyright (c) 2026 Andrew Binder

classdef MonteCarlo
    % Repeated filter runs and the results shared by their analyses.
    methods (Static)

        function monte_carlo = runMonteCarlo(scenario, initial_state, ...
                                             initial_covariance, ...
                                             initial_factor, filter_run, config)
            % Repeat the full EKF against the fixed truth, with a fresh initial estimate
            % sample and measurement noise per run. Every run integrates its
            % own covariance and computes its own gains, so ANEES and NIS test the
            % filter as implemented.
            import tetherFilterDemoSupport.*
            filter_rhs = Scenario.dynamicsFunction(config, config.process_noise);
            abs_tolerance = [1e-10 * ones(18, 1); 1e-24 * ones(324, 1)];
            filter_options = odeset('RelTol', 1e-9, 'AbsTol', abs_tolerance);
            velocity_scale = config.length_scale / config.time_scale;
            measurement_covariances = ...
                coa.meas.measurementCovariances(config.length_scale, velocity_scale);

            blocks = MonteCarlo.neesBlocks();
            updates = filter_run.updates;
            run_count = config.monte_carlo_runs;
            step_count = config.measurement_count;
            errors = zeros(run_count, step_count, 18);
            measurement_count = numel(updates(1).innovation);
            innovations = zeros(run_count, step_count, measurement_count);
            residuals = zeros(size(innovations));
            run_covariances = zeros(18, 18, step_count, run_count);
            nis = zeros(run_count, step_count);
            final_errors = zeros(run_count, 18);

            timer = tic;
            for run = 1:run_count
                truth_solution = scenario.truth_solution;
                sampled_state = Scenario.perturbInitialState(initial_state, ...
                                                             initial_factor);
                packed_state = [sampled_state; initial_covariance(:)];
                previous_time = 0;

                for k = 1:step_count
                    time = config.measurement_times(k);
                    solution = ode113(filter_rhs, [previous_time, time], ...
                                      packed_state, ...
                                      filter_options);
                    truth_state = deval(truth_solution, time);
                    prior_covariance = Filter.unpackCovariance(solution.y(19:end, end));
                    update = Filter.measurementUpdate(solution.y(1:18, end), ...
                                                      prior_covariance, truth_state, ...
                                                      measurement_covariances, config);
                    packed_state = [update.posterior_state; ...
                                    update.posterior_covariance(:)];
                    previous_time = time;

                    errors(run, k, :) = truth_state - update.posterior_state;
                    innovations(run, k, :) = update.innovation;
                    residuals(run, k, :) = update.residual;
                    run_covariances(:, :, k, run) = update.posterior_covariance;
                    nis(run, k) = update.innovation' * ...
                        (update.innovation_covariance \ update.innovation);
                end

                % Coast (no update) to the end of the scenario.
                solution = ode113(filter_rhs, [previous_time, config.duration], ...
                                  packed_state, filter_options);
                final_errors(run, :) = truth_solution.y(:, end) - solution.y(1:18, end);

                if run == run_count || floor(10 * run / run_count) > ...
                        floor(10 * (run - 1) / run_count)
                    elapsed = toc(timer);
                    fprintf(['  Monte Carlo: %d/%d runs (%3.0f%%), %.0f s elapsed, ' ...
                             '~%.0f s remaining\n'], run, run_count, ...
                            100 * run / run_count, elapsed, ...
                            elapsed * (run_count - run) / run);
                end
            end

            % Every analysis derives from these stored arrays and reference-run
            % covariances, so new metrics and plots never need further propagation.
            monte_carlo.errors = errors;           % [run, update, state]: truth - est.
            monte_carlo.innovations = innovations; % [run, update, measurement row]
            monte_carlo.residuals = residuals;            % post-fit, same shape
            monte_carlo.reference_updates = updates;      % reference run shown in plots
            monte_carlo.covariances = cat(3, updates.posterior_covariance);
            monte_carlo.blocks = blocks;
            monte_carlo.nees = coa.filter.computeNees(errors, run_covariances, blocks);
            monte_carlo.anees = squeeze(mean(monte_carlo.nees, 1));
            monte_carlo.nis = nis;                        % [run, update]
            monte_carlo.anis = mean(nis, 1)';
            monte_carlo.nis_dof = measurement_count;
            monte_carlo.time_minutes = ...
                config.measurement_times * config.time_scale / 60;
            monte_carlo.run_count = run_count;
            % Final-time snapshots: the last measurement update, then the
            % prediction-only end of the scenario. Corner-type analyses loop over these.
            to_minutes = config.time_scale / 60;
            update_endpoint = ...
                MonteCarlo.makeEndpoint('final update', ...
                                        config.measurement_times(end) * to_minutes, ...
                                        squeeze(errors(:, end, :)), ...
                                        monte_carlo.covariances(:, :, end), ...
                                        updates(end).posterior_state);
            predicted_endpoint = ...
                MonteCarlo.makeEndpoint('prediction only', ...
                                        config.duration * to_minutes, final_errors, ...
                                        filter_run.final_covariance, ...
                                        filter_run.final_state);
            monte_carlo.endpoints = [update_endpoint, predicted_endpoint];
            [monte_carlo.bounds.lower, monte_carlo.bounds.upper] = ...
                coa.filter.meanChiSquareBounds(run_count, [blocks.dof], ...
                                               config.nees_confidence);
        end

        function endpoint = makeEndpoint(kind, minutes, errors, covariance, state)
            % Final-time snapshot: errors [runs x 18] (truth - estimate), the reference
            % run's predicted covariance, and its estimate at that time.
            endpoint = struct('label', sprintf('t = %.0f min (%s)', minutes, kind), ...
                              'errors', errors, 'covariance', covariance, ...
                              'reference_state', state);
        end

        function blocks = neesBlocks()
            % Whole-state NEES plus one block per body (position and velocity).
            blocks = struct('name', "All states", 'rows', 1:18, 'dof', 18);
            names = ["Counterweight (1)", "Catch mechanism (2)", "Payload (3)"];
            for body = 1:3
                offset = 3 * (body - 1);
                blocks(end + 1) = struct('name', names(body), ...
                                         'rows', [offset + (1:3), ...
                                                  9 + offset + (1:3)], ...
                                         'dof', 6);
            end
        end

    end
end
