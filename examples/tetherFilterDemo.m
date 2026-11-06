% Copyright (c) 2026 Andrew Binder

% tetherFilterDemo
%   Demonstrates a continuous-prediction, discrete-update EKF for a tethered
%   three-body system, with evenly spaced measurement updates using all
%   implemented measurement groups. The last update occurs a configurable
%   margin before the end of the scenario (the payload rendezvous).
%
%   State convention (nondimensional, CR3BP rotating frame):
%     x = [p1; p2; p3; v1; v2; v3], with each block 3x1.
%     Body 1 = counterweight, 2 = catch mechanism, 3 = payload.
%   The tether truth and EKF use the critically damped spring-damper model;
%   process noise is included in covariance propagation, but not simulated as
%   a random realization in the truth trajectory. Measurement updates use an
%   iterated EKF (config.ekf_iterations) because the long-baseline crosslink
%   measurements are strongly nonlinear over the initial position uncertainty.
%
%   Truth is the nominal initial condition; each Monte Carlo run starts its
%   filter from an estimate sampled from the initial covariance.

%% Run the EKF demonstration
clear;
clc;
close('all');
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'toolbox'));
rng(1, 'twister');

config = demoConfig();
preparePlotOutput(config);
[initial_state, initial_covariance, initial_factor] = makeInitialEstimate(config);
scenario = propagateScenario(initial_state, config);
plotTruthTrajectory(scenario, config);
% Truth is the nominal initial condition; each trial's estimate is a sample
% drawn around it with the initial covariance.
filter_run = runFilter(scenario, ...
                       perturbInitialState(initial_state, initial_factor), ...
                       initial_covariance, config);

reportFilterRun(filter_run, config);
plotPredictionHistory(filter_run.history, config);

fprintf(['\nReference filter run complete. Starting %d Monte Carlo ' ...
         'runs...\n'], config.monte_carlo_runs);
monte_carlo = runMonteCarlo(scenario, initial_state, initial_covariance, ...
                            initial_factor, filter_run, config);
fprintf('Monte Carlo runs complete. Generating reports and visualizations...\n');
runMonteCarloAnalyses(monteCarloAnalyses(), monte_carlo, config);

%% Scenario configuration and initialization
function config = demoConfig()
    % Physical, numerical, and plotting parameters for this example.
    config.length_scale = 384400;            % Earth-Moon length unit [km]
    config.time_scale = 375190.258892627;     % Earth-Moon time unit [s]
    config.mu = 0.012150585;                  % Moon's mass fraction

    config.tether_length_km = 200;
    config.spin_rate = 2 * pi / (55 * 60);    % Inertial rate [rad/s]
    config.mass_2_fraction = 0.5;
    config.tether_length_rate = 0;
    config.tether_direction = [1; 0; 0];
    natural_frequency = config.spin_rate * config.time_scale;
    initial_length = config.tether_length_km / config.length_scale;
    config.spring_damper = struct('initial_length', initial_length, ...
                                  'stiffness', natural_frequency^2, ...
                                  'damping', 2 * natural_frequency);
    config.halo_state = [0.831987145310294; 0; 0.12575552414079; ...
                         0; 0.239887839725451; 0];

    config.duration_s = 24 * 3600;
    config.duration = config.duration_s / config.time_scale;
    config.measurement_count = 6;
    % The last update occurs this long before the end, when the payload
    % reaches the catch mechanism and relative measurements are singular-prone.
    config.final_measurement_margin_min = 30;
    last_measurement_s = config.duration_s - ...
                         60 * config.final_measurement_margin_min;
    update_count = config.measurement_count;
    config.measurement_times = (1:update_count) * ...
                               (last_measurement_s / update_count) / ...
                               config.time_scale;
    config.initial_sigma = struct('com_position', 100, ...       % [km]
                                  'com_velocity', 1, ...          % [m/s]
                                  'relative_position', 10 / 1000, ... % [km]
                                  'relative_velocity', 0.1, ...     % [m/s]
                                  'payload_position', 100, ...    % [km]
                                  'payload_velocity', 1);         % [m/s]
    config.process_noise = struct('tension_psd', 1e-6, ...
                                  'transverse_psd', 1e-8, ...
                                  'body_3_acceleration_psd', 1e-8);

    config.measurement_mask = 2^8 - 1;
    config.ekf_iterations = 5;                % Max IEKF relinearizations per update
    config.residual_sigma_floor = 0.1;        % Residual sigma floor [x sensor sigma]
    config.monte_carlo_runs = 100;            % Full-EKF runs for NEES/ANEES/NIS
    config.nees_confidence = 0.95;
    config.corner_ellipse_sigma = 3;          % Mahalanobis radius of corner ellipses
    config.plot_samples_per_interval = 31;

    % Output mode per plot (keys listed in plotKeys):
    %   "off"  - skip the plot, including its data preparation
    %   "save" - render hidden, write an image to plot_output_folder, close it
    %   "show" - display on screen and leave open
    config.plot_modes = struct('prediction_history', "off", ...
                               'trajectory', "show", ...
                               'anees', "off", ...
                               'rmse', "off", ...
                               'innovations', "off", ...
                               'residuals', "off", ...
                               'corner', "off", ...
                               'relative_corner', "off");
    config.plot_output_folder = fullfile(fileparts(mfilename('fullpath')), ...
                                         '..', 'runs', 'tetherFilterDemo');
    config.plot_file_format = "png";
    config.plot_resolution = 150;             % dpi for saved images

    config.position_rows = 1:9;
    config.velocity_rows = 10:18;
    config.body_names = ["Counterweight (1)", ...
                         "Catch mechanism (2)", "Payload (3)"];
end

function [state, covariance, factor] = makeInitialEstimate(config)
    % Initialize the tethered bodies, rendezvous payload, and covariance.
    accel_fn = @(x) coa.orbit.cr3bp(x, config.mu);

    tether_length = config.tether_length_km / config.length_scale;
    spin_rate = config.spin_rate * config.time_scale;
    tether_rate = config.tether_length_rate;
    tether_state = tetherInitialState(config.halo_state, tether_length, ...
                                      config.mass_2_fraction, spin_rate, ...
                                      tether_rate, config.tether_direction);
    state = rendezvousInitialState(tether_state, config, accel_fn);

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
    duration = config.duration;
    tether_rhs = dynamicsFunction(config);

    % Body 3 is a placeholder here; it does not affect the tether dynamics.
    placeholder_state = [tether_state(1:6); tether_state(4:6); ...
                         tether_state(7:12); tether_state(10:12)];
    tether_solution = ode113(tether_rhs, [0, duration], placeholder_state, ...
                             odeset('RelTol', 1e-10, 'AbsTol', 1e-12));
    final_tether_state = tether_solution.y(:, end);

    catch_state = [final_tether_state(4:6); final_tether_state(13:15)];
    payload_state = propagateSingleBody(catch_state, -duration, accel_fn);

    % Reorder [p1; p2; v1; v2; p3; v3] as [p1; p2; p3; v1; v2; v3].
    combined_state = [tether_state; payload_state];
    state = [combined_state(1:6); combined_state(13:15); ...
             combined_state(7:12); combined_state(16:18)];
end

%% Truth trajectory
function scenario = propagateScenario(initial_state, config)
    % Propagate the nominal initial condition; every trial shares this truth.
    scenario.truth_solution = ode113(dynamicsFunction(config), ...
                                     [0, config.duration], initial_state, ...
                                     odeset('RelTol', 1e-10, 'AbsTol', 1e-12));
end

function rhs = dynamicsFunction(config, process_noise)
    % Tethered three-body dynamics; optional process noise gives the packed
    % state-and-covariance dynamics used by the filter. State-only and
    % covariance propagation both use the same spring-damper tether model.
    accel_fn = @(x) coa.orbit.cr3bp(x, config.mu);
    if nargin < 2
        process_noise = struct();
    end
    args = {config.mass_2_fraction, config.tether_length_rate, accel_fn, ...
            process_noise, config.spring_damper};
    rhs = @(t, x) coa.orbit.tetheredThreeBodyEkfPropagation(t, x, args{:});
end

function state = perturbInitialState(initial_state, initial_factor)
    % Draw a zero-mean initial error with the covariance used by the filter.
    state = initial_state + initial_factor * randn(18, 1);
end

%% Filter propagation and measurement updates
function filter_run = runFilter(scenario, initial_state, initial_covariance, config)
    % Propagate to each measurement, update, then continue from the posterior.
    filter_rhs = dynamicsFunction(config, config.process_noise);
    abs_tolerance = [1e-10 * ones(18, 1); 1e-24 * ones(324, 1)];
    filter_options = odeset('RelTol', 1e-9, 'AbsTol', abs_tolerance);

    length_scale = config.length_scale;
    velocity_scale = length_scale / config.time_scale;
    covariance_fn = @coa.meas.measurementCovariances;
    measurement_covariances = covariance_fn(length_scale, velocity_scale);
    packed_state = [initial_state; initial_covariance(:)];
    previous_time = 0;
    filter_run.history = emptyHistory();

    for k = 1:config.measurement_count
        measurement_time = config.measurement_times(k);
        filter_solution = ode113(filter_rhs, [previous_time, measurement_time], ...
                                 packed_state, filter_options);
        sample_interval = linspace(previous_time, measurement_time, ...
                                   config.plot_samples_per_interval);
        filter_samples = deval(filter_solution, sample_interval);
        truth_samples = deval(scenario.truth_solution, sample_interval);
        first_sample = 1 + (k > 1);
        samples = first_sample:numel(sample_interval);
        filter_run.history = appendPredictionHistory(filter_run.history, ...
                                                     sample_interval(samples), ...
                                                     filter_samples(:, samples), ...
                                                     truth_samples(:, samples), config);

        prior_state = filter_solution.y(1:18, end);
        prior_covariance = unpackCovariance(filter_solution.y(19:end, end));
        truth_state = deval(scenario.truth_solution, measurement_time);
        update = measurementUpdate(prior_state, prior_covariance, ...
                                   truth_state, measurement_covariances, config);
        filter_run.updates(k, 1) = update;

        packed_state = [update.posterior_state; ...
                        update.posterior_covariance(:)];
        filter_run.history = appendPredictionHistory(filter_run.history, ...
                                                     measurement_time, packed_state, ...
                                                     truth_state, config);
        previous_time = measurement_time;
    end

    % Predict (no update) from the last measurement to the end of the scenario.
    final_packed = packed_state;
    if previous_time < config.duration
        filter_solution = ode113(filter_rhs, [previous_time, config.duration], ...
                                 packed_state, filter_options);
        sample_interval = linspace(previous_time, config.duration, ...
                                   config.plot_samples_per_interval);
        filter_samples = deval(filter_solution, sample_interval);
        truth_samples = deval(scenario.truth_solution, sample_interval);
        samples = 2:numel(sample_interval);
        filter_run.history = appendPredictionHistory(filter_run.history, ...
                                                     sample_interval(samples), ...
                                                     filter_samples(:, samples), ...
                                                     truth_samples(:, samples), config);
        final_packed = filter_solution.y(:, end);
    end
    filter_run.final_state = final_packed(1:18);
    filter_run.final_covariance = unpackCovariance(final_packed(19:end));
end

function update = measurementUpdate(prior_state, prior_covariance, truth_state, ...
                                    measurement_covariances, config)
    % Simulate one noisy observation from truth, then apply the shared update.
    mask = config.measurement_mask;
    model = @(state) coa.meas.measurementModels(state, mask, ...
                                                measurement_covariances);
    [true_measurement, ~, noise_covariance] = model(truth_state);
    noise_chol = chol(noise_covariance, 'lower');
    observation = true_measurement + noise_chol * randn(size(true_measurement));

    [posterior_state, posterior_covariance, info] = ...
        coa.filter.ekfUpdate(prior_state, prior_covariance, observation, ...
                             model, angleRows(mask), config.ekf_iterations);

    update = info;
    update.truth_state = truth_state;
    update.prior_state = prior_state;
    update.posterior_state = posterior_state;
    update.posterior_covariance = posterior_covariance;
    update.noise_chol = noise_chol;
end

function rows = angleRows(mask)
    % Rows of the stacked measurement vector holding right ascension and
    % declination (groups 1 and 2), which need wrapping to [-pi, pi).
    rows = (1:(6 * bitand(mask, 1) + 4 * bitand(floor(mask / 2), 1)))';
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
    history.position_error = [history.position_error, ...
                              state_error(config.position_rows, :) * position_scale];
    history.velocity_error = [history.velocity_error, ...
                              state_error(config.velocity_rows, :) * velocity_scale];
    position_sigma = sqrt(max(variance(config.position_rows, :), 0));
    velocity_sigma = sqrt(max(variance(config.velocity_rows, :), 0));
    history.position_sigma = [history.position_sigma, position_sigma * position_scale];
    history.velocity_sigma = [history.velocity_sigma, velocity_sigma * velocity_scale];
end

%% Monte Carlo consistency check (NEES / ANEES)
function monte_carlo = runMonteCarlo(scenario, initial_state, initial_covariance, ...
                                     initial_factor, filter_run, config)
    % Repeat the full EKF against the fixed truth, with a fresh initial estimate
    % sample and measurement noise per run. Every run integrates its
    % own covariance and computes its own gains, so ANEES and NIS test the
    % filter as implemented.
    filter_rhs = dynamicsFunction(config, config.process_noise);
    abs_tolerance = [1e-10 * ones(18, 1); 1e-24 * ones(324, 1)];
    filter_options = odeset('RelTol', 1e-9, 'AbsTol', abs_tolerance);
    velocity_scale = config.length_scale / config.time_scale;
    measurement_covariances = ...
        coa.meas.measurementCovariances(config.length_scale, velocity_scale);

    blocks = neesBlocks();
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
        sampled_state = perturbInitialState(initial_state, initial_factor);
        packed_state = [sampled_state; initial_covariance(:)];
        previous_time = 0;

        for k = 1:step_count
            time = config.measurement_times(k);
            solution = ode113(filter_rhs, [previous_time, time], packed_state, ...
                              filter_options);
            truth_state = deval(truth_solution, time);
            update = measurementUpdate(solution.y(1:18, end), ...
                                       unpackCovariance(solution.y(19:end, end)), ...
                                       truth_state, measurement_covariances, config);
            packed_state = [update.posterior_state; update.posterior_covariance(:)];
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
    monte_carlo.errors = errors;                  % [run, update, state]: truth - est.
    monte_carlo.innovations = innovations;        % [run, update, measurement row]
    monte_carlo.residuals = residuals;            % post-fit, same shape
    monte_carlo.reference_updates = updates;      % the single filter run shown in plots
    monte_carlo.covariances = cat(3, updates.posterior_covariance);
    monte_carlo.blocks = blocks;
    monte_carlo.nees = computeNees(errors, run_covariances, blocks);
    monte_carlo.anees = squeeze(mean(monte_carlo.nees, 1));
    monte_carlo.nis = nis;                        % [run, update]
    monte_carlo.anis = mean(nis, 1)';
    monte_carlo.nis_dof = measurement_count;
    monte_carlo.time_minutes = config.measurement_times * config.time_scale / 60;
    monte_carlo.run_count = run_count;
    % Final-time snapshots: the last measurement update, then the prediction-only
    % end of the scenario. Corner-type analyses loop over these.
    to_minutes = config.time_scale / 60;
    update_endpoint = makeEndpoint('final update', ...
                                   config.measurement_times(end) * to_minutes, ...
                                   squeeze(errors(:, end, :)), ...
                                   monte_carlo.covariances(:, :, end), ...
                                   updates(end).posterior_state);
    predicted_endpoint = makeEndpoint('prediction only', ...
                                      config.duration * to_minutes, ...
                                      final_errors, filter_run.final_covariance, ...
                                      filter_run.final_state);
    monte_carlo.endpoints = [update_endpoint, predicted_endpoint];
    monte_carlo.bounds = neesBounds(run_count, [blocks.dof], ...
                                    config.nees_confidence);
end

function endpoint = makeEndpoint(kind, minutes, errors, covariance, state)
    % Final-time snapshot: errors [runs x 18] (truth - estimate), the reference
    % run's predicted covariance, and its estimate at that time.
    endpoint = struct('label', sprintf('t = %.0f min (%s)', minutes, kind), ...
                      'errors', errors, 'covariance', covariance, ...
                      'reference_state', state);
end

function nees = computeNees(errors, covariances, blocks)
    % NEES per run, update, and block from stored errors and covariances.
    % covariances is [18, 18, update, run] (each run's own posterior).
    [run_count, step_count, ~] = size(errors);
    nees = zeros(run_count, step_count, numel(blocks));
    for run = 1:run_count
        for k = 1:step_count
            for b = 1:numel(blocks)
                rows = blocks(b).rows;
                e = squeeze(errors(run, k, rows));
                nees(run, k, b) = e' * (covariances(rows, rows, k, run) \ e);
            end
        end
    end
end

function analyses = monteCarloAnalyses()
    % Registry of Monte Carlo analyses, run in order. To add one, write a
    % function taking (monte_carlo, config) -- a report, a plot, or both -- and
    % append its handle(s) here.
    analyses = {@reportAnees, @reportNis, @plotAnees, @plotCorner, ...
                @plotInnovations, @plotResiduals, @plotRmse, ...
                @reportRelativeState, @plotRelativeCorner};
end

function runMonteCarloAnalyses(analyses, monte_carlo, config)
    % Run every registered analysis on the shared Monte Carlo results.
    for analysis = analyses
        analysis{1}(monte_carlo, config);
    end
end

function blocks = neesBlocks()
    % Whole-state NEES plus one block per body (position and velocity).
    blocks = struct('name', "All states", 'rows', 1:18, 'dof', 18);
    names = ["Counterweight (1)", "Catch mechanism (2)", "Payload (3)"];
    for body = 1:3
        offset = 3 * (body - 1);
        blocks(end + 1) = struct('name', names(body), ...
                                 'rows', [offset + (1:3), 9 + offset + (1:3)], ...
                                 'dof', 6);
    end
end

function bounds = neesBounds(run_count, dof, confidence)
    % Two-sided chi-square acceptance region for ANEES, via gammaincinv.
    tail = (1 - confidence) / 2;
    total_dof = run_count * dof;
    bounds.lower = 2 * gammaincinv(tail, total_dof / 2) / run_count;
    bounds.upper = 2 * gammaincinv(1 - tail, total_dof / 2) / run_count;
end

function reportAnees(monte_carlo, config)
    % Summarize ANEES against its acceptance region for each block.
    fprintf('\nANEES over %d full-EKF Monte Carlo runs (%.0f%% region)\n', ...
            monte_carlo.run_count, 100 * config.nees_confidence);
    fprintf('%-20s %5s  %9s  %-19s  %s\n', 'Block', 'dof', 'Mean ANEES', ...
            'Acceptance region', 'Steps inside');

    for b = 1:numel(monte_carlo.blocks)
        anees = monte_carlo.anees(:, b);
        lower = monte_carlo.bounds.lower(b);
        upper = monte_carlo.bounds.upper(b);
        inside = mean(anees >= lower & anees <= upper);
        fprintf('%-20s %5d  %9.3g  [%7.3g, %7.3g]  %5.1f%%\n', ...
                monte_carlo.blocks(b).name, monte_carlo.blocks(b).dof, ...
                mean(anees), lower, upper, 100 * inside);
    end
end

function reportNis(monte_carlo, config)
    % Average NIS per update against its chi-square acceptance region.
    bounds = neesBounds(monte_carlo.run_count, monte_carlo.nis_dof, ...
                        config.nees_confidence);
    fprintf('\nANIS per update (dof %d, region [%.1f, %.1f])\n', ...
            monte_carlo.nis_dof, bounds.lower, bounds.upper);
    for k = 1:numel(monte_carlo.anis)
        fprintf('  t = %7.1f min: ANIS = %9.3g\n', monte_carlo.time_minutes(k), ...
                monte_carlo.anis(k));
    end
end

function plotAnees(monte_carlo, config)
    % Plot per-body ANEES versus time with the chi-square acceptance region.
    if plotMode(config, 'anees') == "off"
        return
    end
    fig = newFigure(config, 'anees', 'Name', 'ANEES consistency', 'Color', 'w');
    body_blocks = 2:numel(monte_carlo.blocks);
    layout = tiledlayout(1, numel(body_blocks), ...
                         'TileSpacing', 'compact', 'Padding', 'compact');
    time = monte_carlo.time_minutes;
    colors = lines(numel(body_blocks));

    for panel = 1:numel(body_blocks)
        b = body_blocks(panel);
        ax = nexttile(layout);
        hold(ax, 'on');
        % Normalize by the block dimension so the expected value is 1.
        dof = monte_carlo.blocks(b).dof;
        lower = monte_carlo.bounds.lower(b) / dof;
        upper = monte_carlo.bounds.upper(b) / dof;
        fill(ax, [time, fliplr(time)], ...
             [lower * ones(size(time)), upper * ones(size(time))], ...
             [0.5, 0.5, 0.5], 'FaceAlpha', 0.2, 'EdgeColor', [0.45, 0.45, 0.45], ...
             'LineWidth', 0.5);
        yline(ax, 1, 'k--');
        plot(ax, time, monte_carlo.anees(:, b) / dof, '-o', ...
             'Color', colors(panel, :), 'LineWidth', 1.1, 'MarkerSize', 3);
        set(ax, 'YScale', 'log');
        grid(ax, 'on');
        title(ax, sprintf('%s (n = %d)', monte_carlo.blocks(b).name, ...
                          monte_carlo.blocks(b).dof));
        ylabel(ax, 'ANEES / n');
    end

    xlabel(layout, 'Time [min]');
    title(layout, sprintf(['Normalized ANEES over %d runs with %.0f%% acceptance ' ...
                           'region (dashed: expected value 1)'],  ...
                          monte_carlo.run_count, 100 * config.nees_confidence));
    finishFigure(fig, config, 'anees');
end

%% RMSE versus predicted standard deviation
function plotRmse(monte_carlo, config)
    % Dimensional RMSE over the Monte Carlo runs vs. the filter's predicted
    % sigma at each update. Columns are bodies; rows are position (km) and
    % velocity (m/s) axes. RMSE above the predicted sigma means overconfident.
    if plotMode(config, 'rmse') == "off"
        return
    end
    rmse = squeeze(sqrt(mean(monte_carlo.errors.^2, 1)));   % [update, state]
    sigma = sqrt(max(reshape(diagBlocks(monte_carlo.covariances), 18, []), 0))';
    position_scale = config.length_scale;                      % km
    velocity_scale = config.length_scale / config.time_scale * 1000;  % m/s
    scale = [position_scale * ones(1, 9), velocity_scale * ones(1, 9)];
    rmse = rmse .* scale;
    sigma = sigma .* scale;

    time = monte_carlo.time_minutes;
    axis_names = ["x", "y", "z"];
    fig = newFigure(config, 'rmse', 'Name', 'RMSE vs predicted sigma', ...
                    'Color', 'w', 'Position', [80, 40, 1000, 900]);
    layout = tiledlayout(6, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for row = 1:6
        for body = 1:3
            is_velocity = row > 3;
            axis_index = row - 3 * is_velocity;
            state = 9 * is_velocity + 3 * (body - 1) + axis_index;
            ax = nexttile(layout);
            hold(ax, 'on');
            plot(ax, time, sigma(:, state), 'k--', 'LineWidth', 1.1);
            plot(ax, time, rmse(:, state), '-o', 'LineWidth', 1.1, ...
                 'MarkerSize', 3, 'Color', [0, 0.447, 0.741]);
            set(ax, 'FontSize', 7);
            grid(ax, 'on');
            if body == 1
                prefix = ["p", "v"];
                unit = ["km", "m/s"];
                ylabel(ax, sprintf('%s%s [%s]', prefix(1 + is_velocity), ...
                                   axis_names(axis_index), unit(1 + is_velocity)));
            end
            if row == 1
                title(ax, config.body_names(body), 'FontSize', 9);
            end
            if row < 6
                set(ax, 'XTickLabel', []);
            end
        end
    end
    lg = legend(nexttile(1), {'predicted \sigma', 'RMSE'});
    lg.Layout.Tile = 'east';
    lg.FontSize = 8;
    xlabel(layout, 'Time [min]');
    title(layout, sprintf(['RMSE over %d runs vs. predicted \\sigma ' ...
                           '(RMSE above \\sigma: overconfident)'], ...
                          monte_carlo.run_count));
    finishFigure(fig, config, 'rmse');
end

function d = diagBlocks(covariances)
    % Diagonals of every 18x18 page, returned as [18, 1, pages].
    pages = size(covariances, 3);
    d = zeros(18, 1, pages);
    for k = 1:pages
        d(:, 1, k) = diag(covariances(:, :, k));
    end
end

%% Measurement innovation and residual histories
function plotInnovations(monte_carlo, config)
    % Pre-fit innovations, normalized by sqrt(diag(S)).
    plotMeasurementHistory(monte_carlo, config, 'innovation', ...
                           'innovation_covariance', "Pre-fit innovations");
end

function plotResiduals(monte_carlo, config)
    % Post-fit residuals, normalized by sqrt(diag(R - H P+ H')), floored at
    % config.residual_sigma_floor times the sensor sigma.
    plotMeasurementHistory(monte_carlo, config, 'residual', ...
                           'residual_covariance', "Post-fit residuals");
end

function plotMeasurementHistory(monte_carlo, config, value_field, covariance_field, ...
                                title_text)
    % One tile per active sensor group. Values are normalized by their predicted
    % sigma, so every group shares the same axis and the +/-3 bounds mean the
    % same thing everywhere. Gray dots are all Monte Carlo runs (reused from
    % the consistency check); colored lines are the single reference filter run.
    % The plot key is the plural of value_field ('innovations' or 'residuals').
    key = plural(value_field);
    if plotMode(config, key) == "off"
        return
    end
    updates = monte_carlo.reference_updates;
    sigma = sqrt(max(cell2mat(arrayfun(@(u) diag(u.(covariance_field)), ...
                                       updates', 'UniformOutput', false)), 0));
    % Post-fit variances R - H P+ H' can vanish (or round negative) for
    % near-exact measurements; floor sigma at a fraction of the sensor sigma so
    % the normalized values stay finite. Innovations are never affected, since
    % S >= R.
    noise_sigma = sqrt(cell2mat(arrayfun(@(u) sum(u.noise_chol.^2, 2), ...
                                         updates', 'UniformOutput', false)));
    sigma = max(sigma, config.residual_sigma_floor * noise_sigma);
    reference = cell2mat(arrayfun(@(u) u.(value_field), updates', ...
                                  'UniformOutput', false)) ./ sigma;
    samples = monte_carlo.(plural(value_field));
    time = monte_carlo.time_minutes;
    groups = activeMeasurementGroups(config.measurement_mask);

    fig = newFigure(config, key, 'Name', title_text, 'Color', 'w', ...
                    'Position', [60, 60, 1100, 800]);
    layout = tiledlayout(ceil(numel(groups) / 2), 2, ...
                         'TileSpacing', 'compact', 'Padding', 'compact');
    for g = 1:numel(groups)
        ax = nexttile(layout);
        hold(ax, 'on');
        rows = groups(g).rows;
        colors = lines(numel(rows));
        for r = 1:numel(rows)
            normalized = squeeze(samples(:, :, rows(r))) ./ sigma(rows(r), :);
            plot(ax, repmat(time, size(normalized, 1), 1), normalized, '.', ...
                 'Color', [0.7, 0.7, 0.7], 'MarkerSize', 4);
        end
        yline(ax, [-3, 3], 'k--');
        for r = 1:numel(rows)
            plot(ax, time, reference(rows(r), :), 'o', 'Color', colors(r, :), ...
                 'LineWidth', 1, 'MarkerSize', 3);
        end
        grid(ax, 'on');
        title(ax, sprintf('%s (%d rows)', groups(g).name, numel(rows)), ...
              'FontSize', 9);
        ylim(ax, [-5, 5]);
    end
    xlabel(layout, 'Time [min]');
    ylabel(layout, 'Value / predicted \sigma');
    title(layout, sprintf(['%s: reference run (colored, one per row) over %d ' ...
                           'Monte Carlo runs (gray); dashed: \\pm3\\sigma'], ...
                          title_text, monte_carlo.run_count));
    finishFigure(fig, config, key);
end

function name = plural(field)
    % Name of the Monte Carlo array holding a per-update field.
    name = [field, 's'];
end

function groups = activeMeasurementGroups(mask)
    % Selected sensor groups with their rows in the stacked measurement vector,
    % in the order used by coa.meas.measurementModels.
    names = ["Ground optical angles", "Camera / star-tracker angles", ...
             "Camera angular size", "DSN range", "DSN range rate", ...
             "Crosslink range", "Crosslink range rate", "GNSS position"];
    sizes = [6, 4, 1, 3, 3, 2, 2, 9];
    groups = struct('name', {}, 'rows', {});
    offset = 0;
    for group = 1:numel(names)
        if bitand(mask, 2^(group - 1)) > 0
            groups(end + 1) = struct('name', names(group), ...
                                     'rows', offset + (1:sizes(group))); %#ok<AGROW>
            offset = offset + sizes(group);
        end
    end
end

%% Results and visualizations
function plotCorner(monte_carlo, config)
    % 18-D corner plot of the full-state errors at each final-time endpoint.
    if plotMode(config, 'corner') == "off"
        return
    end
    colors = lines(3);
    colors = colors(arrayfun(@bodyOf, 1:18), :);
    for endpoint = monte_carlo.endpoints
        drawCornerPlot(config, 'corner', ...
                       "Error corner plot, " + endpoint.label, endpoint.errors, ...
                       endpoint.covariance, stateLabels(), colors, ...
                       sprintf('Error corner plot at %s, %d runs', endpoint.label, ...
                               monte_carlo.run_count));
    end
end

function drawCornerPlot(config, key, name, errors, covariance, labels, colors, ...
                        title_text)
    % Generic corner plot of [runs x n] errors against an n x n covariance,
    % standardized by predicted sigma. Standardizing gives every axis the same
    % +/-limit range, so cells can be packed edge to edge with tick labels only
    % on the outer border. Labels and colors are per state. The figure is shown
    % or saved according to config.plot_modes.(key).
    %   Diagonal: marginal histogram vs. the predicted N(0,1).
    %   Lower:    samples with the predicted covariance ellipse.
    %   Upper:    sample minus predicted correlation (blue/red = under/over).
    radius = config.corner_ellipse_sigma;
    sigma = sqrt(diag(covariance));
    z = errors ./ sigma';
    predicted_corr = covariance ./ (sigma * sigma');
    sample_corr = corrcoef(z);
    n = size(z, 2);
    font_size = max(6, min(10, 60 / n));
    marker_size = max(3, min(10, 54 / n));

    limit = 4;
    corr_limit = 0.5;
    cmap = [linspace(0, 1, 32)', linspace(0, 1, 32)', ones(32, 1); ...
            ones(32, 1), linspace(1, 0, 32)', linspace(1, 0, 32)'];
    angle = linspace(0, 2 * pi, 100);
    unit_circle = [cos(angle); sin(angle)];
    edges = linspace(-limit, limit, 17);
    grid_points = linspace(-limit, limit, 200);
    bar_x = [edges(1:end - 1); edges(2:end); edges(2:end); edges(1:end - 1)];

    screen = get(0, 'ScreenSize');
    side = 0.9 * min(screen(3:4));
    fig = newFigure(config, key, 'Name', name, 'Color', 'w', ...
                    'Position', [screen(1) + 20, screen(2) + 40, side, side]);
    extra = 0.1 * (size(labels, 1) > 1);
    margin = struct('left', 0.06 + extra, 'right', 0.10, ...
                    'bottom', 0.06 + 0.5 * extra, 'top', 0.09);
    cell_w = (1 - margin.left - margin.right) / n;
    cell_h = (1 - margin.bottom - margin.top) / n;

    for i = 1:n
        for j = 1:n
            position = [margin.left + (j - 1) * cell_w, ...
                        margin.bottom + (n - i) * cell_h, cell_w, cell_h];
            ax = axes(fig, 'Position', position, 'Box', 'on', ...
                      'FontSize', font_size, 'TickLength', [0.03, 0.03], ...
                      'XLim', [-limit, limit], 'YLim', [-limit, limit], ...
                      'XTick', [-2, 0, 2], 'YTick', [-2, 0, 2], ...
                      'NextPlot', 'add'); %#ok<LAXES>

            if i == j
                height = histcounts(z(:, j), edges, 'Normalization', 'pdf');
                patch(ax, bar_x, [0; 0; 1; 1] * height, colors(j, :), ...
                      'EdgeColor', 'none', 'FaceAlpha', 0.8);
                plot(ax, grid_points, exp(-grid_points.^2 / 2) / sqrt(2 * pi), ...
                     'k', 'LineWidth', 0.8);
                set(ax, 'YLim', [0, max(0.6, 1.05 * max(height))], 'YTick', []);
            elseif i > j
                scatter(ax, z(:, j), z(:, i), marker_size, colors(i, :), ...
                        'filled', 'MarkerFaceAlpha', 0.5);
                pair = [j, i];
                ellipse = radius * chol(predicted_corr(pair, pair), 'lower') * ...
                          unit_circle;
                plot(ax, ellipse(1, :), ellipse(2, :), 'k-', 'LineWidth', 0.7);
            else
                difference = sample_corr(i, j) - predicted_corr(i, j);
                level = (min(max(difference, -corr_limit), corr_limit) + ...
                         corr_limit) / (2 * corr_limit);
                set(ax, 'Color', cmap(1 + round(level * (size(cmap, 1) - 1)), :), ...
                    'XTick', [], 'YTick', []);
                text(ax, 0, 0, sprintf('%+.2f', difference), ...
                     'FontSize', font_size, 'HorizontalAlignment', 'center');
            end

            if i < n
                set(ax, 'XTickLabel', []);
            end
            if j > 1 || i == 1
                set(ax, 'YTickLabel', []);
            end
            if i == n
                xlabel(ax, labels(:, j), 'FontSize', font_size + 1, ...
                       'Color', colors(j, :));
            end
            if j == 1
                ylabel(ax, labels(:, i), 'FontSize', font_size + 1, 'Rotation', 0, ...
                       'HorizontalAlignment', 'right', 'Color', colors(i, :));
            end
        end
    end

    colormap(fig, cmap);
    colorbar_ax = axes(fig, 'Position', [0, 0, 1e-3, 1e-3], 'Visible', 'off');
    caxis(colorbar_ax, [-corr_limit, corr_limit]);
    bar = colorbar(colorbar_ax, 'Position', ...
                   [1 - margin.right + 0.03, margin.bottom + 0.3, 0.012, 0.3]);
    bar.Label.String = 'Sample - predicted correlation (upper triangle)';
    bar.FontSize = 7;

    sgtitle(fig, {sprintf('%s, axes in predicted \\sigma', title_text), ...
                  sprintf(['lower: samples with predicted %g\\sigma ellipse | ' ...
                           'diagonal: histogram vs. N(0,1) | ' ...
                           'upper: correlation mismatch'], radius)}, ...
            'FontSize', 10);
    finishFigure(fig, config, key);
end

%% Relative payload-to-catch-mechanism state in the tether frame
function relative = relativeStateDistribution(endpoint, config)
    % Errors and covariance of d = [r3 - r2; v3 - v2], expressed in the tether
    % frame of the reference run's estimate at the endpoint, in meters and mm/s.
    % d is an affine (here linear) function of the state, so the final-time errors
    % and covariance map exactly: e_d = M e, P_d = M P M'. No new samples needed.
    frame = tetherFrame(endpoint.reference_state);

    difference = zeros(6, 18);
    difference(1:3, 4:6) = -eye(3);
    difference(1:3, 7:9) = eye(3);
    difference(4:6, 13:15) = -eye(3);
    difference(4:6, 16:18) = eye(3);

    units = [config.length_scale * 1e3 * ones(3, 1); ...
             config.length_scale / config.time_scale * 1e6 * ones(3, 1)];
    map = diag(units) * blkdiag(frame, frame) * difference;
    covariance = map * endpoint.covariance * map';

    relative.frame = frame;
    relative.errors = endpoint.errors * map';
    relative.covariance = (covariance + covariance') / 2;
    relative.names = ["\Delta r_{e1}", "\Delta r_{e2}", "\Delta r_{e3}", ...
                      "\Delta v_{e1}", "\Delta v_{e2}", "\Delta v_{e3}"];
    relative.units = ["m", "m", "m", "mm/s", "mm/s", "mm/s"];
end

function frame = tetherFrame(state)
    % Rows are the frame axes in the state's coordinates: e1 along r2 - r1,
    % e3 along (r2 - r1) x (v2 - v1), and e2 = e3 x e1 completing the triad.
    tether = state(4:6) - state(1:3);
    relative_velocity = state(13:15) - state(10:12);
    normal = cross(tether, relative_velocity);
    e1 = tether / norm(tether);
    e3 = normal / norm(normal);
    frame = [e1, cross(e3, e1), e3]';
end

function reportRelativeState(monte_carlo, config)
    % Predicted vs. sample spread of the relative state in the tether frame.
    for endpoint = monte_carlo.endpoints
        relative = relativeStateDistribution(endpoint, config);
        predicted = sqrt(diag(relative.covariance));
        sample = sqrt(mean(relative.errors.^2, 1))';
        fprintf('\nRelative state r3-r2, v3-v2 in tether frame at %s\n', ...
                endpoint.label);
        fprintf('%-8s %-6s  %12s  %12s  %s\n', 'Comp.', 'Unit', 'Predicted sd', ...
                'Sample RMS', 'Ratio');
        for k = 1:6
            fprintf('%-8s %-6s  %12.4g  %12.4g  %5.2f\n', ...
                    char(erase(relative.names(k), {'\', '{', '}'})), ...
                    relative.units(k), predicted(k), sample(k), ...
                    sample(k) / predicted(k));
        end
    end
end

function plotRelativeCorner(monte_carlo, config)
    % 6-D corner plot of the relative state in the tether frame, per endpoint.
    if plotMode(config, 'relative_corner') == "off"
        return
    end
    colors = [repmat([0, 0.447, 0.741], 3, 1); repmat([0.85, 0.325, 0.098], 3, 1)];
    for endpoint = monte_carlo.endpoints
        relative = relativeStateDistribution(endpoint, config);
        sigma = sqrt(diag(relative.covariance));
        labels = strings(2, 6);
        for k = 1:6
            labels(1, k) = relative.names(k);
            labels(2, k) = sprintf('\\sigma=%.2g %s', sigma(k), relative.units(k));
        end
        name = "Relative state corner plot (tether frame), " + endpoint.label;
        drawCornerPlot(config, 'relative_corner', name, ...
                       relative.errors, relative.covariance, labels, colors, ...
                       sprintf(['Error of [r_3-r_2; v_3-v_2] in the tether frame ' ...
                                'at %s, %d runs'], endpoint.label, ...
                               monte_carlo.run_count));
    end
end

function body = bodyOf(index)
    % Body (1-3) owning a state index; positions then velocities, 3 per body.
    body = mod(floor((index - 1) / 3), 3) + 1;
end

function labels = stateLabels()
    % Short labels p<body><axis> and v<body><axis> in state-vector order.
    axes_names = ["x", "y", "z"];
    labels = strings(1, 18);
    for index = 1:18
        kind = "p";
        if index > 9
            kind = "v";
        end
        labels(index) = kind + bodyOf(index) + axes_names(mod(index - 1, 3) + 1);
    end
end

function reportFilterRun(filter_run, config)
    % Report prior/posterior errors and uncertainty for every measurement.
    updates = filter_run.updates;
    fprintf(['EKF demonstration: %.0f-second propagation with %d evenly ' ...
             'spaced measurement updates\n'], config.duration_s, ...
            config.measurement_count);
    fprintf('Measurements: all implemented coa.meas measurement groups\n');
    fprintf('Update  Time [min]  Innovation norm  Pos. error [km]  Vel. error [m/s]\n');

    for k = 1:numel(updates)
        prior_error = updates(k).prior_state - updates(k).truth_state;
        posterior_error = updates(k).posterior_state - updates(k).truth_state;
        pos_rows = config.position_rows;
        vel_rows = config.velocity_rows;
        fprintf('%6d  %10.2f  %15.3g  %14.3g -> %-8.3g  %14.3g -> %-8.3g\n', ...
                k, config.measurement_times(k) * config.time_scale / 60, ...
                norm(updates(k).innovation), ...
                stateNorm(prior_error, pos_rows, config, 'position'), ...
                stateNorm(posterior_error, pos_rows, config, 'position'), ...
                stateNorm(prior_error, vel_rows, config, 'velocity'), ...
                stateNorm(posterior_error, vel_rows, config, 'velocity'));
    end
end

function plotPredictionHistory(history, config)
    % Plot signed XYZ errors and shaded componentwise +/-3-sigma bounds.
    if plotMode(config, 'prediction_history') == "off"
        return
    end
    time_minutes = history.time * config.time_scale / 60;
    component_colors = lines(3);
    component_names = ["x", "y", "z"];
    position_limit = symmetricPlotLimit( ...
                                        history.position_error, history.position_sigma);
    velocity_limit = symmetricPlotLimit( ...
                                        history.velocity_error, history.velocity_sigma);

    for body = 1:3
        body_rows = 3 * (body - 1) + (1:3);
        fig = newFigure(config, 'prediction_history', 'Name', ...
                        config.body_names(body) + " signed errors", 'Color', 'w');
        layout = tiledlayout(3, 2, 'TileSpacing', 'compact', ...
                             'Padding', 'compact');

        for component = 1:3
            row = body_rows(component);
            position_axes = nexttile(layout, 2 * component - 1);
            plotSignedComponent(position_axes, time_minutes, ...
                                history.position_error(row, :), ...
                                history.position_sigma(row, :), config, ...
                                component_colors(component, :));
            ylim(position_axes, position_limit);
            title(position_axes, component_names(component) + " position");
            ylabel(position_axes, 'Position error [km]');

            velocity_axes = nexttile(layout, 2 * component);
            plotSignedComponent(velocity_axes, time_minutes, ...
                                history.velocity_error(row, :), ...
                                history.velocity_sigma(row, :), config, ...
                                component_colors(component, :));
            ylim(velocity_axes, velocity_limit);
            title(velocity_axes, component_names(component) + " velocity");
            ylabel(velocity_axes, 'Velocity error [m/s]');
        end

        xlabel(layout, 'Time [min]');
        title(layout, config.body_names(body) + ...
              " signed component errors with shaded \pm3\sigma bounds");
        finishFigure(fig, config, 'prediction_history');
    end
end

function limits = symmetricPlotLimit(errors, sigmas)
    % Include all signed errors and three-sigma envelopes in a shared range.
    extent = max([abs(errors(:)); 3 * sigmas(:)]);
    extent = max(extent, eps(extent));
    limits = 1.05 * [-extent, extent];
end

function plotSignedComponent(ax, time, error, sigma, config, color)
    % Draw one signed component error and its +/-3-sigma envelope.
    hold(ax, 'on');

    bound = 3 * sigma;
    fill(ax, [time, fliplr(time)], [bound, -fliplr(bound)], ...
         color, 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off');
    boundary_color = [0.45, 0.45, 0.45];
    plot(ax, time, bound, '-', 'Color', boundary_color, 'LineWidth', 0.5, ...
         'HandleVisibility', 'off');
    plot(ax, time, -bound, '-', 'Color', boundary_color, 'LineWidth', 0.5, ...
         'HandleVisibility', 'off');

    update_times = config.measurement_times * config.time_scale / 60;
    xline(ax, update_times, ':', 'Color', [0.65, 0.65, 0.65], ...
          'HandleVisibility', 'off');

    plot(ax, time, error, 'Color', color, 'LineWidth', 1.1);

    yline(ax, 0, 'k-', 'HandleVisibility', 'off');
    grid(ax, 'on');
end

function plotTruthTrajectory(scenario, config)
    % Plot the nominal truth paths, payload-to-catch range, and inbound speed.
    if plotMode(config, 'trajectory') == "off"
        return
    end

    sample_count = config.measurement_count * ...
                   config.plot_samples_per_interval + 1;
    sample_time = linspace(0, config.duration, sample_count);
    truth = deval(scenario.truth_solution, sample_time);
    position = truth(1:9, :) * config.length_scale;
    velocity = truth(10:18, :) * config.length_scale / config.time_scale;
    time_minutes = sample_time * config.time_scale / 60;

    relative_position = position(7:9, :) - position(4:6, :);
    relative_velocity = velocity(7:9, :) - velocity(4:6, :);
    payload_to_catch_range = vecnorm(relative_position);
    inbound_speed = -sum(relative_position .* relative_velocity, 1) ./ ...
                    payload_to_catch_range * 1000;

    fig = newFigure(config, 'trajectory', 'Name', ...
                    'Truth motion and payload range', 'Color', 'w', ...
                    'Position', [60, 60, 1200, 720]);
    layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', ...
                         'Padding', 'compact');

    trajectory_axes = nexttile(layout, 1, [2, 1]);
    hold(trajectory_axes, 'on');
    colors = lines(3);
    bodies = {1:3, 4:6, 7:9};
    names = ["Counterweight (m1)", "Catch mechanism (m2)", ...
             "Payload (m3)"];
    for body = 1:numel(bodies)
        rows = bodies{body};
        plot3(trajectory_axes, position(rows(1), :), position(rows(2), :), ...
              position(rows(3), :), 'Color', colors(body, :), ...
              'LineWidth', 1.2);
    end
    grid(trajectory_axes, 'on');
    axis(trajectory_axes, 'equal');
    view(trajectory_axes, 3);
    xlabel(trajectory_axes, 'X [km]');
    ylabel(trajectory_axes, 'Y [km]');
    zlabel(trajectory_axes, 'Z [km]');
    title(trajectory_axes, 'Truth trajectories in the rotating frame');
    legend(trajectory_axes, names, 'Location', 'best');

    range_axes = nexttile(layout, 2);
    plot(range_axes, time_minutes, payload_to_catch_range, 'LineWidth', 1.2);
    grid(range_axes, 'on');
    ylabel(range_axes, 'Range |r_3 - r_2| [km]');
    title(range_axes, 'Payload (m3) to catch mechanism (m2)');

    speed_axes = nexttile(layout, 4);
    plot(speed_axes, time_minutes, inbound_speed, 'LineWidth', 1.2);
    yline(speed_axes, 0, 'k-', 'HandleVisibility', 'off');
    grid(speed_axes, 'on');
    xlabel(speed_axes, 'Time [min]');
    ylabel(speed_axes, 'Inbound speed [m/s]');
    title(speed_axes, 'Positive when closing');

    sgtitle(layout, 'Truth motion and payload-to-catch range history');
    finishFigure(fig, config, 'trajectory');
end

%% Plot output control
function keys = plotKeys()
    % Every plot that can be switched between "off", "save", and "show".
    keys = ["prediction_history", "trajectory", "anees", "rmse", "innovations", ...
            "residuals", "corner", "relative_corner"];
end

function preparePlotOutput(config)
    % Fail before the long run on a mistyped plot setting; create the output
    % folder if any plot is saved.
    given = string(fieldnames(config.plot_modes))';
    unknown = setdiff(given, plotKeys());
    missing = setdiff(plotKeys(), given);
    if ~isempty(unknown) || ~isempty(missing)
        error('tetherFilterDemo:PlotModeKeys', ...
              'config.plot_modes has unknown keys [%s] and lacks keys [%s].', ...
              strjoin(unknown, ', '), strjoin(missing, ', '));
    end

    modes = string(struct2cell(config.plot_modes));
    invalid = ~ismember(modes, ["off", "save", "show"]);
    if any(invalid)
        error('tetherFilterDemo:PlotModeValue', ...
              'Plot modes must be "off", "save", or "show"; invalid: %s.', ...
              strjoin(given(invalid'), ', '));
    end
    if any(modes == "save") && ~isfolder(config.plot_output_folder)
        mkdir(config.plot_output_folder);
    end
end

function mode = plotMode(config, key)
    % Output mode ("off", "save", or "show") configured for one plot.
    mode = string(config.plot_modes.(key));
end

function fig = newFigure(config, key, varargin)
    % Create a figure, hidden when it will only be saved. Rendering offscreen
    % avoids the cost of drawing a large figure to the display.
    visibility = 'on';
    if plotMode(config, key) == "save"
        visibility = 'off';
    end
    fig = figure(varargin{:}, 'Visible', visibility);
end

function finishFigure(fig, config, key)
    % Write a saved figure to file and close it; leave shown figures open.
    if plotMode(config, key) ~= "save"
        return
    end
    stem = strip(regexprep(string(fig.Name), '[^A-Za-z0-9]+', '_'), '_');
    file = fullfile(config.plot_output_folder, ...
                    key + "_" + stem + "." + config.plot_file_format);
    exportgraphics(fig, file, 'Resolution', config.plot_resolution);
    close(fig);
    fprintf('  Saved %s\n', file);
end

%% Shared unit and covariance helpers
function value = stateNorm(state_error, rows, config, quantity)
    % Convert a state-error norm to km or m/s.
    value = norm(state_error(rows));
    if strcmp(quantity, 'position')
        value = value * config.length_scale;
    else
        value = value * config.length_scale / config.time_scale * 1000;
    end
end

function covariance = unpackCovariance(packed_covariance)
    % Reshape a column-stacked covariance and remove integration asymmetry.
    covariance = reshape(packed_covariance, [18, 18]);
    covariance = (covariance + covariance') / 2;
end

%% Dynamics and initial-condition helpers
function state = propagateSingleBody(state_0, duration, accel_fn)
    % Coast a single body [p; v]; negative duration propagates backward.
    rhs = @(~, x) [x(4:6); accel_fn(x)];
    solution = ode113(rhs, [0, duration], state_0, ...
                      odeset('AbsTol', 1e-10, 'RelTol', 1e-10));
    state = solution.y(:, end);
end

function state = tetherInitialState(com_state, length_nd, mass_2_fraction, ...
                                    spin_rate_nd, length_dot, direction)
    % Initialize [p1; p2; v1; v2] around the tether center of mass.
    frame_rate = [0; 0; 1];
    relative_spin = spin_rate_nd * [0; 0; 1] - frame_rate;

    offset_1 = direction * (mass_2_fraction * length_nd);
    offset_2 = -direction * ((1 - mass_2_fraction) * length_nd);
    velocity_1 = cross(relative_spin, offset_1) + ...
                 mass_2_fraction * length_dot * direction;
    velocity_2 = cross(relative_spin, offset_2) - ...
                 (1 - mass_2_fraction) * length_dot * direction;

    state = [com_state(1:3) + offset_1; com_state(1:3) + offset_2; ...
             com_state(4:6) + velocity_1; com_state(4:6) + velocity_2];
end
