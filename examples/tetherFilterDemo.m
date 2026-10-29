% Copyright (c) 2026 Andrew Binder

% tetherFilterDemo
%   Demonstrates a continuous-prediction, discrete-update EKF for a tethered
%   three-body system over 24 hours, with 72 evenly spaced measurement
%   updates (one every 20 minutes) using all implemented measurement groups.
%
%   State convention (nondimensional, CR3BP rotating frame):
%     x = [p1; p2; p3; v1; v2; v3], with each block 3x1.
%     Body 1 = counterweight, 2 = catch mechanism, 3 = payload.
%
%   Truth begins with a random perturbation from the estimated initial state.
%   Process noise is included in covariance propagation, but not simulated as
%   a random realization in the truth trajectory.

%% Run the EKF demonstration
clear;
clc;
close('all');
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'toolbox'));
rng(1, 'twister');

config = demoConfig();
[initial_state, initial_covariance] = makeInitialEstimate(config);
scenario = propagateScenario(initial_state, initial_covariance, config);
filter_run = runFilter(scenario, initial_state, initial_covariance, config);

reportFilterRun(filter_run, config);
plotPredictionHistory(filter_run.history, filter_run.updates, config);

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
    config.halo_state = [0.831987145310294; 0; 0.12575552414079; ...
                         0; 0.239887839725451; 0];

    config.duration_s = 24 * 60 * 60;
    config.duration = config.duration_s / config.time_scale;
    config.measurement_count = 24 * 3;
    config.measurement_times = linspace( ...
                                        config.duration / config.measurement_count, ...
                                        config.duration, config.measurement_count);
    config.position_sigma_km = 0.01;
    config.velocity_sigma_mps = 1e-4;
    config.process_noise = struct( ...
                                  'tension_psd', 1e-8, ...
                                  'transverse_psd', 1e-8, ...
                                  'body_3_acceleration_psd', 1e-8);

    config.measurement_mask = 2^8 - 1;
    config.plot_samples_per_interval = 31;
    config.position_rows = 1:9;
    config.velocity_rows = 10:18;
    config.body_names = ["Counterweight (1)", ...
                         "Catch mechanism (2)", "Payload (3)"];
end

function [state, covariance] = makeInitialEstimate(config)
    % Initialize the tethered bodies, rendezvous payload, and covariance.
    accel_fn = @(x) coa.orbit.cr3bp(x, config.mu);
    propagate = @coa.orbit.tetheredThreeBodyEkfPropagation;

    tether_length = config.tether_length_km / config.length_scale;
    spin_rate = config.spin_rate * config.time_scale;
    tether_rate = config.tether_length_rate;
    tether_state = tetherInitialState(config.halo_state, tether_length, ...
                                      config.mass_2_fraction, spin_rate, ...
                                      tether_rate, config.tether_direction);
    state = rendezvousInitialState(tether_state, config, accel_fn, propagate);

    position_sigma = config.position_sigma_km / config.length_scale;
    velocity_sigma = config.velocity_sigma_mps / 1000 * ...
                     config.time_scale / config.length_scale;
    covariance = diag([repmat(position_sigma^2, 9, 1); ...
                       repmat(velocity_sigma^2, 9, 1)]);
end

function state = rendezvousInitialState(tether_state, config, accel_fn, propagate)
    % Back-propagate the catch mechanism to initialize a rendezvous payload.
    duration = config.duration;
    tether_rhs = @(t, x) propagate(t, x, config.mass_2_fraction, ...
                                   config.tether_length_rate, accel_fn);

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
function scenario = propagateScenario(initial_state, initial_covariance, config)
    % Propagate one perturbed truth trajectory through the measurement times.
    accel_fn = @(x) coa.orbit.cr3bp(x, config.mu);
    propagate = @coa.orbit.tetheredThreeBodyEkfPropagation;
    truth_state = perturbInitialState(initial_state, config);
    time_span = [0, config.duration];

    state_rhs = @(t, x) propagate(t, x, config.mass_2_fraction, ...
                                  config.tether_length_rate, accel_fn);
    scenario.truth_solution = ode113( ...
                                     state_rhs, time_span, truth_state, ...
                                     odeset('RelTol', 1e-10, 'AbsTol', 1e-12));
end

function state = perturbInitialState(initial_state, config)
    % Draw truth's initial error with the same componentwise sigmas as P0.
    position_sigma = config.position_sigma_km / config.length_scale;
    velocity_sigma = config.velocity_sigma_mps / 1000 * ...
                     config.time_scale / config.length_scale;
    initial_error = [position_sigma * randn(9, 1); ...
                     velocity_sigma * randn(9, 1)];
    state = initial_state + initial_error;
end

%% Filter propagation and measurement updates
function filter_run = runFilter(scenario, initial_state, initial_covariance, config)
    % Propagate to each measurement, update, then continue from the posterior.
    accel_fn = @(x) coa.orbit.cr3bp(x, config.mu);
    propagate = @coa.orbit.tetheredThreeBodyEkfPropagation;
    filter_rhs = @(t, packed) propagate(t, packed, config.mass_2_fraction, ...
                                        config.tether_length_rate, accel_fn, ...
                                        config.process_noise);
    abs_tolerance = [1e-10 * ones(18, 1); 1e-24 * ones(324, 1)];
    filter_options = odeset('RelTol', 1e-9, 'AbsTol', abs_tolerance);

    length_scale = config.length_scale;
    velocity_scale = length_scale / config.time_scale;
    covariance_fn = @coa.meas.measurementCovariances;
    measurement_covariances = covariance_fn(length_scale, velocity_scale);
    packed_state = [initial_state; initial_covariance(:)];
    previous_time = 0;
    filter_run.updates = repmat(emptyUpdate(), config.measurement_count, 1);
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
        update.time = measurement_time;
        update.integration_steps = numel(filter_solution.x);
        filter_run.updates(k) = update;

        packed_state = [update.posterior_state; ...
                        update.posterior_covariance(:)];
        filter_run.history = appendPredictionHistory(filter_run.history, ...
                                                     measurement_time, packed_state, ...
                                                     truth_state, config);
        previous_time = measurement_time;
    end
end

function update = measurementUpdate(prior_state, prior_covariance, truth_state, ...
                                    measurement_covariances, config)
    % Generate one noisy position observation and apply a Joseph-form update.
    model = @coa.meas.measurementModels;
    [predicted_measurement, H, R] = model(prior_state, config.measurement_mask, ...
                                          measurement_covariances);
    true_measurement = model(truth_state, config.measurement_mask);
    observation = true_measurement + ...
        chol(R, 'lower') * randn(size(true_measurement));

    innovation = observation - predicted_measurement;
    innovation_covariance = H * prior_covariance * H' + R;
    kalman_gain = (innovation_covariance \ (H * prior_covariance))';
    posterior_state = prior_state + kalman_gain * innovation;

    correction = eye(18) - kalman_gain * H;
    posterior_covariance = correction * prior_covariance * correction' + ...
                           kalman_gain * R * kalman_gain';

    update.truth_state = truth_state;
    update.prior_state = prior_state;
    update.posterior_state = posterior_state;
    update.prior_covariance = prior_covariance;
    update.posterior_covariance = ...
        (posterior_covariance + posterior_covariance') / 2;
    update.observation = observation;
    update.innovation = innovation;
end

function update = emptyUpdate()
    % Preallocate a measurement-update record with stable field shapes.
    update = struct('time', 0, ...
                    'truth_state', zeros(18, 1), ...
                    'prior_state', zeros(18, 1), ...
                    'posterior_state', zeros(18, 1), ...
                    'prior_covariance', zeros(18), ...
                    'posterior_covariance', zeros(18), ...
                    'observation', zeros(9, 1), ...
                    'innovation', zeros(9, 1), ...
                    'integration_steps', 0);
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
    for k = 1:numel(time)
        state = packed(1:18, k);
        covariance = unpackCovariance(packed(19:end, k));
        state_error = state - truth(:, k);

        history.time(end + 1) = time(k);
        history.position_error(:, end + 1) = ...
            state_error(config.position_rows) * config.length_scale;
        history.velocity_error(:, end + 1) = ...
            state_error(config.velocity_rows) * ...
            config.length_scale / config.time_scale * 1000;

        position_variance = diag(covariance(config.position_rows, ...
                                            config.position_rows));
        velocity_variance = diag(covariance(config.velocity_rows, ...
                                            config.velocity_rows));
        history.position_sigma(:, end + 1) = ...
            sqrt(max(position_variance, 0)) * config.length_scale;
        history.velocity_sigma(:, end + 1) = ...
            sqrt(max(velocity_variance, 0)) * ...
            config.length_scale / config.time_scale * 1000;
    end
end

%% Results and visualizations
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
                k, updates(k).time * config.time_scale / 60, ...
                norm(updates(k).innovation), ...
                stateNorm(prior_error, pos_rows, config, 'position'), ...
                stateNorm(posterior_error, pos_rows, config, 'position'), ...
                stateNorm(prior_error, vel_rows, config, 'velocity'), ...
                stateNorm(posterior_error, vel_rows, config, 'velocity'));
    end
end

function plotPredictionHistory(history, updates, config)
    % Plot signed XYZ errors and shaded componentwise +/-3-sigma bounds.
    time_minutes = history.time * config.time_scale / 60;
    component_colors = lines(3);
    component_names = ["x", "y", "z"];
    position_limit = symmetricPlotLimit( ...
                                        history.position_error, history.position_sigma);
    velocity_limit = symmetricPlotLimit( ...
                                        history.velocity_error, history.velocity_sigma);

    for body = 1:3
        body_rows = 3 * (body - 1) + (1:3);
        figure('Name', config.body_names(body) + " signed errors", ...
               'Color', 'w');
        layout = tiledlayout(3, 2, 'TileSpacing', 'compact', ...
                             'Padding', 'compact');

        for component = 1:3
            row = body_rows(component);
            position_axes = nexttile(layout, 2 * component - 1);
            plotSignedComponent(position_axes, time_minutes, ...
                                history.position_error(row, :), ...
                                history.position_sigma(row, :), updates, config, ...
                                component_colors(component, :));
            ylim(position_axes, position_limit);
            title(position_axes, component_names(component) + " position");
            ylabel(position_axes, 'Position error [km]');

            velocity_axes = nexttile(layout, 2 * component);
            plotSignedComponent(velocity_axes, time_minutes, ...
                                history.velocity_error(row, :), ...
                                history.velocity_sigma(row, :), updates, config, ...
                                component_colors(component, :));
            ylim(velocity_axes, velocity_limit);
            title(velocity_axes, component_names(component) + " velocity");
            ylabel(velocity_axes, 'Velocity error [m/s]');
        end

        xlabel(layout, 'Time [min]');
        title(layout, config.body_names(body) + ...
              " signed component errors with shaded \pm3\sigma bounds");
    end
end

function limits = symmetricPlotLimit(errors, sigmas)
    % Include all signed errors and three-sigma envelopes in a shared range.
    extent = max([abs(errors(:)); 3 * sigmas(:)]);
    extent = max(extent, eps(extent));
    limits = 1.05 * [-extent, extent];
end

function plotSignedComponent(ax, time, error, sigma, updates, config, color)
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

    for k = 1:numel(updates)
        update_time = updates(k).time * config.time_scale / 60;
        xline(ax, update_time, ':', 'Color', [0.65, 0.65, 0.65], ...
              'HandleVisibility', 'off');
    end

    plot(ax, time, error, 'Color', color, 'LineWidth', 1.1);

    yline(ax, 0, 'k-', 'HandleVisibility', 'off');
    grid(ax, 'on');
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
