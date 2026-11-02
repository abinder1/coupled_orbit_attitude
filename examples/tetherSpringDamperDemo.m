% Copyright (c) 2026 Andrew Binder

% tetherSpringDamperDemo
%   State-only comparison of the rigid (constraint-based) tether model,
%   coa.orbit.tetherAccel, with the spring-damper tether model,
%   coa.orbit.tetherSpringDamperDynamics, for the m1-m2 pair in the CR3BP.
%
%   The spring is initialized at equilibrium: its rest length is chosen so
%   that, at t = 0, the spring force equals the tension the rigid model needs
%   to hold the separation at exactly tether_length_km under the prescribed
%   spin. A stiffness (k) and damping-ratio (zeta) sweep then reports how
%   closely the spring model tracks the rigid one, and how each model handles
%   initial states that violate the length constraint.
%
%   Studies
%     1. Nominal sweep: tracking error against the rigid trajectory, length
%        deviation, tension, and solver cost versus (k, zeta).
%     2. Perturbed starts: the same random length-violating initial errors are
%        given to the rigid model and to every (k, zeta); counts integration
%        failures and the worst length excursion.
%     3. Time histories of length and tension for selected stiffnesses.

%% Run the comparison
clear;
clc;
close('all');
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'toolbox'));

config = demoConfig();
problem = makeProblem(config);

fprintf('Rigid reference and spring equilibrium\n');
rigid = integrateRigid(problem, problem.initial_state);
fprintf('  Rigid: %d steps, %.2f s\n', numel(rigid.solution.x), rigid.wall_time);
fprintf('  Equilibrium tension per unit total mass: %.4g\n', ...
        problem.equilibrium_force);

mode_count = numel(config.feedforward_modes);
nominal = cell(1, mode_count);
perturbed = cell(1, mode_count);
for mode = 1:mode_count
    feedforward = config.feedforward_modes(mode);
    label = modeLabel(feedforward);
    fprintf('\n=== %s ===\n', label);
    nominal{mode} = runNominalSweep(config, problem, rigid, feedforward);
    printNominalTable(config, problem, nominal{mode}, feedforward);

    fprintf('\nPerturbed starts (%d draws, %g km, %g m/s per component)\n', ...
            config.perturbed_draws, config.perturbed_position_km, ...
            config.perturbed_velocity_mps);
    perturbed{mode} = runPerturbedSweep(config, problem, feedforward);
    printPerturbedTable(config, perturbed{mode});
end

if config.show_plots
    for mode = 1:mode_count
        plotSweepMaps(config, nominal{mode}, perturbed{mode}, ...
                      modeLabel(config.feedforward_modes(mode)));
    end
    plotTimeHistories(config, problem, rigid);
end

%% Configuration
function config = demoConfig()
    % Physical, numerical, and sweep parameters for this example.
    config.length_scale = 384400;             % Earth-Moon length unit [km]
    config.time_scale = 375190.258892627;     % Earth-Moon time unit [s]
    config.mu = 0.012150585;                  % Moon's mass fraction

    config.tether_length_km = 200;
    config.spin_rate = 2 * pi / (55 * 60);    % Inertial rate [rad/s]
    config.mass_2_fraction = 0.5;
    config.tether_length_rate = 0;
    config.tether_direction = [1; 0; 0];
    config.halo_state = [0.831987145310294; 0; 0.12575552414079; ...
                         0; 0.239887839725451; 0];
    config.duration_s = 2 * 3600;
    config.duration = config.duration_s / config.time_scale;

    % Stiffness is force per unit total mass per unit length (nondimensional);
    % damping is expressed as a ratio of the critical damping of the pair.
    config.stiffness_values = [1e2, 1e3, 1e4, 1e5, 3e5, 1e6];
    config.damping_ratios = [0.05, 0.2, 0.5, 1, 2];
    config.history_stiffness_values = [1e3, 1e4, 1e5];
    config.history_damping_ratio = 0.5;

    config.perturbed_draws = 40;
    config.perturbed_position_km = 10;
    config.perturbed_velocity_mps = 1;
    config.seed = 1;                          % Same draws for every mode
    config.history_samples = 601;

    % false: rest length offset by the static stretch under the spin load.
    % true: centripetal/gravity-gradient feedforward, rest length = tether length.
    config.feedforward_modes = true;

    config.rigid_options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
    config.spring_options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
    config.show_plots = true;
end

function problem = makeProblem(config)
    % Initial state, dynamics handles, and equilibrium spring force.
    problem.mass_2_fraction = config.mass_2_fraction;
    problem.length = config.tether_length_km / config.length_scale;
    problem.length_rate = config.tether_length_rate;
    problem.duration = config.duration;
    problem.accel_fn = @(x) coa.orbit.cr3bp(x, config.mu);

    spin_rate = config.spin_rate * config.time_scale;
    initial_args = {config.halo_state, problem.length, ...
                    config.mass_2_fraction, spin_rate, ...
                    config.tether_length_rate, config.tether_direction};
    problem.initial_state = tetherInitialState(initial_args{:});

    problem.rigid_rhs = @(t, x) rigidDynamics(x, problem);
    problem.equilibrium_force = rigidTension(problem.initial_state, problem);
end

%% Dynamics and integration
function state_derivative = rigidDynamics(state, problem)
    % Rigid-constraint dynamics of [p1; p2; v1; v2] from coa.orbit.tetherAccel.
    body_1 = [state(1:3); state(7:9)];
    body_2 = [state(4:6); state(10:12)];
    args = {problem.mass_2_fraction, problem.length_rate, problem.accel_fn};
    acceleration = coa.orbit.tetherAccel(body_1, body_2, args{:});
    state_derivative = [state(7:12); acceleration];
end

function force = rigidTension(state, problem)
    % Spring-equivalent force per unit total mass that the constraint applies:
    % eta * (1 - eta) * (u . c), with c the correction in tetherAccel.
    eta = problem.mass_2_fraction;
    body_1 = [state(1:3); state(7:9)];
    body_2 = [state(4:6); state(10:12)];
    args = {eta, problem.length_rate, problem.accel_fn};
    constrained = coa.orbit.tetherAccel(body_1, body_2, args{:});
    correction = (constrained(1:3) - problem.accel_fn(body_1)) / eta;
    separation = state(4:6) - state(1:3);
    force = eta * (1 - eta) * dot(separation / norm(separation), correction);
end

function result = integrateRigid(problem, initial_state)
    % Integrate the rigid model; flag a stall instead of erroring.
    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
    result = integrateModel(@ode113, problem.rigid_rhs, initial_state, ...
                            problem.duration, options);
end

function result = integrateSpring(problem, initial_state, parameters, options)
    % Integrate the spring-damper model for one (stiffness, damping) pair.
    args = {problem.mass_2_fraction, parameters.rest_length, ...
            problem.length_rate, parameters.stiffness, parameters.damping, ...
            problem.accel_fn, parameters.feedforward};
    rhs = @(t, x) coa.orbit.tetherSpringDamperDynamics(t, x, args{:});
    result = integrateModel(@ode15s, rhs, initial_state, problem.duration, ...
                            options);
end

function result = integrateModel(solver, rhs, initial_state, duration, options)
    % Run one solver and report whether it reached the end of the scenario.
    warning_state = warning('off', 'all');
    cleanup = onCleanup(@() warning(warning_state));
    timer = tic;
    try
        result.solution = solver(rhs, [0, duration], initial_state, options);
        result.success = result.solution.x(end) >= duration * (1 - 1e-12) && ...
                         all(isfinite(result.solution.y(:)));
    catch error
        if strcmp(error.identifier, 'MATLAB:UndefinedFunction') || ...
                contains(error.identifier, 'validators') || ...
                contains(error.identifier, 'invalidInput')
            rethrow(error);
        end
        result.solution = struct('x', 0, 'y', initial_state);
        result.success = false;
    end
    result.wall_time = toc(timer);
end

function parameters = springParameters(problem, stiffness, damping_ratio, ...
                                       feedforward)
    % Rest length and damping coefficient for one sweep point. Without
    % feedforward, the rest length is offset so the spring carries the spin
    % tension at t = 0; with it, the rest length is the tether length.
    eta = problem.mass_2_fraction;
    parameters.stiffness = stiffness;
    parameters.damping = 2 * damping_ratio * sqrt(stiffness * eta * (1 - eta));
    parameters.feedforward = feedforward;
    parameters.static_stretch = ~feedforward * problem.equilibrium_force / ...
                                stiffness;
    parameters.rest_length = problem.length - parameters.static_stretch;
    parameters.frequency = sqrt(stiffness / (eta * (1 - eta)));
    parameters.valid = parameters.rest_length > 0;
end

function label = modeLabel(feedforward)
    % Human-readable name of a feedforward setting.
    if feedforward
        label = 'with centripetal feedforward';
    else
        label = 'spring-damper only (equilibrium rest length)';
    end
end

%% Study 1: nominal sweep
function nominal = runNominalSweep(config, problem, rigid, feedforward)
    % Compare every (k, zeta) with the rigid trajectory from the same start.
    time = linspace(0, problem.duration, config.history_samples);
    rigid_states = deval(rigid.solution, time);
    rigid_length = separationLength(rigid_states);

    shape = [numel(config.damping_ratios), numel(config.stiffness_values)];
    fields = {'max_length_dev_km', 'rms_length_dev_km', 'max_position_err_km', ...
              'final_position_err_km', 'steps', 'wall_time'};
    for name = fields
        nominal.(name{1}) = nan(shape);
    end

    for column = 1:numel(config.stiffness_values)
        for row = 1:numel(config.damping_ratios)
            stiffness = config.stiffness_values(column);
            damping_ratio = config.damping_ratios(row);
            parameters = springParameters(problem, stiffness, damping_ratio, ...
                                          feedforward);
            if ~parameters.valid
                continue
            end
            run = integrateSpring(problem, problem.initial_state, parameters, ...
                                  config.spring_options);
            if ~run.success
                continue
            end
            states = deval(run.solution, time);
            length_dev = (separationLength(states) - problem.length) * ...
                         config.length_scale;
            position_err = positionError(states, rigid_states, config);

            nominal.max_length_dev_km(row, column) = max(abs(length_dev));
            nominal.rms_length_dev_km(row, column) = rms(length_dev);
            nominal.max_position_err_km(row, column) = max(position_err);
            nominal.final_position_err_km(row, column) = position_err(end);
            nominal.steps(row, column) = numel(run.solution.x);
            nominal.wall_time(row, column) = run.wall_time;
        end
    end
    nominal.time = time;
    nominal.rigid_length_dev_km = (rigid_length - problem.length) * ...
                                  config.length_scale;
end

function length = separationLength(states)
    % Distance between bodies 1 and 2 for each column of states.
    length = vecnorm(states(4:6, :) - states(1:3, :), 2, 1);
end

function error_km = positionError(states, reference_states, config)
    % Largest body-position difference (bodies 1 and 2) at each time.
    error_1 = vecnorm(states(1:3, :) - reference_states(1:3, :), 2, 1);
    error_2 = vecnorm(states(4:6, :) - reference_states(4:6, :), 2, 1);
    error_km = max(error_1, error_2) * config.length_scale;
end

function printNominalTable(config, problem, nominal, feedforward)
    % One line per (k, zeta): static stretch, frequency, tracking, and cost.
    fprintf(['Nominal sweep. Columns: k, wn/spin, stretch [km], zeta, ', ...
             'max |L-L0| [km], max pos err [km], final pos err [km], steps\n']);
    spin = config.spin_rate * config.time_scale;
    for column = 1:numel(config.stiffness_values)
        parameters = springParameters(problem, ...
                                      config.stiffness_values(column), 1, feedforward);
        for row = 1:numel(config.damping_ratios)
            if ~parameters.valid
                fprintf('  k=%-8.0e rest length <= 0 (stretch %.0f km); skipped\n', ...
                        parameters.stiffness, ...
                        parameters.static_stretch * config.length_scale);
                break
            end
            fprintf(['  k=%-8.0e wn/spin=%-7.1f stretch=%-8.3g zeta=%-5g ', ...
                     'maxdev=%-9.3g maxpos=%-9.3g finalpos=%-9.3g steps=%d\n'], ...
                    parameters.stiffness, parameters.frequency / spin, ...
                    parameters.static_stretch * config.length_scale, ...
                    config.damping_ratios(row), ...
                    nominal.max_length_dev_km(row, column), ...
                    nominal.max_position_err_km(row, column), ...
                    nominal.final_position_err_km(row, column), ...
                    nominal.steps(row, column));
        end
    end
end

%% Study 2: perturbed starts
function perturbed = runPerturbedSweep(config, problem, feedforward)
    % Give the rigid model and every (k, zeta) identical length-violating starts.
    rng(config.seed, 'twister');
    position_sigma = config.perturbed_position_km / config.length_scale;
    velocity_sigma = config.perturbed_velocity_mps / 1000 * ...
                     config.time_scale / config.length_scale;
    errors = [position_sigma * randn(6, config.perturbed_draws); ...
              velocity_sigma * randn(6, config.perturbed_draws)];

    perturbed.rigid_failures = 0;
    for draw = 1:config.perturbed_draws
        run = integrateRigid(problem, problem.initial_state + errors(:, draw));
        perturbed.rigid_failures = perturbed.rigid_failures + ~run.success;
    end

    shape = [numel(config.damping_ratios), numel(config.stiffness_values)];
    perturbed.failures = nan(shape);
    perturbed.worst_length_dev_km = nan(shape);
    for column = 1:numel(config.stiffness_values)
        for row = 1:numel(config.damping_ratios)
            stiffness = config.stiffness_values(column);
            damping_ratio = config.damping_ratios(row);
            parameters = springParameters(problem, stiffness, damping_ratio, ...
                                          feedforward);
            if ~parameters.valid
                continue
            end
            failures = 0;
            worst = 0;
            for draw = 1:config.perturbed_draws
                start_state = problem.initial_state + errors(:, draw);
                run = integrateSpring(problem, start_state, parameters, ...
                                      config.spring_options);
                failures = failures + ~run.success;
                if run.success
                    distance = separationLength(run.solution.y);
                    worst = max(worst, max(abs(distance - problem.length)));
                end
            end
            perturbed.failures(row, column) = failures;
            perturbed.worst_length_dev_km(row, column) = ...
                worst * config.length_scale;
        end
    end
end

function printPerturbedTable(config, perturbed)
    % Failure counts for the rigid model and for each spring setting.
    variable_names = compose('k_%.0e', config.stiffness_values);
    row_names = compose('zeta_%g', config.damping_ratios);
    fprintf('  Rigid model: %d/%d integrations failed\n', ...
            perturbed.rigid_failures, config.perturbed_draws);
    fprintf('  Spring-damper failures (rows zeta, columns k):\n');
    disp(array2table(perturbed.failures, 'VariableNames', variable_names, ...
                     'RowNames', row_names));
    fprintf('  Worst |L - L0| [km] among successful spring runs:\n');
    disp(array2table(perturbed.worst_length_dev_km, ...
                     'VariableNames', variable_names, 'RowNames', row_names));
end

%% Plots
function plotSweepMaps(config, nominal, perturbed, mode_label)
    % Heat maps over (k, zeta) of tracking error, length, cost, and robustness.
    figure('Name', ['Spring-damper sweep: ' mode_label], 'NumberTitle', 'off');
    layout = tiledlayout(2, 3, 'TileSpacing', 'compact');
    title(layout, sprintf('Spring-damper vs rigid tether (%g km), %s', ...
                          config.tether_length_km, mode_label));

    sweepMap(config, nominal.max_length_dev_km, 'log10 max |L - L0| [km]', true);
    sweepMap(config, nominal.max_position_err_km, ...
             'log10 max body position error vs rigid [km]', true);
    sweepMap(config, nominal.final_position_err_km, ...
             'log10 final position error vs rigid [km]', true);
    sweepMap(config, nominal.steps, 'ode15s steps', false);
    sweepMap(config, perturbed.failures, ...
             sprintf('Perturbed failures out of %d (rigid: %d)', ...
                     config.perturbed_draws, perturbed.rigid_failures), false);
    sweepMap(config, perturbed.worst_length_dev_km, ...
             'log10 worst |L - L0| in perturbed runs [km]', true);
    drawnow;
end

function sweepMap(config, values, label, use_log)
    % One heat map with stiffness on x and damping ratio on y.
    nexttile;
    if use_log
        values = log10(values);
    end
    imagesc(values);
    set(gca, 'YDir', 'normal');
    xticks(1:numel(config.stiffness_values));
    xticklabels(compose('%.0e', config.stiffness_values));
    yticks(1:numel(config.damping_ratios));
    yticklabels(compose('%g', config.damping_ratios));
    xlabel('stiffness k');
    ylabel('damping ratio \zeta');
    title(label, 'FontSize', 9);
    colorbar;
end

function plotTimeHistories(config, problem, rigid)
    % Length and tension versus time for selected stiffnesses against rigid.
    time = linspace(0, problem.duration, config.history_samples);
    minutes = time * config.time_scale / 60;
    rigid_states = deval(rigid.solution, time);
    rigid_force = arrayfun(@(i) rigidTension(rigid_states(:, i), problem), ...
                           1:numel(time));

    figure('Name', 'Spring-damper time histories', 'NumberTitle', 'off');
    layout = tiledlayout(2, 1, 'TileSpacing', 'compact');
    length_axes = nexttile;
    hold(length_axes, 'on');
    force_axes = nexttile;
    hold(force_axes, 'on');
    rigid_length_dev = (separationLength(rigid_states) - problem.length) * ...
                       config.length_scale;
    plot(length_axes, minutes, rigid_length_dev, 'k', 'LineWidth', 1.5, ...
         'DisplayName', 'rigid');
    plot(force_axes, minutes, rigid_force, 'k', 'LineWidth', 1.5, ...
         'DisplayName', 'rigid');

    for stiffness = config.history_stiffness_values
        for feedforward = config.feedforward_modes
            damping_ratio = config.history_damping_ratio;
            parameters = springParameters(problem, stiffness, damping_ratio, ...
                                          feedforward);
            if ~parameters.valid
                continue
            end
            run = integrateSpring(problem, problem.initial_state, parameters, ...
                                  config.spring_options);
            if ~run.success
                continue
            end
            states = deval(run.solution, time);
            label = sprintf('k = %.0e', stiffness);
            length_dev = (separationLength(states) - problem.length) * ...
                         config.length_scale;
            plot(length_axes, minutes, length_dev, 'DisplayName', label);
            plot(force_axes, minutes, ...
                 springForce(states, time, parameters, problem), ...
                 'DisplayName', label);
        end
    end

    ylabel(length_axes, 'L - L_0 [km]');
    title(length_axes, sprintf('Length deviation, \\zeta = %g', ...
                               config.history_damping_ratio));
    legend(length_axes, 'Location', 'best');
    ylabel(force_axes, 'Tension per unit total mass');
    xlabel(force_axes, 'Time [min]');
    title(force_axes, 'Tether force');
    grid(length_axes, 'on');
    grid(force_axes, 'on');
    title(layout, 'Rigid vs spring-damper tether');
    drawnow;
end

function force = springForce(states, time, parameters, problem)
    % Total tether force per unit total mass at each column of states,
    % including the feedforward term when the parameters enable it.
    separation = states(4:6, :) - states(1:3, :);
    distance = vecnorm(separation, 2, 1);
    stretch_rate = sum(separation ./ distance .* (states(10:12, :) - ...
                                                  states(7:9, :)), 1);
    rest = parameters.rest_length + problem.length_rate * time;
    force = parameters.stiffness * (distance - rest) + ...
            parameters.damping * (stretch_rate - problem.length_rate);

    if parameters.feedforward
        eta = problem.mass_2_fraction;
        for i = 1:numel(time)
            a_1 = problem.accel_fn(states([1:3, 7:9], i));
            a_2 = problem.accel_fn(states([4:6, 10:12], i));
            relative_velocity = states(10:12, i) - states(7:9, i);
            kappa = dot(relative_velocity, relative_velocity) - ...
                    problem.length_rate^2;
            load = dot(separation(:, i), a_2 - a_1) + kappa;
            force(i) = force(i) + eta * (1 - eta) * load / distance(i);
        end
    end
end

%% Initial state
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
