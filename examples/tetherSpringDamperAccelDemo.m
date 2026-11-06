% Copyright (c) 2026 Andrew Binder

% tetherSpringDamperAccelDemo
%   States-only comparison of a critically damped spring-damper tether
%   (coa.orbit.tetherSpringDamperAccel, zeta = 1, c = 2 * sqrt(k)) against the
%   same tether with no spring or damper (k = c = 0).
%
%   Both tethers start with the same broken constraint: extra length and
%   length rate relative to the prescribed history L(t) = L_0 + Ldot * t. They
%   are integrated for one day in the CR3BP for every Ldot in
%   [-0.5, 0.5] m/s. With k = c = 0 only the feedforward terms act (free
%   accelerations and centripetal), so the length error never decays:
%   e(t) = e_0 + e_dot_0 * t. With the critically damped gains the error obeys
%   e_ddot + 2 * omega_n * e_dot + omega_n^2 * e = 0, whose exact solution is
%   e(t) = (e_0 + (e_dot_0 + omega_n * e_0) * t) * exp(-omega_n * t).
%
%   Figures compare the two tethers through
%     - length and length-rate error against the exact theory,
%     - tension and the spring, damper, and feedforward contributions,
%     - how far apart the two trajectories, spin rates, and orientations drift,
%     - the separation vector in the plane, and
%     - solver cost (ode113 steps, function evaluations, wall time).

%% Run the comparison
clear;
clc;
close('all');
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'toolbox'));

config = demoConfig();
problem = makeProblem(config);
runs = runCases(config, problem);
printComparison(config, runs);
if config.show_plots
    plotErrorHistories(config, problem, runs);
    plotForces(config, runs);
    plotTrajectoryDifference(config, runs);
    plotGeometry(config, runs);
    plotCost(config, runs);
end

%% Configuration
function config = demoConfig()
    % Physical, numerical, and comparison parameters for this example.
    config.length_scale = 384400;             % Earth-Moon length unit [km]
    config.time_scale = 375190.258892627;     % Earth-Moon time unit [s]
    config.mu = 0.012150585;                  % Moon's mass fraction

    config.tether_length_km = 200;
    config.spin_rate = 2 * pi / (55 * 60);    % Inertial rate [rad/s]
    config.mass_2_fraction = 0.5;
    config.tether_direction = [1; 0; 0];
    config.halo_state = [0.831987145310294; 0; 0.12575552414079; ...
                         0; 0.239887839725451; 0];
    config.duration_s = 24 * 3600;
    config.duration = config.duration_s / config.time_scale;

    % Broken constraint: the tether starts with this much extra length and
    % length rate relative to the reference L_0 and Ldot.
    config.initial_length_error_km = 5;
    config.initial_length_rate_error_mps = 0.1;

    % Critical damping: k = omega_n^2 and c = 2 * omega_n.
    config.frequency_ratio = 1;               % omega_n / spin rate
    config.length_rate_values_mps = [-0.5, -0.25, 0, 0.25, 0.5];
    config.settle_tolerance_m = 1;

    config.early_hours = 3;                   % Window for transient plots
    config.sample_count = 1441;
    config.ode_options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);

    % Case drawn in the single-Ldot panels.
    config.plot_length_rate_mps = 0.5;
    config.show_plots = true;
end

function problem = makeProblem(config)
    % Nondimensional constants shared by every integration.
    problem.eta = config.mass_2_fraction;
    problem.length = config.tether_length_km / config.length_scale;
    problem.accel_fn = @(x) coa.orbit.cr3bp(x, config.mu);
    problem.spin_rate = config.spin_rate * config.time_scale;

    problem.rate_scale = 1e-3 * config.time_scale / config.length_scale;
    problem.accel_scale = 1e3 * config.length_scale / config.time_scale^2;
    problem.length_error = config.initial_length_error_km / config.length_scale;
    problem.rate_error = config.initial_length_rate_error_mps * problem.rate_scale;
    problem.omega_n = config.frequency_ratio * problem.spin_rate;
end

%% Integration
function runs = runCases(config, problem)
    % Integrate the uncorrected (k = c = 0) and critically damped tethers.
    rates = config.length_rate_values_mps;
    stiffness = problem.omega_n^2;
    damping = 2 * problem.omega_n;
    runs.uncorrected = cell(numel(rates), 1);
    runs.critical = cell(numel(rates), 1);

    fprintf('Integrating %d cases over %g hours...\n', 2 * numel(rates), ...
            config.duration_s / 3600);
    for r = 1:numel(rates)
        runs.uncorrected{r} = integrateCase(config, problem, rates(r), 0, 0);
        runs.critical{r} = integrateCase(config, problem, rates(r), ...
                                         stiffness, damping);
    end
end

function run = integrateCase(config, problem, rate_mps, stiffness, damping)
    % Integrate one gain pair and sample the error, force, and geometry history.
    rate = rate_mps * problem.rate_scale;
    case_args = struct('rate', rate, 'stiffness', stiffness, 'damping', damping);
    rhs = @(t, x) springDamperDynamics(t, x, problem, case_args);
    state_0 = initialState(config, problem, rate);

    warning_state = warning('off', 'all');
    restore_warnings = onCleanup(@() warning(warning_state));
    timer = tic;
    solution = ode113(rhs, [0, config.duration], state_0, config.ode_options);
    run.wall_time = toc(timer);
    clear restore_warnings;

    run.completed = solution.x(end) >= config.duration;
    run.steps = numel(solution.x) - 1;
    run.function_evaluations = solution.stats.nfevals;
    run.stiffness = stiffness;
    run.damping = damping;

    end_time = min(solution.x(end), config.duration);
    % Fine early samples resolve the correction transient, whose time scale is
    % 1/omega_n and can be well under the coarse sample spacing.
    early_time = min(end_time, config.early_hours * 3600 / config.time_scale);
    times = unique([linspace(0, end_time, config.sample_count), ...
                    linspace(0, early_time, 721)]);
    states = deval(solution, times);
    run = sampleHistory(run, config, problem, times, states, case_args);
end

function run = sampleHistory(run, config, problem, times, states, case_args)
    % Length error, rate error, tension, and geometry at the sample times.
    sample_count = numel(times);
    length_error = zeros(1, sample_count);
    rate_error = zeros(1, sample_count);
    length_value = zeros(1, sample_count);
    tension = zeros(1, sample_count);
    spin = zeros(1, sample_count);
    for k = 1:sample_count
        state = states(:, k);
        separation = state(4:6) - state(1:3);
        relative_velocity = state(10:12) - state(7:9);
        length_value(k) = norm(separation);
        line_of_sight = separation / length_value(k);
        reference = problem.length + case_args.rate * times(k);
        length_error(k) = length_value(k) - reference;
        rate_error(k) = dot(line_of_sight, relative_velocity) - case_args.rate;
        [~, tension(k)] = springDamperDynamics(times(k), state, problem, case_args);

        % Inertial spin about z: rotating-frame rate plus the frame rate of 1.
        angular_momentum = cross(separation, relative_velocity);
        spin(k) = angular_momentum(3) / sum(separation(1:2).^2) + 1;
    end

    run.hours = times * config.time_scale / 3600;
    run.length_km = length_value * config.length_scale;
    run.separation_km = (states(4:6, :) - states(1:3, :)) * config.length_scale;
    run.position_1_km = states(1:3, :) * config.length_scale;
    run.position_2_km = states(4:6, :) * config.length_scale;
    run.center_km = ((1 - problem.eta) * states(1:3, :) + ...
                     problem.eta * states(4:6, :)) * config.length_scale;
    run.spin_ratio = spin / problem.spin_rate;
    run.length_error_m = length_error * config.length_scale * 1e3;
    run.rate_error_mps = rate_error / problem.rate_scale;
    run.tension_mps2 = tension * problem.accel_scale;

    % Tension = feedforward (free accelerations and centripetal) + feedback.
    spring = case_args.stiffness * length_error;
    damper = case_args.damping * rate_error;
    run.spring_mps2 = spring * problem.accel_scale;
    run.damper_mps2 = damper * problem.accel_scale;
    run.feedforward_mps2 = (tension - spring - damper) * problem.accel_scale;
    run.correction_mps2 = run.spring_mps2 + run.damper_mps2;

    initial_error = run.length_error_m(1);
    same_side = run.length_error_m * sign(initial_error);
    run.overshoot_percent = 100 * max(0, -min(same_side)) / abs(initial_error);
    run.final_error_m = run.length_error_m(end);
    run.final_rate_error_mps = run.rate_error_mps(end);
    run.peak_error_m = max(abs(run.length_error_m));
    run.peak_rate_error_mps = max(abs(run.rate_error_mps));
    run.peak_tension_mps2 = max(abs(run.tension_mps2));
    run.peak_correction_mps2 = max(abs(run.correction_mps2));
    run.settle_hours = settleTime(run.hours, run.length_error_m, ...
                                  config.settle_tolerance_m);
end

function settle = settleTime(hours, error_m, tolerance)
    % First time after which the error stays inside the tolerance (NaN if never).
    outside = find(abs(error_m) > tolerance, 1, 'last');
    if isempty(outside)
        settle = 0;
    elseif outside == numel(hours)
        settle = NaN;
    else
        settle = hours(outside + 1);
    end
end

function error_m = theoryError(run, omega_n_per_s)
    % Exact length error of e_ddot + c * e_dot + k * e = 0 from the sampled
    % initial error. A zero omega_n_per_s means no spring and no damper.
    seconds = run.hours * 3600;
    e_0 = run.length_error_m(1);
    e_dot_0 = run.rate_error_mps(1);
    if omega_n_per_s == 0
        error_m = e_0 + e_dot_0 * seconds;
    else
        error_m = (e_0 + (e_dot_0 + omega_n_per_s * e_0) * seconds) .* ...
                  exp(-omega_n_per_s * seconds);
    end
end

%% Dynamics
function [state_derivative, tension] = springDamperDynamics(t, state, problem, ...
                                                            case_args)
    % Dynamics of [p1; p2; v1; v2] with the spring-damper tether.
    body_1 = [state(1:3); state(7:9)];
    body_2 = [state(4:6); state(10:12)];
    args = {problem.eta, problem.length, case_args.rate, case_args.stiffness, ...
            case_args.damping, problem.accel_fn};
    spring_args = [{t, body_1, body_2}, args];
    [acceleration, ~, ~, tension] = ...
        coa.orbit.tetherSpringDamperAccel(spring_args{:});
    state_derivative = [state(7:12); acceleration];
end

function state = initialState(config, problem, rate)
    % Tether state whose length and length rate are off from the reference.
    state = coa.orbit.tetherInitialState(config.halo_state, ...
                                         problem.length + problem.length_error, ...
                                         problem.eta, problem.spin_rate, ...
                                         rate + problem.rate_error, ...
                                         config.tether_direction);
end

%% Report
function printComparison(config, runs)
    % Side-by-side metrics for every Ldot.
    rates = config.length_rate_values_mps;
    omega = config.frequency_ratio * config.spin_rate;
    fprintf('\nInitial error: %+g km length, %+g m/s length rate\n', ...
            config.initial_length_error_km, config.initial_length_rate_error_mps);
    fprintf(['Critically damped: omega_n = %g * spin = %.4g rad/s, ' ...
             'k = %.4g 1/s^2, c = %.4g 1/s\n'], ...
            config.frequency_ratio, omega, omega^2, 2 * omega);
    fprintf('Settle tolerance: %g m\n\n', config.settle_tolerance_m);
    fprintf(['  Ldot   case          settle [h]  peak |e| [m]  final e [m]' ...
             '  final edot [m/s]  peak T [m/s^2]   steps  fevals  wall [s]\n']);

    max_theory_gap = 0;
    for r = 1:numel(rates)
        printRow(rates(r), 'uncorrected', runs.uncorrected{r});
        printRow(rates(r), 'critical', runs.critical{r});
        gap = max(abs(runs.critical{r}.length_error_m - ...
                      theoryError(runs.critical{r}, omega)));
        max_theory_gap = max(max_theory_gap, gap);
    end
    fprintf('\nMax |simulation - exact critical-damping theory|: %.3g m\n', ...
            max_theory_gap);
    fprintf('Critical-damping overshoot: %.3g %% at worst\n', ...
            max(cellfun(@(run) run.overshoot_percent, runs.critical)));
    failed = sum(~cellfun(@(run) run.completed, [runs.uncorrected; runs.critical]));
    fprintf('Integration failures: %d\n', failed);
end

function printRow(rate, name, run)
    fprintf(['  %+5.2f  %-12s  %10.3g  %12.4g  %12.4g  %16.4g  %14.4g' ...
             '  %6d  %6d  %8.3g\n'], ...
            rate, name, run.settle_hours, run.peak_error_m, run.final_error_m, ...
            run.final_rate_error_mps, run.peak_tension_mps2, run.steps, ...
            run.function_evaluations, run.wall_time);
end

%% Plots
function [uncorrected, critical, index] = selectedRuns(config, runs)
    % Runs for the Ldot drawn in single-case panels.
    index = find(config.length_rate_values_mps == config.plot_length_rate_mps, 1);
    uncorrected = runs.uncorrected{index};
    critical = runs.critical{index};
end

function style = plotStyle()
    % Shared colors: gray dashed for uncorrected, blue for critically damped.
    style.uncorrected = {'Color', [0.35, 0.35, 0.35], 'LineStyle', '--', ...
                         'LineWidth', 1.4};
    style.critical = {'Color', [0, 0.447, 0.741], 'LineWidth', 1.6};
    style.theory = {'k.', 'MarkerSize', 6};
end

function plotErrorHistories(config, problem, runs)
    % Length and length-rate error of both tethers against exact theory.
    [uncorrected, critical, ~] = selectedRuns(config, runs);
    style = plotStyle();
    omega = problem.omega_n / config.time_scale;
    early = critical.hours <= config.early_hours;
    rates = config.length_rate_values_mps;
    colors = lines(numel(rates));

    figure('Name', 'Length and rate error');
    layout = tiledlayout(2, 2, 'TileSpacing', 'compact');
    title(layout, sprintf('Prescribed Ldot = %g m/s (initial error %g km, %g m/s)', ...
                          config.plot_length_rate_mps, ...
                          config.initial_length_error_km, ...
                          config.initial_length_rate_error_mps));

    nexttile;
    plot(uncorrected.hours(early), uncorrected.length_error_m(early), ...
         style.uncorrected{:});
    hold on;
    plot(critical.hours(early), critical.length_error_m(early), style.critical{:});
    theory = theoryError(critical, omega);
    plot(critical.hours(early), theory(early), style.theory{:});
    xlabel('time [h]');
    ylabel('L - L_{ref} [m]');
    title('Length error, early (dots: exact theory)');
    legend('no spring/damper', 'critically damped', 'theory', ...
           'Location', 'east');
    grid on;

    nexttile;
    plot(uncorrected.hours(early), uncorrected.rate_error_mps(early), ...
         style.uncorrected{:});
    hold on;
    plot(critical.hours(early), critical.rate_error_mps(early), style.critical{:});
    xlabel('time [h]');
    ylabel('Ldot - Ldot_{ref} [m/s]');
    title('Length-rate error, early');
    grid on;

    nexttile;
    for r = 1:numel(rates)
        semilogy(runs.uncorrected{r}.hours, abs(runs.uncorrected{r}.length_error_m), ...
                 '--', 'Color', colors(r, :));
        hold on;
        semilogy(runs.critical{r}.hours, abs(runs.critical{r}.length_error_m), ...
                 'Color', colors(r, :), 'LineWidth', 1.4);
    end
    yline(config.settle_tolerance_m, ':', 'tolerance');
    ylim([1e-6, 1e5]);
    xlabel('time [h]');
    ylabel('|L - L_{ref}| [m]');
    title('Length error, full day (dashed: no spring/damper)');
    grid on;
    colormap(gca, colors);
    colorbar('Ticks', linspace(0, 1, numel(rates)), ...
             'TickLabels', compose('%+g', rates));

    nexttile;
    reference = problem.length * config.length_scale + ...
                config.plot_length_rate_mps * 1e-3 * critical.hours * 3600;
    plot(critical.hours, reference, 'k:', 'LineWidth', 1.5);
    hold on;
    plot(uncorrected.hours, uncorrected.length_km, style.uncorrected{:});
    plot(critical.hours, critical.length_km, style.critical{:});
    xlabel('time [h]');
    ylabel('separation [km]');
    title('Tether length against the prescription');
    legend('L_{ref}(t)', 'no spring/damper', 'critically damped', ...
           'Location', 'northwest');
    grid on;
end

function plotForces(config, runs)
    % Tension, its feedback decomposition, and the size of the correction.
    [uncorrected, critical, ~] = selectedRuns(config, runs);
    style = plotStyle();
    early = critical.hours <= config.early_hours;
    rates = config.length_rate_values_mps;
    colors = lines(numel(rates));

    figure('Name', 'Tension and correction');
    layout = tiledlayout(2, 2, 'TileSpacing', 'compact');
    title(layout, sprintf('Prescribed Ldot = %g m/s', config.plot_length_rate_mps));

    nexttile;
    plot(uncorrected.hours(early), uncorrected.tension_mps2(early), ...
         style.uncorrected{:});
    hold on;
    plot(critical.hours(early), critical.tension_mps2(early), style.critical{:});
    xlabel('time [h]');
    ylabel('tension [m/s^2]');
    title('Tension (acceleration along the tether), early');
    legend('no spring/damper', 'critically damped', 'Location', 'best');
    grid on;

    nexttile;
    plot(critical.hours(early), critical.feedforward_mps2(early), 'LineWidth', 1.3);
    hold on;
    plot(critical.hours(early), critical.spring_mps2(early), 'LineWidth', 1.3);
    plot(critical.hours(early), critical.damper_mps2(early), 'LineWidth', 1.3);
    plot(critical.hours(early), critical.tension_mps2(early), 'k', 'LineWidth', 1.6);
    xlabel('time [h]');
    ylabel('[m/s^2]');
    title('Critically damped tension split');
    legend('feedforward', 'spring -k e', 'damper -c e_{dot}', 'total', ...
           'Location', 'best');
    grid on;

    nexttile;
    semilogy(critical.hours, abs(critical.tension_mps2 - uncorrected.tension_mps2), ...
             'Color', style.critical{2}, 'LineWidth', 1.4);
    hold on;
    semilogy(critical.hours, abs(critical.tension_mps2), 'k:', 'LineWidth', 1.2);
    xlabel('time [h]');
    ylabel('[m/s^2]');
    title('Full day: |T_{crit} - T_{uncorr}| against |T_{crit}|');
    legend('difference', 'critically damped tension', 'Location', 'best');
    grid on;

    nexttile;
    for r = 1:numel(rates)
        semilogy(runs.critical{r}.hours, abs(runs.critical{r}.correction_mps2), ...
                 'Color', colors(r, :), 'LineWidth', 1.3);
        hold on;
    end
    xlabel('time [h]');
    ylabel('|spring + damper| [m/s^2]');
    title('Feedback effort for every Ldot');
    legend(compose('Ldot = %+g m/s', rates), 'Location', 'northeast');
    grid on;
end

function plotTrajectoryDifference(config, runs)
    % How far the correction moves the trajectory and the tether's rotation.
    [uncorrected, critical, ~] = selectedRuns(config, runs);
    style = plotStyle();
    rates = config.length_rate_values_mps;
    colors = lines(numel(rates));

    figure('Name', 'Trajectory difference');
    layout = tiledlayout(2, 2, 'TileSpacing', 'compact');
    title(layout, sprintf('Critically damped minus uncorrected, Ldot = %g m/s', ...
                          config.plot_length_rate_mps));

    nexttile;
    plot(critical.hours, vecnorm(critical.position_1_km - ...
                                 uncorrected.position_1_km), 'LineWidth', 1.4);
    hold on;
    plot(critical.hours, vecnorm(critical.position_2_km - ...
                                 uncorrected.position_2_km), 'LineWidth', 1.4);
    plot(critical.hours, vecnorm(critical.center_km - uncorrected.center_km), ...
         'k', 'LineWidth', 1.4);
    xlabel('time [h]');
    ylabel('position difference [km]');
    title('Mass and center-of-mass position');
    legend('m_1', 'm_2', 'center of mass', 'Location', 'northwest');
    grid on;

    nexttile;
    for r = 1:numel(rates)
        difference = runs.critical{r}.center_km - runs.uncorrected{r}.center_km;
        plot(runs.critical{r}.hours, vecnorm(difference), 'Color', colors(r, :), ...
             'LineWidth', 1.3);
        hold on;
    end
    xlabel('time [h]');
    ylabel('center-of-mass difference [km]');
    title('Center of mass for every Ldot');
    legend(compose('Ldot = %+g m/s', rates), 'Location', 'northwest');
    grid on;

    nexttile;
    plot(uncorrected.hours, uncorrected.spin_ratio, style.uncorrected{:});
    hold on;
    plot(critical.hours, critical.spin_ratio, style.critical{:});
    xlabel('time [h]');
    ylabel('spin rate / prescribed spin');
    title('Inertial spin rate about z');
    legend('no spring/damper', 'critically damped', 'Location', 'best');
    grid on;

    nexttile;
    cosine = sum(critical.separation_km .* uncorrected.separation_km) ./ ...
             (critical.length_km .* uncorrected.length_km);
    plot(critical.hours, acosd(min(1, max(-1, cosine))), 'LineWidth', 1.4);
    xlabel('time [h]');
    ylabel('angle [deg]');
    title('Angle between the two tether directions');
    grid on;
end

function plotGeometry(config, runs)
    % Separation vector in the x-y plane of the rotating frame.
    [uncorrected, critical, ~] = selectedRuns(config, runs);
    style = plotStyle();
    early = critical.hours <= config.early_hours;
    angle = linspace(0, 2 * pi, 361);
    reference_km = config.tether_length_km;

    figure('Name', 'Separation vector geometry');
    layout = tiledlayout(1, 2, 'TileSpacing', 'compact');
    title(layout, sprintf('Prescribed Ldot = %g m/s', config.plot_length_rate_mps));

    nexttile;
    plot(reference_km * cos(angle), reference_km * sin(angle), 'k:', ...
         'LineWidth', 1.2);
    hold on;
    plot(uncorrected.separation_km(1, early), uncorrected.separation_km(2, early), ...
         style.uncorrected{:});
    plot(critical.separation_km(1, early), critical.separation_km(2, early), ...
         style.critical{:});
    plot(uncorrected.separation_km(1, 1), uncorrected.separation_km(2, 1), ...
         'ko', 'MarkerFaceColor', 'k');
    axis equal;
    xlabel('x [km]');
    ylabel('y [km]');
    title(sprintf('First %g h (dotted: L_0)', config.early_hours));
    legend('L_0 circle', 'no spring/damper', 'critically damped', 'start', ...
           'Location', 'northeast');
    grid on;

    nexttile;
    late = critical.hours >= config.duration_s / 3600 - config.early_hours;
    plot(uncorrected.separation_km(1, late), uncorrected.separation_km(2, late), ...
         style.uncorrected{:});
    hold on;
    plot(critical.separation_km(1, late), critical.separation_km(2, late), ...
         style.critical{:});
    plot(reference_km * cos(angle), reference_km * sin(angle), 'k:', ...
         'LineWidth', 1.2);
    axis equal;
    xlabel('x [km]');
    ylabel('y [km]');
    title(sprintf('Last %g h (dotted: L_0, not L_{ref}(t))', config.early_hours));
    grid on;
end

function plotCost(config, runs)
    % Solver cost of the two tethers for every Ldot.
    rates = config.length_rate_values_mps;
    metrics = {'steps', 'function_evaluations', 'wall_time'};
    labels = {'ode113 steps', 'function evaluations', 'wall time [s]'};

    figure('Name', 'Solver cost');
    layout = tiledlayout(1, numel(metrics), 'TileSpacing', 'compact');
    title(layout, 'Solver cost (no spring/damper vs critically damped)');
    for m = 1:numel(metrics)
        nexttile;
        values = [cellfun(@(run) run.(metrics{m}), runs.uncorrected), ...
                  cellfun(@(run) run.(metrics{m}), runs.critical)];
        bar(rates, values);
        xlabel('Ldot [m/s]');
        ylabel(labels{m});
        grid on;
        if m == 1
            legend('no spring/damper', 'critically damped', 'Location', 'best');
        end
    end
end
