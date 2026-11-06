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

%% Setup
clear;
clc;
close('all');
example_folder = fileparts(mfilename('fullpath'));
addpath(example_folder);
addpath(fullfile(example_folder, '..', 'toolbox'));
import tetherFilterDemoSupport.*
rng(1, 'twister');

config = demoConfig();
PlotOutput.preparePlotOutput(config);

%% Initialize and propagate the shared truth
[initial_state, initial_covariance, initial_factor] = ...
    Scenario.makeInitialEstimate(config);
scenario = Scenario.propagateScenario(initial_state, config);
Visualization.plotTruthTrajectory(scenario, config);

%% Run and display one filter realization
% Truth is the nominal initial condition; each trial's estimate is a sample
% drawn around it with the initial covariance.
sampled_state = Scenario.perturbInitialState(initial_state, initial_factor);
filter_run = Filter.runFilter(scenario, sampled_state, initial_covariance, config);
Visualization.reportFilterRun(filter_run, config);
Visualization.plotPredictionHistory(filter_run.history, config);

%% Repeat the filter and analyze its consistency
fprintf(['\nReference filter run complete. Starting %d Monte Carlo ' ...
         'runs...\n'], config.monte_carlo_runs);
monte_carlo = MonteCarlo.runMonteCarlo(scenario, initial_state, initial_covariance, ...
                                       initial_factor, filter_run, config);
fprintf('Monte Carlo runs complete. Generating reports and visualizations...\n');
Analysis.runMonteCarloAnalyses(Analysis.monteCarloAnalyses(), monte_carlo, config);

%% Scenario configuration
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
    config.plot_modes = struct('prediction_history', "show", ...
                               'trajectory', "show", ...
                               'anees', "show", ...
                               'rmse', "show", ...
                               'innovations', "show", ...
                               'residuals', "show", ...
                               'corner', "off", ...
                               'relative_corner', "show");
    config.plot_output_folder = fullfile(fileparts(mfilename('fullpath')), ...
                                         '..', 'runs', 'tetherFilterDemo');
    config.plot_file_format = "png";
    config.plot_resolution = 150;             % dpi for saved images

    config.position_rows = 1:9;
    config.velocity_rows = 10:18;
    config.body_names = ["Counterweight (1)", ...
                         "Catch mechanism (2)", "Payload (3)"];
end
