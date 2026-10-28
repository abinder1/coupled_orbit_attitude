% Copyright (c) 2026 Andrew Binder

% tetherCovarianceDemo
%   Demonstrates coa.orbit.tetheredThreeBodyEkfPropagation, the ODE right-hand
%   side for the continuous-time (prediction) step of an EKF, on a spinning
%   two-mass tether (bodies 1 and 2) plus a free payload (body 3) in the
%   Earth-Moon CR3BP. Equations: derivations/ekf_equations.tex.
%
%   What the script does, in order:
%     1. Builds an initial tether state with its center of mass (CM) on a halo
%        orbit, spinning about +z at a fixed length.
%     2. Builds the payload state so that the payload ARRIVES AT the catch
%        mechanism at the end of the span (the "rendezvous" construction):
%        propagate the tether forward, copy the catch mechanism's final state,
%        and back-propagate it as a free body to t = 0.
%     3. Propagates state and covariance together, twice (process noise on and
%        off), using the packed vector [x; P(:)] (18 + 324 = 342 elements).
%        P_dot = F*P + P*F' + G*Q*G'  (Eq. 5 of the derivation).
%     4. Re-propagates the state alone (18-element mode, no covariance) for
%        the point-mass visualization.
%     5. Sweeps the two tether process noises separately from zero initial P
%        (Figure 4).
%
%   Figures:
%     1. 'Uncertainty growth'  - 1-sigma position / velocity per body vs time,
%        process noise on (solid) and off (dashed).
%     2. 'Point masses'        - rotating-frame trajectories of the tether CM
%        and payload, and the three bodies relative to the tether CM.
%     3. 'Final correlation'   - correlation matrix of P at the final time
%        (process noise on).
%     4. 'Axial vs transverse noise' - check that tension_psd drives only the
%        along-tether direction and transverse_psd only the perpendicular
%        ones: sweeps both PSDs from zero initial P over a ~2 min span and
%        plots the relative-velocity 1-sigma split into axial and
%        perpendicular parts, against the analytic prediction (dashed).
%
%   State convention (nondimensional, CR3BP rotating frame, x-y-z axes):
%     x = [p1; p2; p3; v1; v2; v3]   (each pK, vK is 3x1)
%     body 1 = counterweight, 2 = catch mechanism, 3 = payload.
%     P is the 18x18 covariance in the same ordering, stored as P(:).
%
%   Units: everything passed to the propagator is nondimensional (length
%   L_star, time T_star). Plots convert to km, m/s, and hours.
%
%   Things to remember when revisiting:
%     - The initial covariance is diagonal (bodies uncorrelated), which does
%       not respect the tether constraints, and the spinning-tether dynamics
%       amplify it strongly. Expect large uncertainty growth for bodies 1
%       and 2; the payload (free body) grows slowly.
%     - The covariance entries are tiny (~1e-14), so ODE AbsTol must be
%       scaled separately for the covariance block.
%     - The tiny negative min eigenvalue printed in the summary is integrator
%       roundoff, not a bug.
%     - Process-noise PSDs are placeholders chosen so the effect is visible.

%% Script setup
clear;
clc;
close('all');
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'toolbox'));

%% System parameters
% Earth-Moon CR3BP characteristic scales
L_star = 384400;            % Length unit [km]
T_star = 375190.258892627;  % Time unit [s]
mu = 0.012150585;           % Moon's mass fraction

% Tether, spun about +z at a fixed length
tether_length_km = 200;
spin_rate = 2 * pi / (55 * 60);  % Inertial spin rate [rad/s]
mass_2_fraction = 0.5;           % m2 / (m1 + m2); 0.5 makes bodies 1 and 2 symmetric
tether_length_rate = 0;          % Prescribed length rate (nondimensional)
lhat_0 = [1; 0; 0];              % Initial unit vector from the CM toward body 1

% Initial tether CM state [position; velocity] on a halo orbit (rotating frame,
% nondimensional)
halo_ic = [0.831987145310294; 0; 0.12575552414079; 0; 0.239887839725451; 0];

% Propagation span (kept short; sensitivity of these dynamics is large)
prop_hours = 24;
t_end = prop_hours * 3600 / T_star;  % Nondimensional

% Initial 1-sigma uncertainties, applied independently to every body
initial_position_sigma_km = 0.01;
initial_velocity_sigma_ms = 1e-4;

% Process-noise spectral densities (nondimensional), see the propagator header:
%   tension_psd             radial (tension) noise on the tether, Eq. 20
%   transverse_psd          transverse noise on the tether, Eq. 26
%   body_3_acceleration_psd white acceleration noise on the payload
process_noise = struct('tension_psd', 1e-8, ...
                       'transverse_psd', 1e-8, ...
                       'body_3_acceleration_psd', 1e-8);
no_process_noise = structfun(@(~) 0, process_noise, 'UniformOutput', false);

% Free-body dynamics: [a, da_dp, da_dv] = cr3bp(state, mu), state = [p; v]
accel_fn = @(state) coa.orbit.cr3bp(state, mu);

% Short alias for the propagator (keeps the ODE wrappers within 88 columns)
propagate = @coa.orbit.tetheredThreeBodyEkfPropagation;

%% Initial state
% Tether bodies 1 and 2 at t = 0. Layout: [p1; p2; v1; v2] (12x1).
terminals_state = tetherInitialState(halo_ic, ...
                                     tether_length_km / L_star, ...
                                     mass_2_fraction, spin_rate * T_star, ...
                                     tether_length_rate, lhat_0);

% Payload (rendezvous construction):
%   (1) propagate the tether forward to t_end, state-only,
%   (2) take the catch mechanism's final state as the payload's final state,
%   (3) back-propagate that state to t = 0 as a free body.
% The payload therefore meets the catch mechanism at t_end.
%
% The propagator needs a full 18-element state, so body 3 is given a
% placeholder copy of the catch mechanism. It does not influence bodies 1 and
% 2 and is discarded.
tether_only_rhs = @(t, state) propagate(t, state, mass_2_fraction, ...
                                        tether_length_rate, accel_fn);
tether_only_state_0 = [terminals_state(1:6); terminals_state(4:6); ...
                       terminals_state(7:12); terminals_state(10:12)];
tether_only_sol = ode113(tether_only_rhs, [0, t_end], tether_only_state_0, ...
                         odeset('RelTol', 1e-10, 'AbsTol', 1e-12));
tether_final_state = tether_only_sol.y(:, end);

% Catch mechanism (body 2) final [p2; v2] from the [p1;p2;p3;v1;v2;v3] layout
catch_mechanism_final = [tether_final_state(4:6); tether_final_state(13:15)];
payload_state = propagateSingleBody(catch_mechanism_final, -t_end, accel_fn);

% Assemble [p1; p2; v1; v2; p3; v3], then re-order to the propagator's
% [p1; p2; p3; v1; v2; v3] layout.
initial_state = [terminals_state; payload_state];
initial_state = [initial_state(1:6); initial_state(13:15); ...
                 initial_state(7:12); initial_state(16:18)];

%% Initial covariance
% Diagonal: every position and velocity component independent, with equal
% sigma per component (nondimensionalized from the km and m/s inputs).
position_sigma = initial_position_sigma_km / L_star;
velocity_sigma = initial_velocity_sigma_ms / 1000 * T_star / L_star;
initial_covariance = diag([repmat(position_sigma^2, 9, 1); ...
                           repmat(velocity_sigma^2, 9, 1)]);

%% Propagate the state and covariance, with and without process noise
% Output times for the stored solution (deval interpolates the ODE solution)
t_eval = linspace(0, t_end, 1000);

% The covariance entries are ~1e-14, so the tolerances need separate scales:
% AbsTol is 1e-10 for the 18 states and 1e-24 for the 324 covariance entries.
ode_opts = odeset('RelTol', 1e-9, ...
                  'AbsTol', [1e-10 * ones(18, 1); 1e-24 * ones(324, 1)]);

% results(k).packed is 342 x numel(t_eval): rows 1:18 are the state, rows
% 19:342 are P(:) (reshape to 18x18 to use).
results = struct('name', {"Process noise on", "Process noise off"}, ...
                 'noise', {process_noise, no_process_noise});
for k = 1:numel(results)
    rhs = @(t, packed_state) propagate(t, packed_state, mass_2_fraction, ...
                                       tether_length_rate, accel_fn, ...
                                       results(k).noise);

    timer = tic;
    sol = ode113(rhs, [0, t_end], [initial_state; initial_covariance(:)], ...
                 ode_opts);
    results(k).elapsed_s = toc(timer);
    results(k).packed = deval(sol, t_eval);
    results(k).n_steps = numel(sol.x);
end

%% Summary
% Sanity check of the final covariance: min/max eigenvalue ratio should be ~0
% (slightly negative is roundoff) and the asymmetry should be ~0.
fprintf('Propagated %.1f h (%.4f nondimensional time units)\n', ...
        prop_hours, t_end);
for k = 1:numel(results)
    final_covariance = reshape(results(k).packed(19:end, end), [18, 18]);
    eigenvalues = eig((final_covariance + final_covariance') / 2);
    fprintf(['  %-18s %5d steps, %5.2f s, min/max eig(P) = %+.1e, ', ...
             'asymmetry = %.1e\n'], ...
            results(k).name, results(k).n_steps, results(k).elapsed_s, ...
            min(eigenvalues) / max(eigenvalues), ...
            norm(final_covariance - final_covariance', 'fro'));
end

%% Figure 1: uncertainty growth
% Left panel is position (km), right panel is velocity (m/s). Each curve is
% sqrt(trace) of the 3x3 diagonal block of P belonging to one body.
t_hours = t_eval * T_star / 3600;
body_names = ["Counterweight (1)", "Catch mechanism (2)", "Payload (3)"];
body_colors = [0.835, 0.369, 0.000; ...  % vermillion
               0.000, 0.447, 0.698; ...  % blue
               0.000, 0.619, 0.451];     % green

figure('Name', 'Uncertainty growth', 'Color', 'w');
line_styles = ["-", "--"];
% Bodies 1 and 2 are symmetric when mass_2_fraction = 0.5, so draw body 1
% wider to keep it visible underneath body 2.
line_widths = [4, 1.5, 1.5];
tiledlayout(1, 2, 'TileSpacing', 'compact');
for panel = 1:2
    nexttile;
    hold('on');
    for body = 1:3
        % rows: indices of this body's position (1:9) or velocity (10:18)
        if panel == 1
            rows = 3 * (body - 1) + (1:3);
            scale = L_star;                           % km
        else
            rows = 9 + 3 * (body - 1) + (1:3);
            scale = L_star / T_star * 1000;           % m/s
        end
        for k = 1:numel(results)
            sigma = blockSigma(results(k).packed, rows) * scale;
            plot(t_hours, sigma, 'Color', body_colors(body, :), ...
                 'LineStyle', line_styles(k), 'LineWidth', line_widths(body), ...
                 'DisplayName', body_names(body) + " - " + results(k).name);
        end
    end
    set(gca, 'YScale', 'log');
    grid('on');
    xlabel('Time [h]');
    if panel == 1
        ylabel('Position 1$\sigma$, $\sqrt{\mathrm{tr}\,P_{pp}}$ [km]', ...
               'Interpreter', 'latex');
        title('Position uncertainty');
        legend('Location', 'northwest', 'FontSize', 7);
    else
        ylabel('Velocity 1$\sigma$, $\sqrt{\mathrm{tr}\,P_{vv}}$ [m/s]', ...
               'Interpreter', 'latex');
        title('Velocity uncertainty');
    end
end

%% Figure 2: point-mass visualization (state-only propagation)
% An 18-element packed state selects state-only mode (no covariance, no
% process_noise needed). Tighter tolerances than above since no covariance
% is carried.
t_plot = linspace(0, t_end, 1500);
state_rhs = @(t, state) propagate(t, state, mass_2_fraction, ...
                                  tether_length_rate, accel_fn);
state_sol = ode113(state_rhs, [0, t_end], initial_state, ...
                   odeset('RelTol', 1e-10, 'AbsTol', 1e-12));
states = deval(state_sol, t_plot);
p1 = states(1:3, :) * L_star;   % [km]
p2 = states(4:6, :) * L_star;
p3 = states(7:9, :) * L_star;
% Tether CM is the mass-weighted mean (body 2 carries mass_2_fraction)
p_com = (1 - mass_2_fraction) * p1 + mass_2_fraction * p2;

figure('Name', 'Point masses', 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact');

% Left: rotating-frame trajectories (the tether is far too small to resolve)
nexttile;
hold('on');
plot3(p_com(1, :), p_com(2, :), p_com(3, :), 'Color', body_colors(1, :), ...
      'LineWidth', 1.5, 'DisplayName', 'Tether CM');
plot3(p3(1, :), p3(2, :), p3(3, :), 'Color', body_colors(3, :), ...
      'LineWidth', 1.5, 'DisplayName', 'Payload (3)');
plot3(p3(1, end), p3(2, end), p3(3, end), 'o', ...
      'Color', body_colors(3, :), 'MarkerFaceColor', body_colors(3, :), ...
      'HandleVisibility', 'off');
plot3(p_com(1, end), p_com(2, end), p_com(3, end), 'o', ...
      'Color', body_colors(1, :), 'MarkerFaceColor', body_colors(1, :), ...
      'HandleVisibility', 'off');
grid('on');
axis('equal');
view([0, 90]);
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
title('Rotating-frame trajectories');
legend('Location', 'best');

% Right: bodies relative to the tether CM, with tether snapshots
nexttile;
hold('on');
rel_1 = p1 - p_com;
rel_2 = p2 - p_com;
rel_3 = p3 - p_com;
plot3(rel_3(1, :), rel_3(2, :), rel_3(3, :), 'Color', body_colors(3, :), ...
      'DisplayName', 'Payload (3)');
snapshots = round(linspace(1, numel(t_plot), 12));
for idx = snapshots
    plot3([rel_1(1, idx), rel_2(1, idx)], [rel_1(2, idx), rel_2(2, idx)], ...
          [rel_1(3, idx), rel_2(3, idx)], '-', 'Color', [0.6, 0.6, 0.6], ...
          'HandleVisibility', 'off');
end
plot3(rel_1(1, snapshots), rel_1(2, snapshots), rel_1(3, snapshots), 'o', ...
      'Color', body_colors(1, :), 'MarkerFaceColor', body_colors(1, :), ...
      'DisplayName', 'Counterweight (1)');
plot3(rel_2(1, snapshots), rel_2(2, snapshots), rel_2(3, snapshots), 's', ...
      'Color', body_colors(2, :), 'MarkerFaceColor', body_colors(2, :), ...
      'DisplayName', 'Catch mechanism (2)');
plot3(0, 0, 0, 'k+', 'DisplayName', 'Tether CM');
plot3(rel_3(1, end), rel_3(2, end), rel_3(3, end), 'o', ...
      'Color', body_colors(3, :), 'MarkerFaceColor', body_colors(3, :), ...
      'HandleVisibility', 'off');
grid('on');
axis('equal');
view([0, 90]);
xlabel('\Deltax [km]');
ylabel('\Deltay [km]');
zlabel('\Deltaz [km]');
title('Relative to tether CM (snapshots)');
legend('Location', 'best');

%% Figure 3: final correlation matrix
% corr_ij = P_ij / sqrt(P_ii * P_jj), process-noise-on run. Axis blocks are
% the 3 components of p1, p2, p3, v1, v2, v3 in order.
final_covariance = reshape(results(1).packed(19:end, end), [18, 18]);
correlation = final_covariance ./ sqrt(diag(final_covariance) * ...
                                       diag(final_covariance)');
figure('Name', 'Final correlation', 'Color', 'w');
imagesc(correlation, [-1, 1]);
axis('image');
colorbar;
colormap(parula);
labels = ["p1", "p2", "p3", "v1", "v2", "v3"];
set(gca, 'XTick', 2:3:18, 'XTickLabel', labels, ...
    'YTick', 2:3:18, 'YTickLabel', labels);
title(sprintf('Correlation of P after %.0f h (process noise on)', prop_hours));

%% Figure 4: axial vs. transverse process-noise check
% Goal: confirm tension_psd feeds only the along-tether direction and
% transverse_psd feeds only the two perpendicular directions.
%
% Method: start from ZERO covariance so P contains only injected process noise
% (no initial-uncertainty amplification), run a short span, and decompose the
% covariance of the relative velocity dv = v2 - v1 into the component along
% the tether direction e = (p2 - p1)/|p2 - p1| (axial) and the remaining two
% (perpendicular, reported as the sum). Noise enters as acceleration, so
% velocity is the cleanest metric; position mixes the directions over time as
% the tether spins.
%
% Prediction: relative acceleration noise has gain g = 1/(eta*(1 - eta)),
% where eta = mass_2_fraction, so for small t
%   axial var = tension_psd * g^2 * t
%   perp  var = 2 * transverse_psd * g^2 * t  (2 perpendicular directions)
% These are drawn as dashed lines; the span is short enough (~13 deg of tether
% spin) that the spin-induced mixing between the two is small.
sweep_hours = 2 / 60;
sweep_end = sweep_hours * 3600 / T_star;
sweep_times = linspace(0, sweep_end, 200);
base_psd = 1e-8;

% [tension multiplier, transverse multiplier] applied to base_psd
sweep_cases = [1, 0; 0, 1; 1, 1; 4, 1; 1, 4];
gain = 1 / (mass_2_fraction * (1 - mass_2_fraction));
sweep_ode_opts = odeset('RelTol', 1e-9, ...
                        'AbsTol', [1e-10 * ones(18, 1); 1e-24 * ones(324, 1)]);

n_cases = size(sweep_cases, 1);
axial_sigma = zeros(n_cases, numel(sweep_times));
perp_sigma = zeros(n_cases, numel(sweep_times));
case_labels = strings(n_cases, 1);
relative_velocity_map = [-eye(3), eye(3)];  % dv = v2 - v1 from [v1; v2]
for c = 1:n_cases
    sweep_noise = struct('tension_psd', sweep_cases(c, 1) * base_psd, ...
                         'transverse_psd', sweep_cases(c, 2) * base_psd, ...
                         'body_3_acceleration_psd', 0);
    rhs = @(t, packed_state) propagate(t, packed_state, mass_2_fraction, ...
                                       tether_length_rate, accel_fn, ...
                                       sweep_noise);
    sol = ode113(rhs, [0, sweep_end], [initial_state; zeros(324, 1)], ...
                 sweep_ode_opts);
    packed = deval(sol, sweep_times);

    for k = 1:numel(sweep_times)
        separation = packed(4:6, k) - packed(1:3, k);
        tether_direction = separation / norm(separation);
        covariance = reshape(packed(19:end, k), [18, 18]);
        relative_velocity_cov = relative_velocity_map * ...
                                covariance(10:15, 10:15) * ...
                                relative_velocity_map';
        axial_variance = tether_direction' * relative_velocity_cov * ...
                         tether_direction;
        axial_sigma(c, k) = sqrt(max(axial_variance, 0));
        perp_sigma(c, k) = sqrt(max(trace(relative_velocity_cov) - ...
                                    axial_variance, 0));
    end
    case_labels(c) = sprintf('q_T = %g, q_\\perp = %g (x base)', ...
                             sweep_cases(c, 1), sweep_cases(c, 2));
end

velocity_scale = L_star / T_star * 1000;  % nondimensional -> m/s
sweep_minutes = sweep_times * T_star / 60;
case_colors = lines(n_cases);

figure('Name', 'Axial vs transverse noise', 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact');
panel_titles = ["Axial (along tether)", "Perpendicular (2 directions, summed)"];
panel_data = {axial_sigma, perp_sigma};
for panel = 1:2
    nexttile;
    hold('on');
    for c = 1:n_cases
        plot(sweep_minutes, panel_data{panel}(c, :) * velocity_scale, ...
             'Color', case_colors(c, :), 'LineWidth', 1.5, ...
             'DisplayName', case_labels(c));
        if panel == 1
            predicted_variance = sweep_cases(c, 1) * base_psd * gain^2 * ...
                                 sweep_times;
        else
            predicted_variance = 2 * sweep_cases(c, 2) * base_psd * gain^2 * ...
                                 sweep_times;
        end
        plot(sweep_minutes, sqrt(predicted_variance) * velocity_scale, '--', ...
             'Color', case_colors(c, :), 'HandleVisibility', 'off');
    end
    set(gca, 'YScale', 'log');
    grid('on');
    xlabel('Time [min]');
    ylabel('Relative velocity 1\sigma [m/s]');
    title(panel_titles(panel));
    if panel == 1
        legend('Location', 'southeast', 'FontSize', 7);
    end
end
sgtitle('Process-noise decomposition (zero initial P; dashed = prediction)');

%% Local functions
function sigma = blockSigma(packed, rows)
    % sqrt(trace) of the sub-block of P over the given state rows, per column.
    n_times = size(packed, 2);
    sigma = zeros(1, n_times);
    for k = 1:n_times
        covariance = reshape(packed(19:end, k), [18, 18]);
        sigma(k) = sqrt(trace(covariance(rows, rows)));
    end
end

function state = propagateSingleBody(state_0, duration, accel_fn)
    % Coast a single body [p; v] under accel_fn for 'duration' (nondimensional;
    % negative propagates backward) and return the final [p; v].
    rhs = @(~, x) [x(4:6); accel_fn(x)];
    sol = ode113(rhs, [0, duration], state_0, ...
                 odeset('AbsTol', 1e-10, 'RelTol', 1e-10));
    state = sol.y(:, end);
end

function state = tetherInitialState(com_state, length_nd, mass_2_fraction, ...
                                    spin_rate_nd, length_dot, lhat_0)
    % Initial [p1; p2; v1; v2] for a tether whose center of mass starts at
    % com_state = [p; v] in the rotating frame (rigid spin plus extension).
    % spin_rate_nd is the inertial spin rate about +z; the CR3BP frame itself
    % rotates at 1 (nondimensional) about +z, so the rotating-frame spin is
    % the difference. Body 1 sits at +lhat_0 and body 2 at -lhat_0, with the
    % CM-relative offsets weighted by the mass fractions.
    omega_frame = [0; 0; 1];
    spin_rel = spin_rate_nd * [0; 0; 1] - omega_frame;

    rel_1 = lhat_0 * (mass_2_fraction * length_nd);
    rel_2 = -lhat_0 * ((1 - mass_2_fraction) * length_nd);

    vel_1 = cross(spin_rel, rel_1) + mass_2_fraction * length_dot * lhat_0;
    vel_2 = cross(spin_rel, rel_2) - (1 - mass_2_fraction) * length_dot * lhat_0;

    state = [com_state(1:3) + rel_1; com_state(1:3) + rel_2; ...
             com_state(4:6) + vel_1; com_state(4:6) + vel_2];
end
