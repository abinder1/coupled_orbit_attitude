% Copyright (c) 2026 Andrew Binder

% tetherCovarianceDemo
%   Demonstrates coa.orbit.tetheredThreeBodyEkfPropagation, the ODE right-hand
%   side for the continuous-time (prediction) step of an EKF, on a spinning
%   two-mass tether (bodies 1 and 2) plus a free payload (body 3) in the
%   Earth-Moon CR3BP. Equations: derivations/ekf_equations.tex.
%
%   The tether is the critically damped spring-damper of
%   coa.orbit.tetherSpringDamperAccel (k = wn^2, c = 2 * wn, wn = spin rate),
%   which pulls |p2 - p1| toward L(t) = L_0 + Ldot * t. The rigid tether
%   (coa.orbit.tetherAccel) appears only in two reference checks.
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
%        P_dot = F*P + P*F' + G*Q*G'  (Eq. 5 of the derivation), with F from
%        the analytic spring-damper Jacobian.
%     4. Re-propagates the state alone (18-element mode, no covariance) for
%        the point-mass visualization.
%     5. Sweeps the two tether process noises separately from zero initial P
%        (Figure 4).
%     6. Starts the tether with a broken length constraint to show the
%        spring-damper recovering it in both the state and P, with two checks
%        against the rigid tether.
%     7. Checks the linearized covariance, including the length variance,
%        against a Monte Carlo run of the nonlinear dynamics at zero length
%        rate.
%     8. Repeats the covariance/Monte Carlo check while reeling in and out
%        at nonzero length rates.
%     9. Sweeps tension process noise and compares steady-state tether length
%        variance with the stationary-oscillator prediction.
%
%   Figures:
%     1. 'Uncertainty growth'  - 1-sigma position / velocity per body vs time,
%        process noise on (solid) and off (dashed).
%     2. 'Point masses'        - rotating-frame trajectories of the tether CM
%        and payload, and the three bodies relative to the tether CM.
%     3. 'Final correlation'   - correlation matrix of P at the final time
%        (process noise on).
%     4. 'Axial vs transverse noise' - sweeps both PSDs from zero initial P
%        over a ~2 min span and plots the relative-velocity 1-sigma split into
%        axial and perpendicular parts. The dashed prediction is the
%        rigid-tether result; the spring-damper's damper bleeds axial noise,
%        so a gap in the axial panel is expected.
%     5. 'Spring-damper constraint recovery' - length error, length-rate error,
%        and the length 1-sigma implied by P from a broken start (check 1,
%        rigid in gray), plus the state difference from the rigid tether
%        when starting on the constraint (check 2).
%     6. 'Linear covariance vs Monte Carlo' - linearized P against the sample
%        covariance of 100 nonlinear state-only runs, per body and for the
%        tether length.
%     7. 'Nonzero length rate' - linearized P against Monte Carlo for reel-in
%        and reel-out cases, with nominal length and rate tracking.
%     8. 'Steady-state length covariance' - tension-noise sweep, full
%        covariance propagation vs the stationary-oscillator prediction.
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
%     - The initial covariance is built from a loose COM and a tight relative
%       state (coa.filter.tetherInitialCovariance), so it is nearly consistent
%       with the tether constraint. The old diagonal form (every body
%       independent) put a length error into P0 that the spin dynamics
%       turned into a large spin-phase error.
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
prop_hours = 1 * 24;
t_end = prop_hours * 3600 / T_star;  % Nondimensional

% Initial 1-sigma uncertainties. The tether pair is described by its center of
% mass (COM) and its relative state: a loose COM with a tight relative state
% keeps the errors almost parallel to the tether constraint.
initial_sigma = struct('com_position', 100, ...          % [km]
                       'com_velocity', 1, ...            % [m/s]
                       'relative_position', 10 / 1000, ...    % [km]
                       'relative_velocity', 1, ...    % [m/s]
                       'payload_position', 100, ...      % [km]
                       'payload_velocity', 1);           % [m/s]

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

% Critically damped spring-damper tether, k = wn^2 and c = 2 * wn, with wn
% equal to the spin rate. The length reference is
% L(t) = L_0 + tether_length_rate * t, with t measured from the start.
natural_frequency = spin_rate * T_star;
spring_damper = struct('initial_length', tether_length_km / L_star, ...
                       'stiffness', natural_frequency^2, ...
                       'damping', 2 * natural_frequency);
rigid_model = struct();  % An empty struct selects the rigid tether

%% Initial state
% Tether bodies 1 and 2 at t = 0. Layout: [p1; p2; v1; v2] (12x1).
terminals_state = coa.orbit.tetherInitialState(halo_ic, ...
                                               tether_length_km / L_star, ...
                                               mass_2_fraction, ...
                                               spin_rate * T_star, ...
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
                                        tether_length_rate, accel_fn, ...
                                        struct(), spring_damper);
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
% P0 = A * A' with independent isotropic errors in the tether COM, the
% relative state, and the payload (see coa.filter.tetherInitialCovariance).
% The sigmas are converted from km and m/s to nondimensional units.
position_to_nd = 1 / L_star;
velocity_to_nd = T_star / L_star / 1000;
nd_sigma = initial_sigma;
for field = ["com_position", "relative_position", "payload_position"]
    nd_sigma.(field) = initial_sigma.(field) * position_to_nd;
end
for field = ["com_velocity", "relative_velocity", "payload_velocity"]
    nd_sigma.(field) = initial_sigma.(field) * velocity_to_nd;
end
[initial_covariance, initial_factor] = ...
    coa.filter.tetherInitialCovariance(mass_2_fraction, nd_sigma);

% Reference only: the old diagonal P0 with the COM sigma on every component.
% It ignores the constraint, so it carries a large length error.
diagonal_covariance = diag([repmat(nd_sigma.com_position^2, 9, 1); ...
                            repmat(nd_sigma.com_velocity^2, 9, 1)]);

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
                                       results(k).noise, spring_damper);

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
sgtitle('Spring-damper tether, 1\sigma per body (solid: noise on, dashed: off)');

%% Figure 2: point-mass visualization (state-only propagation)
% An 18-element packed state selects state-only mode (no covariance, no
% process_noise needed). Tighter tolerances than above since no covariance
% is carried.
t_plot = linspace(0, t_end, 1500);
state_rhs = @(t, state) propagate(t, state, mass_2_fraction, ...
                                  tether_length_rate, accel_fn, struct(), ...
                                  spring_damper);
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
legend('Location', 'northwest');

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
legend('Location', 'southwest');

%% Figure 3: final correlation matrix
% corr_ij = P_ij / sqrt(P_ii * P_jj), process-noise-on run. Axis blocks are
% the 3 components of p1, p2, p3, v1, v2, v3 in order.
final_covariance = reshape(results(1).packed(19:end, end), [18, 18]);
correlation = final_covariance ./ sqrt(diag(final_covariance) * ...
                                       diag(final_covariance)');
figure('Name', 'Final correlation', 'Color', 'w');
imagesc(correlation, [-1, 1]);
axis('image');
colorbar_handle = colorbar;
colorbar_handle.Label.String = 'Correlation coefficient';
colormap(parula);
xlabel('State block (3 components each)');
ylabel('State block (3 components each)');
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
                                       sweep_noise, spring_damper);
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
    case_labels(c) = sprintf('q_{axial} = %g, q_{perp} = %g (x 1e-8)', ...
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
             'Color', case_colors(c, :), 'LineWidth', 4 - 0.7 * (c - 1), ...
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
        plot(NaN, NaN, 'k--', 'DisplayName', 'rigid-tether prediction');
        legend('Location', 'southeast', 'FontSize', 7);
    end
end
sgtitle(['Process-noise decomposition, spring-damper tether (zero initial P; ' ...
         'dashed = rigid-tether prediction)']);

%% Constraint recovery and rigid-tether reference checks
% The spring-damper tether is the model above. Here it is started with a
% broken constraint (extra length and length rate relative to L(t)) and
% compared with the rigid tether (tetherAccel) in two ways:
%   Check 1 (broken start): the rigid tether keeps the length error forever,
%     and P keeps the length variance, while the spring-damper removes both.
%   Check 2 (on-constraint start): with the constraint already satisfied the
%     two models should give nearly the same state.
broken_length_error_km = 1;
broken_rate_error_ms = 0.01;

broken_rate_nd = tether_length_rate + broken_rate_error_ms / 1000 * T_star / L_star;
broken_terminals = coa.orbit.tetherInitialState(halo_ic, ...
                                                (tether_length_km + ...
                                                 broken_length_error_km) / L_star, ...
                                                mass_2_fraction, spin_rate * T_star, ...
                                                broken_rate_nd, lhat_0);
broken_state = [broken_terminals(1:6); initial_state(7:9); ...
                broken_terminals(7:12); initial_state(16:18)];

tether_runs = struct('name', {"Spring-damper", "Spring-damper, broken start", ...
                              "Rigid", "Rigid, broken start"}, ...
                     'state', {initial_state, broken_state, ...
                               initial_state, broken_state}, ...
                     'model', {spring_damper, spring_damper, ...
                               rigid_model, rigid_model});
for k = 1:numel(tether_runs)
    run_args = {mass_2_fraction, tether_length_rate, accel_fn, ...
                process_noise, tether_runs(k).model};
    timer = tic;
    packed = propagatePacked(tether_runs(k).state, initial_covariance, ...
                             t_end, t_eval, run_args, ode_opts);
    tether_runs(k).elapsed_s = toc(timer);
    tether_runs(k).packed = packed;
    reference_length = spring_damper.initial_length + tether_length_rate * t_eval;
    tether_runs(k).history = lengthHistory(packed, reference_length, ...
                                           tether_length_rate, ode_opts.RelTol);
end

fprintf('\nSpring-damper (wn = %.3g * spin rate) vs rigid, process noise on\n', ...
        natural_frequency / (spin_rate * T_star));
for k = 1:numel(tether_runs)
    final_covariance = reshape(tether_runs(k).packed(19:end, end), [18, 18]);
    eigenvalues = eig((final_covariance + final_covariance') / 2);
    history = tether_runs(k).history;
    fprintf(['  %-28s %5.2f s, min/max eig(P) = %+.1e, final L error = ' ...
             '%+.3g m, final L sigma = %.3g m\n'], ...
            tether_runs(k).name, tether_runs(k).elapsed_s, ...
            min(eigenvalues) / max(eigenvalues), ...
            history.length_error(end) * L_star * 1000, ...
            history.length_sigma(end) * L_star * 1000);
end

% Check 2: on the constraint, the two tethers should agree closely.
state_difference = tether_runs(1).packed(1:18, :) - tether_runs(3).packed(1:18, :);
fprintf('  Check 2, max |state(SD) - state(rigid)|: %.3g km, %.3g m/s\n', ...
        max(abs(state_difference(1:9, :)), [], 'all') * L_star, ...
        max(abs(state_difference(10:18, :)), [], 'all') * velocity_scale);

gray = [0.5, 0.5, 0.5];
run_colors = [0.000, 0.447, 0.698; 0.835, 0.369, 0.000; gray; gray];
run_styles = ["-", "-", "-", "--"];
length_scale_m = L_star * 1000;

%% Figure 5: constraint recovery
figure('Name', 'Spring-damper constraint recovery', 'Color', 'w');
tiledlayout(2, 2, 'TileSpacing', 'compact');
panels = ["length_error", "rate_error", "length_sigma"];
panel_titles = ["Length error L - L_{ref}", "Length-rate error", ...
                "Length 1\sigma from P"];
panel_labels = ["[m]", "[m/s]", "[m]"];
panel_scales = [length_scale_m, velocity_scale, length_scale_m];
for panel = 1:3
    nexttile;
    hold('on');
    for k = 1:numel(tether_runs)
        values = abs(tether_runs(k).history.(panels(panel))) * panel_scales(panel);
        semilogy(t_hours, max(values, 1e-12), 'Color', run_colors(k, :), ...
                 'LineStyle', run_styles(k), 'LineWidth', 1.5, ...
                 'DisplayName', tether_runs(k).name);
    end
    set(gca, 'YScale', 'log');
    grid('on');
    xlabel('Time [h]');
    ylabel("|" + panel_titles(panel) + "| " + panel_labels(panel));
    title(panel_titles(panel));
    if panel == 3
        floor_sigma = tether_runs(1).history.length_sigma_floor * length_scale_m;
        floor_line = plot(t_hours, floor_sigma, 'k:', 'LineWidth', 1.5, ...
                          'DisplayName', 'P integration-error floor');
        legend(floor_line, 'Location', 'southeast', 'FontSize', 7);
    end
    if panel == 1
        run_legend = legend('FontSize', 8, 'NumColumns', 4);
        run_legend.Layout.Tile = 'south';
    end
end
nexttile;
hold('on');
for body = 1:3
    rows = 3 * (body - 1) + (1:3);
    plot(t_hours, vecnorm(state_difference(rows, :)) * L_star, ...
         'Color', body_colors(body, :), 'LineWidth', line_widths(body), ...
         'DisplayName', body_names(body));
end
set(gca, 'YScale', 'log');
grid('on');
xlabel('Time [h]');
ylabel('Position difference [km]');
title('Check 2: |p_{SD} - p_{rigid}|, on-constraint start');
legend('Location', 'northwest', 'FontSize', 7);
sgtitle(['Check 1 (panels 1-3, process noise on): spring-damper (color) vs ' ...
         'rigid (gray); check 2 (panel 4)']);

%% Figure 6: linear covariance vs Monte Carlo (spring-damper, no process noise)
% The sampled states are propagated state-only through the nonlinear
% spring-damper dynamics. Agreement with the linearized P supports the
% spring-damper Jacobian over this span.
monte_carlo_count = 100;
monte_carlo_times = linspace(0, t_end, 25);
linear_args = {mass_2_fraction, tether_length_rate, accel_fn, ...
               no_process_noise, spring_damper};
linear_packed = propagatePacked(initial_state, initial_covariance, t_end, ...
                                monte_carlo_times, linear_args, ode_opts);

rng(1);
samples = zeros(monte_carlo_count, 18, numel(monte_carlo_times));
sample_args = {mass_2_fraction, tether_length_rate, accel_fn, struct(), ...
               spring_damper};
fprintf('\nMonte Carlo check (%d samples)...\n', monte_carlo_count);
for n = 1:monte_carlo_count
    sample_state = initial_state + initial_factor * randn(18, 1);
    sample_rhs = @(t, state) propagate(t, state, sample_args{:});
    sample_sol = ode113(sample_rhs, [0, t_end], sample_state, ...
                        odeset('RelTol', 1e-10, 'AbsTol', 1e-12));
    samples(n, :, :) = deval(sample_sol, monte_carlo_times);
end

% Spin-phase check: the linear phase 1-sigma (new P0 and the old diagonal P0)
% against the Monte Carlo spread of the wrapped phase error.
diagonal_packed = propagatePacked(initial_state, diagonal_covariance, t_end, ...
                                  monte_carlo_times, linear_args, ode_opts);
phase_new = phaseSigma(linear_packed);
phase_diagonal = phaseSigma(diagonal_packed);
separation_samples = samples(:, 4:6, :) - samples(:, 1:3, :);
nominal_phase = atan2(linear_packed(5, :) - linear_packed(2, :), ...
                      linear_packed(4, :) - linear_packed(1, :));
phase_error = angle(exp(1i * (squeeze(atan2(separation_samples(:, 2, :), ...
                                            separation_samples(:, 1, :))) - ...
                              nominal_phase)));
phase_sample = std(phase_error, 0, 1);
phase_resultant = abs(mean(exp(1i * phase_error), 1));
fprintf('\nSpin-phase 1-sigma [rad]: linear (new P0), linear (diagonal P0), MC\n');
for k = round(linspace(1, numel(monte_carlo_times), 7))
    fprintf('  t = %5.1f h   %8.4f   %8.4f   %8.4f   (R = %.3f)\n', ...
            monte_carlo_times(k) * T_star / 3600, phase_new(k), ...
            phase_diagonal(k), phase_sample(k), phase_resultant(k));
end

figure('Name', 'Linear covariance vs Monte Carlo', 'Color', 'w');
tiledlayout(1, 3, 'TileSpacing', 'compact');
monte_carlo_hours = monte_carlo_times * T_star / 3600;
quantity_names = ["position", "velocity"];
for panel = 1:2
    nexttile;
    hold('on');
    for body = 1:3
        if panel == 1
            rows = 3 * (body - 1) + (1:3);
            scale = L_star;
        else
            rows = 9 + 3 * (body - 1) + (1:3);
            scale = velocity_scale;
        end
        linear_sigma = blockSigma(linear_packed, rows) * scale;
        sample_sigma = zeros(1, numel(monte_carlo_times));
        for k = 1:numel(monte_carlo_times)
            sample_covariance = cov(samples(:, rows, k));
            sample_sigma(k) = sqrt(trace(sample_covariance)) * scale;
        end
        plot(monte_carlo_hours, linear_sigma, '-', ...
             'Color', body_colors(body, :), 'LineWidth', line_widths(body), ...
             'DisplayName', body_names(body) + " - linear P");
        plot(monte_carlo_hours, sample_sigma, 'o', ...
             'Color', body_colors(body, :), 'MarkerSize', 4 + 4 * (body == 1), ...
             'DisplayName', body_names(body) + " - Monte Carlo");
        fprintf('  %-20s %-8s final sigma ratio (MC / linear) = %.3f\n', ...
                body_names(body), quantity_names(panel), ...
                sample_sigma(end) / linear_sigma(end));
    end
    set(gca, 'YScale', 'log');
    grid('on');
    xlabel('Time [h]');
    if panel == 1
        ylabel('Position 1\sigma, sqrt(tr P_{pp}) [km]');
        title('Position uncertainty');
        legend('Location', 'northwest', 'FontSize', 7);
    else
        ylabel('Velocity 1\sigma, sqrt(tr P_{vv}) [m/s]');
        title('Velocity uncertainty');
    end
end
sgtitle(['Spring-damper tether (no process noise): linear covariance ' ...
         '(lines) vs Monte Carlo (markers)']);

% Third panel: length 1-sigma, the quantity the spring-damper is meant to
% hold down. The linear value comes from P, the sample value from |p2 - p1|.
nexttile;
hold('on');
dense_packed = propagatePacked(initial_state, initial_covariance, t_end, ...
                               t_eval, linear_args, ode_opts);
dense_history = lengthHistory(dense_packed, zeros(size(t_eval)), 0, ...
                              ode_opts.RelTol);
linear_history = lengthHistory(linear_packed, zeros(size(monte_carlo_times)), ...
                               0, ode_opts.RelTol);
sample_length = squeeze(vecnorm(samples(:, 4:6, :) - samples(:, 1:3, :), 2, 2));
plot(t_hours, dense_history.length_sigma * length_scale_m, '-', ...
     'LineWidth', 1.5, 'DisplayName', 'linear P');
plot(t_hours, dense_history.length_sigma_floor * length_scale_m, 'k:', ...
     'LineWidth', 1.5, 'DisplayName', 'P integration-error floor');
plot(monte_carlo_hours, std(sample_length, 0, 1) * length_scale_m, 'o', ...
     'MarkerSize', 4, 'DisplayName', 'Monte Carlo');
% Below the floor the linear value is integration error, not signal. The
% Monte Carlo spread (state-only runs at tighter tolerance) shows the true
% length spread has collapsed to ~1e-6 m.
fprintf('  Length (t = %.2f h)   sigma ratio (MC / linear) = %.3f\n', ...
        monte_carlo_hours(2), ...
        std(sample_length(:, 2)) / linear_history.length_sigma(2));
set(gca, 'YScale', 'log');
grid('on');
xlabel('Time [h]');
ylabel('Tether length 1\sigma [m]');
title('Tether length uncertainty');
legend('Location', 'northeast', 'FontSize', 7);

%% Figure 7: nonzero length rate, linear covariance vs Monte Carlo
% Repeats the Monte Carlo check with the tether reeling in and out. The
% reference length is L(t) = L_0 + Ldot * t, so the Jacobian and the
% covariance propagation are exercised with a nonzero length-rate command.
% The same initial covariance, payload, and random draws are used for each rate.
reel_rates_ms = [-0.5, 0.5];
reel_colors = [0.835, 0.369, 0.000; 0.000, 0.447, 0.698];
reel_count = 100;
rng(2);
reel_draws = randn(18, reel_count);
reel_times = linspace(0, t_end, 25);
reel_hours = reel_times * T_star / 3600;

figure('Name', 'Nonzero length rate', 'Color', 'w');
tiledlayout(1, 3, 'TileSpacing', 'compact');
reel_axes = [nexttile, nexttile, nexttile];
hold(reel_axes(1), 'on');
hold(reel_axes(2), 'on');
hold(reel_axes(3), 'on');
fprintf('\nNonzero length rate (%d samples each)\n', reel_count);
for r = 1:numel(reel_rates_ms)
    reel_rate_nd = reel_rates_ms(r) / 1000 * T_star / L_star;
    reel_terminals = coa.orbit.tetherInitialState(halo_ic, ...
                                                  tether_length_km / L_star, ...
                                                  mass_2_fraction, ...
                                                  spin_rate * T_star, ...
                                                  reel_rate_nd, lhat_0);
    reel_state = [reel_terminals(1:6); initial_state(7:9); ...
                  reel_terminals(7:12); initial_state(16:18)];
    reel_args = {mass_2_fraction, reel_rate_nd, accel_fn, no_process_noise, ...
                 spring_damper};
    reel_packed = propagatePacked(reel_state, initial_covariance, t_end, ...
                                  reel_times, reel_args, ode_opts);

    reel_samples = zeros(reel_count, 18, numel(reel_times));
    reel_sample_args = {mass_2_fraction, reel_rate_nd, accel_fn, struct(), ...
                        spring_damper};
    for n = 1:reel_count
        sample_rhs = @(t, state) propagate(t, state, reel_sample_args{:});
        sample_sol = ode113(sample_rhs, [0, t_end], ...
                            reel_state + initial_factor * reel_draws(:, n), ...
                            odeset('RelTol', 1e-10, 'AbsTol', 1e-12));
        reel_samples(n, :, :) = deval(sample_sol, reel_times);
    end

    reference_length = spring_damper.initial_length + reel_rate_nd * reel_times;
    reel_history = lengthHistory(reel_packed, reference_length, reel_rate_nd, ...
                                 ode_opts.RelTol);
    reel_name = sprintf('Ldot = %+.1f m/s', reel_rates_ms(r));

    for panel = 1:2
        if panel == 1
            rows = 1:3;
            scale = L_star;
        else
            rows = 10:12;
            scale = velocity_scale;
        end
        linear_sigma = blockSigma(reel_packed, rows) * scale;
        sample_sigma = zeros(1, numel(reel_times));
        for k = 1:numel(reel_times)
            sample_sigma(k) = sqrt(trace(cov(reel_samples(:, rows, k)))) * scale;
        end
        plot(reel_axes(panel), reel_hours, linear_sigma, '-', ...
             'Color', reel_colors(r, :), 'LineWidth', 1.5, ...
             'DisplayName', reel_name + " - linear P");
        plot(reel_axes(panel), reel_hours, sample_sigma, 'o', ...
             'Color', reel_colors(r, :), 'MarkerSize', 4, ...
             'DisplayName', reel_name + " - Monte Carlo");
        fprintf('  %s  body 1 %-8s final sigma ratio (MC / linear) = %.3f\n', ...
                reel_name, quantity_names(panel), ...
                sample_sigma(end) / linear_sigma(end));
    end
    rows = 7:9;
    payload_linear = blockSigma(reel_packed, rows);
    payload_sample = sqrt(trace(cov(reel_samples(:, rows, end))));
    fprintf('  %s  payload position final sigma ratio = %.3f\n', reel_name, ...
            payload_sample / payload_linear(end));

    plot(reel_axes(3), reel_hours, ...
         max(abs(reel_history.length_error) * length_scale_m, 1e-12), '-', ...
         'Color', reel_colors(r, :), 'LineWidth', 1.5, ...
         'DisplayName', reel_name + " - |L - L_{ref}| [m]");
    plot(reel_axes(3), reel_hours, ...
         max(abs(reel_history.rate_error) * velocity_scale, 1e-12), '--', ...
         'Color', reel_colors(r, :), 'LineWidth', 1.5, ...
         'DisplayName', reel_name + " - |dL/dt - Ldot| [m/s]");
    fprintf('  %s  final length %.3f km (reference %.3f km), error %.2g m\n', ...
            reel_name, norm(reel_packed(4:6, end) - reel_packed(1:3, end)) * ...
            L_star, reference_length(end) * L_star, ...
            reel_history.length_error(end) * length_scale_m);
end
for panel = 1:3
    set(reel_axes(panel), 'YScale', 'log');
    grid(reel_axes(panel), 'on');
    xlabel(reel_axes(panel), 'Time [h]');
end
legend(reel_axes(1), 'Location', 'northwest', 'FontSize', 7);
legend(reel_axes(2), 'Location', 'southeast', 'FontSize', 7);
ylim(reel_axes(3), [1e-12, 1e-2]);
legend(reel_axes(3), 'Location', 'northeast', 'FontSize', 7);
ylabel(reel_axes(1), 'Counterweight position 1\sigma [km]');
title(reel_axes(1), 'Position uncertainty');
ylabel(reel_axes(2), 'Counterweight velocity 1\sigma [m/s]');
title(reel_axes(2), 'Velocity uncertainty');
ylabel(reel_axes(3), 'Error [m] or [m/s]');
title(reel_axes(3), 'Constraint tracking (nominal state)');
sgtitle(['Reeling tether (no process noise): linear covariance (lines) vs ' ...
         'Monte Carlo (markers)']);

%% Figure 8: steady-state length covariance vs tension process noise
% Start from zero covariance and inject tension noise only. The settling span
% is 20 damping time constants, so the transient from P(0) = 0 is negligible.
% The full Cartesian covariance is projected onto the instantaneous length
% direction and compared with the stationary one-dimensional oscillator rule.
tension_psd_sweep = base_psd * 10.^(-2:2);
steady_duration = 20 / natural_frequency;
steady_covariance_rel_tol = 1e-12;
steady_covariance_opts = odeset( ...
                                'RelTol', steady_covariance_rel_tol, ...
                                'AbsTol', [1e-12 * ones(18, 1); 1e-28 * ones(324, 1)]);
steady_length_variance = zeros(size(tension_psd_sweep));
steady_length_variance_floor = zeros(size(tension_psd_sweep));
steady_length_variance_theory = tension_psd_sweep * gain^2 / ...
                                (2 * spring_damper.damping * ...
                                 spring_damper.stiffness);
for c = 1:numel(tension_psd_sweep)
    steady_noise = struct('tension_psd', tension_psd_sweep(c), ...
                          'transverse_psd', 0, ...
                          'body_3_acceleration_psd', 0);
    steady_rhs = @(t, packed_state) propagate( ...
                                              t, packed_state, ...
                                              mass_2_fraction, ...
                                              tether_length_rate, accel_fn, ...
                                              steady_noise, spring_damper);
    steady_sol = ode113(steady_rhs, [0, steady_duration], ...
                        [initial_state; zeros(324, 1)], ...
                        steady_covariance_opts);
    steady_state = steady_sol.y(1:18, end);
    steady_covariance = reshape(steady_sol.y(19:end, end), [18, 18]);
    separation = steady_state(4:6) - steady_state(1:3);
    direction = separation / norm(separation);
    length_gradient = [-direction; direction];
    position_covariance = steady_covariance(1:6, 1:6);
    steady_length_variance(c) = ...
        length_gradient' * position_covariance * length_gradient;
    steady_length_variance_floor(c) = ...
        steady_covariance_rel_tol * norm(position_covariance);
end

length_variance_scale_m2 = length_scale_m^2;
figure('Name', 'Steady-state length covariance', 'Color', 'w');
loglog(tension_psd_sweep, ...
       steady_length_variance * length_variance_scale_m2, 'o', ...
       'MarkerSize', 7, 'LineWidth', 1.5, ...
       'DisplayName', 'Full covariance propagation');
hold('on');
loglog(tension_psd_sweep, ...
       steady_length_variance_theory * length_variance_scale_m2, '--', ...
       'LineWidth', 1.5, ...
       'DisplayName', 'Stationary-oscillator theory');
loglog(tension_psd_sweep, ...
       steady_length_variance_floor * length_variance_scale_m2, ':', ...
       'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.5, ...
       'DisplayName', 'Estimated integration-error floor');
grid('on');
xlabel('Tension process-noise PSD q [nondimensional]');
ylabel('Steady-state tether length variance [m^2]');
title(sprintf('Length covariance after %.1f damping time constants', ...
              steady_duration * natural_frequency));
legend('Location', 'northwest');
for c = 1:numel(tension_psd_sweep)
    fprintf(['  q_T = %.1e: simulated/theory variance = %.3f, ', ...
             'variance/floor = %.2g\n'], ...
            tension_psd_sweep(c), ...
            steady_length_variance(c) / steady_length_variance_theory(c), ...
            steady_length_variance(c) / steady_length_variance_floor(c));
end

%% Local functions
function packed = propagatePacked(state_0, covariance_0, t_end, t_eval, args, ...
                                  ode_opts)
    % Propagate [state; P(:)] with the EKF propagator and sample at t_eval.
    % args are the propagator inputs after packed_state.
    propagator = @coa.orbit.tetheredThreeBodyEkfPropagation;
    rhs = @(t, packed_state) propagator(t, packed_state, args{:});
    sol = ode113(rhs, [0, t_end], [state_0; covariance_0(:)], ode_opts);
    packed = deval(sol, t_eval);
end

function history = lengthHistory(packed, reference_length, length_rate, rel_tol)
    % Length error, length-rate error, and the length 1-sigma implied by P.
    % rel_tol is the ODE RelTol used to propagate P; it sets the noise floor.
    n_times = size(packed, 2);
    history.length_error = zeros(1, n_times);
    history.rate_error = zeros(1, n_times);
    history.length_sigma = zeros(1, n_times);
    % The quadratic form g' * P * g subtracts large unconstrained directions
    % of P, so once the length variance has decayed below roughly
    % rel_tol * norm(P) it is integration error (and can be negative), not
    % signal.
    history.length_sigma_floor = zeros(1, n_times);
    for k = 1:n_times
        separation = packed(4:6, k) - packed(1:3, k);
        separation_rate = packed(13:15, k) - packed(10:12, k);
        length_value = norm(separation);
        direction = separation / length_value;
        history.length_error(k) = length_value - reference_length(k);
        history.rate_error(k) = dot(direction, separation_rate) - length_rate;
        covariance = reshape(packed(19:end, k), [18, 18]);
        gradient = [-direction; direction];
        position_covariance = covariance(1:6, 1:6);
        history.length_sigma(k) = ...
            sqrt(abs(gradient' * position_covariance * gradient));
        history.length_sigma_floor(k) = sqrt(rel_tol * norm(position_covariance));
    end
end

function sigma = phaseSigma(packed)
    % 1-sigma of the tether spin phase atan2(dy, dx), d = p2 - p1, from P.
    n_times = size(packed, 2);
    sigma = zeros(1, n_times);
    for k = 1:n_times
        separation = packed(4:6, k) - packed(1:3, k);
        covariance = reshape(packed(19:end, k), [18, 18]);
        gradient = zeros(18, 1);
        planar_gradient = [-separation(2); separation(1)] / ...
                          (separation(1)^2 + separation(2)^2);
        gradient(4:5) = planar_gradient;
        gradient(1:2) = -planar_gradient;
        sigma(k) = sqrt(gradient' * covariance * gradient);
    end
end

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
