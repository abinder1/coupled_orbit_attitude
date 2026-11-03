% Copyright (c) 2026 Andrew Binder

function packed_derivative = tetheredThreeBodyEkfPropagation(normalized_time, ...
                                                             packed_state, ...
                                                             mass_2_fraction, ...
                                                             tether_length_rate, ...
                                                             accel_fn, ...
                                                             process_noise, ...
                                                             spring_damper)
    % TETHEREDTHREEBODYEKFPROPAGATION ODE right-hand side for EKF propagation.
    %   packed_derivative = coa.orbit.tetheredThreeBodyEkfPropagation( ...
    %       normalized_time, packed_state, mass_2_fraction, ...
    %       tether_length_rate, accel_fn)
    %   propagates the three-body state only. packed_state is the 18x1 state
    %   [p1; p2; p3; v1; v2; v3] (positions, then velocities; each block 3x1)
    %   and packed_derivative is its 18x1 time derivative.
    %
    %   packed_derivative = coa.orbit.tetheredThreeBodyEkfPropagation( ...
    %       normalized_time, packed_state, mass_2_fraction, ...
    %       tether_length_rate, accel_fn, process_noise)
    %   also propagates the 18x18 state covariance P. packed_state is then
    %   342x1, [state; P(:)] (P column-stacked, same state ordering), and
    %   packed_derivative is [state_dot; P_dot(:)]. This implements the
    %   Lyapunov equation (Eq. 5) of derivations/ekf_equations.tex:
    %       P_dot = F*P + P*F' + G*Q*G'.
    %
    %   packed_derivative = coa.orbit.tetheredThreeBodyEkfPropagation( ...
    %       normalized_time, packed_state, mass_2_fraction, ...
    %       tether_length_rate, accel_fn, process_noise, spring_damper)
    %   replaces the rigid tether constraint with the spring-damper tether of
    %   coa.orbit.tetherSpringDamperAccel, which drives |p2 - p1| toward
    %   L(t) = initial_length + tether_length_rate * normalized_time. Here
    %   normalized_time is the time since the instant at which the length equals
    %   initial_length, so it must share the origin of the truth model; a
    %   solver restarted at t = 0 needs initial_length updated to the current
    %   reference. spring_damper is a scalar struct with fields initial_length,
    %   stiffness, and damping (nondimensional, k = wn^2, c = 2 * zeta * wn).
    %   An empty struct (the default) keeps the rigid constraint. The process
    %   noise is unchanged and was derived for the rigid tether.
    %
    %   The mode is selected by the length of packed_state (18 or 342); any
    %   other length is rejected. A 342-element packed_state requires
    %   process_noise to contain all three fields below (otherwise MATLAB
    %   errors on the missing field). For an 18-element packed_state,
    %   process_noise is ignored and the acceleration Jacobians are not
    %   evaluated.
    %
    %   Bodies: 1 is the counterweight, 2 the catch mechanism, and 3 the
    %   external payload. Bodies 1 and 2 are coupled by coa.orbit.tetherAccel;
    %   body 3 moves under accel_fn alone, so F has no coupling between body 3
    %   and bodies 1-2.
    %
    %   Inputs:
    %     normalized_time     Nondimensional time. Unused for the rigid tether
    %                         (accel_fn and the process-noise PSDs are assumed
    %                         time-invariant); sets L(t) for the spring-damper.
    %     packed_state        18x1 or 342x1 column vector, described above.
    %     mass_2_fraction     eta = m2 / (m1 + m2); body 1 has mass fraction
    %                         1 - eta.
    %     tether_length_rate  Prescribed rate of change of |p2 - p1|, in
    %                         nondimensional units.
    %     accel_fn            Single-body acceleration following the
    %                         coa.orbit.cr3bp convention,
    %                         [a, da_dp, da_dv] = accel_fn(state), with state
    %                         6x1 [position; velocity]. The Jacobians are
    %                         requested only when propagating covariance.
    %     process_noise       Scalar struct of constant spectral densities
    %                         (extra fields are ignored):
    %                           tension_psd              q_T, radial tether-
    %                                                    tension noise (Eq. 20)
    %                           transverse_psd           q_perp, transverse
    %                                                    tether-force noise
    %                                                    (Eq. 26)
    %                           body_3_acceleration_psd  white acceleration
    %                                                    noise on body 3 (an
    %                                                    extension beyond
    %                                                    ekf_equations.tex)
    %
    %   Output:
    %     packed_derivative   Column vector the same size as packed_state.
    arguments
        normalized_time (1, 1) double
        packed_state double {mustBeStateVectorLength}
        mass_2_fraction (1, 1) double
        tether_length_rate (1, 1) double
        accel_fn (1, 1) function_handle
        process_noise (1, 1) struct = struct()
        spring_damper (1, 1) struct = struct()
    end

    propagate_covariance = numel(packed_state) > 18;
    use_spring_damper = ~isempty(fieldnames(spring_damper));
    if use_spring_damper
        spring_args = {spring_damper.initial_length, tether_length_rate, ...
                       spring_damper.stiffness, spring_damper.damping};
    end

    %% Unpack the state
    position_1 = packed_state(1:3);
    position_2 = packed_state(4:6);
    position_3 = packed_state(7:9);
    velocity_1 = packed_state(10:12);
    velocity_2 = packed_state(13:15);
    velocity_3 = packed_state(16:18);

    state_1 = [position_1; velocity_1];  % Counterweight
    state_2 = [position_2; velocity_2];  % Catch mechanism
    state_3 = [position_3; velocity_3];  % External payload

    %% State propagation: x_dot = f(x)
    % Acceleration Jacobians are only evaluated when covariance is propagated.
    if propagate_covariance
        if use_spring_damper
            [tether_acceleration, ...
             tether_acceleration_position_jacobian, ...
             tether_acceleration_velocity_jacobian] = ...
                coa.orbit.tetherSpringDamperAccel(normalized_time, state_1, ...
                                                  state_2, mass_2_fraction, ...
                                                  spring_args{:}, accel_fn);
        else
            [tether_acceleration, ...
             tether_acceleration_position_jacobian, ...
             tether_acceleration_velocity_jacobian] = ...
                coa.orbit.tetherAccel(state_1, state_2, mass_2_fraction, ...
                                      tether_length_rate, accel_fn);
        end
        [acceleration_3, ...
         body_3_acceleration_position_jacobian, ...
         body_3_acceleration_velocity_jacobian] = accel_fn(state_3);
    else
        if use_spring_damper
            tether_acceleration = ...
                coa.orbit.tetherSpringDamperAccel(normalized_time, state_1, ...
                                                  state_2, mass_2_fraction, ...
                                                  spring_args{:}, accel_fn);
        else
            tether_acceleration = coa.orbit.tetherAccel(state_1, state_2, ...
                                                        mass_2_fraction, ...
                                                        tether_length_rate, ...
                                                        accel_fn);
        end
        acceleration_3 = accel_fn(state_3);
    end

    state_derivative = [velocity_1; velocity_2; velocity_3; ...
                        tether_acceleration; acceleration_3];

    if ~propagate_covariance
        packed_derivative = state_derivative;
        return
    end

    %% Linearized dynamics: F = df/dx
    % With x = [positions; velocities], F = [0, I; da/dp, da/dv]. Bodies 1-2
    % are coupled through the tether; body 3 is decoupled from them.
    acceleration_position_jacobian = ...
        blkdiag(tether_acceleration_position_jacobian, ...
                body_3_acceleration_position_jacobian);
    acceleration_velocity_jacobian = ...
        blkdiag(tether_acceleration_velocity_jacobian, ...
                body_3_acceleration_velocity_jacobian);

    dynamics_jacobian = [zeros(9), eye(9); ...
                         acceleration_position_jacobian, ...
                         acceleration_velocity_jacobian];

    %% Process noise: G*Q*G'
    tension_psd = process_noise.tension_psd;
    transverse_psd = process_noise.transverse_psd;
    body_3_acceleration_psd = process_noise.body_3_acceleration_psd;

    separation = position_2 - position_1;
    tether_direction = separation / norm(separation);
    transverse_projector = eye(3) - tether_direction * tether_direction';

    mass_1_fraction = 1 - mass_2_fraction;

    % Radial tether-tension noise, Eq. 20. Acceleration noise enters only the
    % velocity-derivative rows: body 1 gets e/(1-eta), body 2 gets -e/eta.
    tension_shape = [zeros(9, 1); ...
                     tether_direction / mass_1_fraction; ...
                     -tether_direction / mass_2_fraction; ...
                     zeros(3, 1)];
    radial_noise_covariance = tension_psd * (tension_shape * tension_shape');

    % Transverse tether-force noise, Eq. 26, with P_perp = I - e*e'.
    transverse_shape = [zeros(9, 3); ...
                        transverse_projector / mass_1_fraction; ...
                        -transverse_projector / mass_2_fraction; ...
                        zeros(3, 3)];
    transverse_noise_covariance = ...
        transverse_psd * (transverse_shape * transverse_shape');

    % White acceleration noise on body 3 (extension, not in ekf_equations.tex).
    body_3_shape = [zeros(15, 3); eye(3)];
    body_3_noise_covariance = ...
        body_3_acceleration_psd * (body_3_shape * body_3_shape');

    process_noise_covariance = radial_noise_covariance + ...
                               transverse_noise_covariance + ...
                               body_3_noise_covariance;

    %% Covariance propagation: P_dot = F*P + P*F' + G*Q*G' (Eq. 5)
    covariance = reshape(packed_state(19:end), [18, 18]);
    covariance_derivative = dynamics_jacobian * covariance + ...
                            covariance * dynamics_jacobian' + ...
                            process_noise_covariance;

    packed_derivative = [state_derivative; covariance_derivative(:)];
end
