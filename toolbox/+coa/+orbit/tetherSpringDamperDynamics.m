% Copyright (c) 2026 Andrew Binder

function state_derivative = tetherSpringDamperDynamics(time, state, ...
                                                       mass_2_fraction, ...
                                                       rest_length, length_rate, ...
                                                       stiffness, damping, accel_fn, ...
                                                       use_feedforward)
    % TETHERSPRINGDAMPERDYNAMICS State-only dynamics of a spring-damper tether.
    %   xdot = coa.orbit.tetherSpringDamperDynamics(t, x, eta, L, L_dot, k, c, ...
    %   accel_fn) returns the derivative of the 12x1 state x = [p1; p2; v1; v2]
    %   of two masses joined by a spring and a damper that both act along
    %   d = p2 - p1, instead of the rigid constraint used by tetherAccel.
    %
    %   With u = d / |d|, the scalar force on body 1 along +u is
    %       F = k * (|d| - L) + c * (u . (v2 - v1) - L_dot) + F_ff,
    %   with L = L(t) the unstretched length. The spring and damper are
    %   linear, always active (they also push when |d| < L), and apply equal
    %   and opposite forces, so the center of mass feels no tether force.
    %
    %   The optional feedforward term F_ff (enabled with use_feedforward) is
    %   the tension the rigid constraint of tetherAccel would apply at the
    %   current state, eta * (1 - eta) * (d . (a2 - a1) + |d_dot|^2 - L_dot^2)
    %   / |d|. It cancels the centripetal and differential-gravity loads that
    %   vary with the current length and spin, so the spring only corrects
    %   deviations. With it, the equilibrium length is the rest length itself
    %   (no static stretch) and, on the constraint manifold, the dynamics
    %   match the rigid model.
    %
    %   Inputs
    %     t                Time, so that L(t) = rest_length + length_rate * t.
    %     x                [p1; p2; v1; v2], 12x1, in the frame of accel_fn.
    %     mass_2_fraction  eta = m2 / (m1 + m2); body 1 has fraction 1 - eta.
    %     rest_length      Unstretched length L(0).
    %     length_rate      Prescribed rate of change of L (0 for a fixed tether).
    %     stiffness        k, force per unit total mass per unit length, i.e.
    %                      F_physical / (m1 + m2) = k * stretch.
    %     damping          c, force per unit total mass per unit speed.
    %     accel_fn         Single-body acceleration, a = accel_fn([p; v]), as in
    %                      coa.orbit.cr3bp, for example @(x) cr3bp(x, mu).
    %     use_feedforward  Optional logical, default false; adds F_ff.
    %
    %   Units follow accel_fn (nondimensional for CR3BP). The undamped
    %   stretch frequency is sqrt(k / (eta * (1 - eta))) and the damping ratio is
    %   c / (2 * sqrt(k * eta * (1 - eta))) for the reduced-mass pair; both
    %   must be resolved by the integrator. A stiff tether (high frequency
    %   relative to the orbit) makes the system stiff, and ode15s may be needed.
    arguments
        time (1, 1) double
        state (12, 1) double
        mass_2_fraction (1, 1) double
        rest_length (1, 1) double
        length_rate (1, 1) double
        stiffness (1, 1) double
        damping (1, 1) double
        accel_fn (1, 1) function_handle
        use_feedforward (1, 1) logical = false
    end

    position_1 = state(1:3);
    position_2 = state(4:6);
    velocity_1 = state(7:9);
    velocity_2 = state(10:12);

    acceleration_1 = accel_fn([position_1; velocity_1]);
    acceleration_2 = accel_fn([position_2; velocity_2]);

    separation = position_2 - position_1;
    distance = norm(separation);
    direction = separation / distance;
    relative_velocity = velocity_2 - velocity_1;
    stretch_rate = dot(direction, relative_velocity);

    force = stiffness * (distance - (rest_length + length_rate * time)) + ...
            damping * (stretch_rate - length_rate);

    if use_feedforward
        kappa = dot(relative_velocity, relative_velocity) - length_rate^2;
        load = dot(separation, acceleration_2 - acceleration_1) + kappa;
        force = force + mass_2_fraction * (1 - mass_2_fraction) * ...
                load / distance;
    end

    % Force is per unit total mass; body 1 is pulled toward body 2.
    acceleration_1 = acceleration_1 + force * direction / (1 - mass_2_fraction);
    acceleration_2 = acceleration_2 - force * direction / mass_2_fraction;

    state_derivative = [velocity_1; velocity_2; acceleration_1; acceleration_2];
end
