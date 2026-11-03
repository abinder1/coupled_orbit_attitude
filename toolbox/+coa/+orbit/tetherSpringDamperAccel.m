% Copyright (c) 2026 Andrew Binder

function [a, da_dp, da_dv, tension] = tetherSpringDamperAccel(t, x1, x2, eta, ...
                                                              initial_length, ...
                                                              length_dot, ...
                                                              stiffness, damping, ...
                                                              accel_fn)
    % TETHERSPRINGDAMPERACCEL Two masses driven toward a prescribed tether length.
    %   a = coa.orbit.tetherSpringDamperAccel(t, x1, x2, eta, initial_length, ...
    %   length_dot, stiffness, damping, accel_fn) returns the stacked 6x1
    %   accelerations [a1; a2] of two masses with states x1 = [p1; v1] and
    %   x2 = [p2; v2]. Instead of enforcing |p2 - p1| = L(t) exactly, as
    %   coa.orbit.tetherAccel does, a tension along the line of sight acts as a
    %   spring and damper on the length and length-rate errors relative to
    %
    %       L(t) = initial_length + length_dot * t.
    %
    %   The tension is chosen so that the length obeys
    %
    %       L_ddot = -stiffness * (L - L_ref) - damping * (L_dot - L_dot_ref),
    %
    %   which includes the centripetal term (|d_dot|^2 - L_dot^2) / L. Without
    %   that term the spring would have to carry the full spin load and the
    %   system becomes very stiff.
    %
    %   [a, da_dp, da_dv] = coa.orbit.tetherSpringDamperAccel(...) also returns
    %   the 6x6 partials of [a1; a2] with respect to [p1; p2] and [v1; v2], in
    %   the same layout as coa.orbit.tetherAccel. accel_fn must then return its
    %   own Jacobians, [a, da_dp, da_dv] = accel_fn(x). Unlike tetherAccel, the
    %   length rate is the measured |d|-dot rather than the constant length_dot,
    %   so the Jacobian is exact for the off-constraint states this model allows.
    %
    %   [a, da_dp, da_dv, tension] = coa.orbit.tetherSpringDamperAccel(...)
    %   also returns the scalar tension, the relative acceleration (force
    %   divided by reduced mass) the tether applies along the line of sight. It
    %   is positive when the tether pulls the masses together. Request fewer
    %   outputs to skip the Jacobian work.
    %
    %   Inputs:
    %     t               scalar time, in the same units as length_dot.
    %     x1, x2          6x1 states [position; velocity] of each mass.
    %     eta             m2 / (m1 + m2), in [0, 1].
    %     initial_length  L(0), positive.
    %     length_dot      constant rate of change of the reference length.
    %     stiffness       k >= 0, length-error gain in 1/time^2 (k = wn^2).
    %     damping         c >= 0, rate-error gain in 1/time (c = 2 * zeta * wn).
    %     accel_fn        function handle returning the unconstrained acceleration
    %                     of one mass, a = accel_fn(x), as for coa.orbit.cr3bp.
    %
    %   Notes:
    %     - The gains act on the relative coordinate, so they do not depend on
    %       the individual masses. The tension is split with weights eta and
    %       -(1 - eta) so the mass-weighted acceleration is unchanged.
    %     - The tension is not clipped at zero, so the tether can push as well
    %       as pull. With stiffness = damping = 0 and a state that already
    %       satisfies the length and rate, the result equals
    %       coa.orbit.tetherAccel.
    %     - The state must not have coincident masses (|p2 - p1| = 0).
    arguments
        t (1, 1) double {mustBeFinite}
        x1 (6, 1) double
        x2 (6, 1) double
        eta (1, 1) double {mustBeGreaterThanOrEqual(eta, 0), ...
                           mustBeLessThanOrEqual(eta, 1)}
        initial_length (1, 1) double {mustBePositive, mustBeFinite}
        length_dot (1, 1) double {mustBeFinite}
        stiffness (1, 1) double {mustBeNonnegative, mustBeFinite}
        damping (1, 1) double {mustBeNonnegative, mustBeFinite}
        accel_fn (1, 1) function_handle
    end

    d = x2(1:3) - x1(1:3);
    d_dot = x2(4:6) - x1(4:6);
    length = norm(d);
    line_of_sight = d / length;
    length_rate = dot(line_of_sight, d_dot);

    length_error = length - (initial_length + length_dot * t);
    rate_error = length_rate - length_dot;
    length_accel = -stiffness * length_error - damping * rate_error;
    centripetal = (dot(d_dot, d_dot) - length_rate^2) / length;

    if nargout > 1
        [a1, a1_dp, a1_dv] = accel_fn(x1);
        [a2, a2_dp, a2_dv] = accel_fn(x2);
    else
        a1 = accel_fn(x1);
        a2 = accel_fn(x2);
    end

    % L_ddot = e . d_ddot + centripetal, with d_ddot = (a2 - a1) - tension * e.
    relative_accel = a2 - a1;
    tension = dot(line_of_sight, relative_accel) + centripetal - length_accel;
    correction = tension * line_of_sight;

    a = [a1 + eta * correction; a2 - (1 - eta) * correction];

    if nargout > 1
        I3 = eye(3);
        transverse = I3 - line_of_sight * line_of_sight';
        rate_dd = (transverse * d_dot)' / length;

        % Gradient of the tension with its explicit dependence on d and d_dot;
        % the dependence through a2 - a1 is added below.
        tension_dd = (transverse * relative_accel)' / length + ...
                     (damping - 2 * length_rate / length) * rate_dd + ...
                     (stiffness - centripetal / length) * line_of_sight';
        tension_dv = 2 * (d_dot - length_rate * line_of_sight)' / length + ...
                     damping * line_of_sight';

        tension_dp1 = -tension_dd - line_of_sight' * a1_dp;
        tension_dp2 = tension_dd + line_of_sight' * a2_dp;
        tension_dv1 = -tension_dv - line_of_sight' * a1_dv;
        tension_dv2 = tension_dv + line_of_sight' * a2_dv;

        % correction = tension * e, and e depends on d only (de/dd = P / L).
        turning = tension * transverse / length;
        correction_dp1 = line_of_sight * tension_dp1 - turning;
        correction_dp2 = line_of_sight * tension_dp2 + turning;
        correction_dv1 = line_of_sight * tension_dv1;
        correction_dv2 = line_of_sight * tension_dv2;

        comp = 1 - eta;
        da_dp = [a1_dp + eta * correction_dp1, eta * correction_dp2; ...
                 -comp * correction_dp1, a2_dp - comp * correction_dp2];
        da_dv = [a1_dv + eta * correction_dv1, eta * correction_dv2; ...
                 -comp * correction_dv1, a2_dv - comp * correction_dv2];
    end
end
