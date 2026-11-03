% Copyright (c) 2026 Andrew Binder

function [a, tension] = tetherSpringDamperAccel(t, x1, x2, eta, initial_length, ...
                                                length_dot, stiffness, damping, ...
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
    %   system becomes very stiff. This is a states-only function; it does not
    %   return Jacobians.
    %
    %   [a, tension] = coa.orbit.tetherSpringDamperAccel(...) also returns the
    %   scalar tension, the relative acceleration (force divided by reduced
    %   mass) the tether applies along the line of sight. It is positive when
    %   the tether pulls the masses together.
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

    a1 = accel_fn(x1);
    a2 = accel_fn(x2);

    % L_ddot = e . d_ddot + centripetal, with d_ddot = (a2 - a1) - tension * e.
    tension = dot(line_of_sight, a2 - a1) + centripetal - length_accel;
    correction = tension * line_of_sight;

    a = [a1 + eta * correction; a2 - (1 - eta) * correction];
end
