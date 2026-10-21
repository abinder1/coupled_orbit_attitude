% Copyright (c) 2026 Andrew Binder

function [a, da_dp, da_dv] = tetherAccel(x1, x2, eta, length_dot, accel_fn)
    % TETHERACCEL Accelerations of two masses held at a prescribed separation.
    %   a = coa.orbit.tetherAccel(x1, x2, eta, length_dot, accel_fn) returns
    %   the stacked 6x1 constrained accelerations [a1; a2] of two masses with
    %   states x1 = [p1; v1] and x2 = [p2; v2]. eta = m2 / (m1 + m2), and
    %   length_dot is the rate of change of |p2 - p1|. accel_fn evaluates the
    %   unconstrained acceleration of a single mass and must follow the
    %   coa.orbit.cr3bp convention, [a, da_dp, da_dv] = accel_fn(x), for
    %   example @(x) coa.orbit.cr3bp(x, mu). Implements Secs. 2.4 and 3 of
    %   derivations/lagrangian_based.tex.
    %
    %   [a, da_dp, da_dv] = coa.orbit.tetherAccel(...) also returns the 6x6
    %   partials of [a1; a2] with respect to [p1; p2] and [v1; v2].
    arguments
        x1 (6, 1) double
        x2 (6, 1) double
        eta (1, 1) double
        length_dot (1, 1) double
        accel_fn (1, 1) function_handle
    end

    p1 = x1(1:3);
    v1 = x1(4:6);
    p2 = x2(1:3);
    v2 = x2(4:6);

    if nargout > 1
        [a1, a1_dp, a1_dv] = accel_fn(x1);
        [a2, a2_dp, a2_dv] = accel_fn(x2);
    else
        a1 = accel_fn(x1);
        a2 = accel_fn(x2);
    end

    d = p2 - p1;
    d_dot = v2 - v1;
    dd = dot(d, d);
    kappa = dot(d_dot, d_dot) - length_dot^2;

    % Correction shared by both masses, D * (a2 - a1) + q, applied with
    % weights eta and -(1 - eta). Since D * x = d * (d . x) / (d . d) and
    % q is parallel to d, no 3x3 matrices are needed for the accelerations.
    c = d * ((dot(d, a2 - a1) + kappa) / dd);

    a = [a1 + eta * c; a2 - (1 - eta) * c];

    if nargout > 1
        I3 = eye(3);
        D = (d * d') / dd;
        I_minus_2D = I3 - 2 * D;

        % Terms from differentiating D * (a2 - a1) and q (Sec. 3.4)
        S = (d * (a2 - a1)' + dot(d, a2 - a1) * I_minus_2D) / dd;
        K = kappa / dd * I_minus_2D;
        W = 2 / dd * (d * d_dot');
        SK = S + K;

        eta_D = eta * D;
        comp_D = (1 - eta) * D;
        comp_SK = (1 - eta) * SK;
        comp_W = (1 - eta) * W;

        da_dp = [(I3 - eta_D) * a1_dp - eta * SK, eta_D * a2_dp + eta * SK; ...
                 comp_D * a1_dp + comp_SK, (I3 - comp_D) * a2_dp - comp_SK];

        da_dv = [(I3 - eta_D) * a1_dv - eta * W, eta_D * a2_dv + eta * W; ...
                 comp_D * a1_dv + comp_W, (I3 - comp_D) * a2_dv - comp_W];
    end
end
