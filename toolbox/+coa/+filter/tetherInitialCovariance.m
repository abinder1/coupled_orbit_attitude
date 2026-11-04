% Copyright (c) 2026 Andrew Binder

function [covariance, factor] = tetherInitialCovariance(mass_2_fraction, sigma)
    % TETHERINITIALCOVARIANCE Covariance with a loose COM and tight relative state.
    %   [P, A] = coa.filter.tetherInitialCovariance(eta, sigma) returns the
    %   18x18 covariance P = A * A' of the state x = [p1; p2; p3; v1; v2; v3]
    %   and its factor A, for a tether pair (bodies 1 and 2) whose center of
    %   mass (COM) is poorly known but whose bodies are known well relative to
    %   the COM. Body 3 (the payload) is independent of the pair.
    %
    %   The pair is described by independent, isotropic, zero-mean errors in
    %       c = (1 - eta) * p1 + eta * p2     (COM position),
    %       r = p2 - p1                       (relative position),
    %   and the same for the velocities, with eta = m2 / (m1 + m2). The bodies
    %   follow as p1 = c - eta * r and p2 = c + (1 - eta) * r. Because the
    %   map is linear, P is exact, and A * randn(18, 1) draws errors that
    %   have exactly covariance P.
    %
    %   A tight relative sigma keeps errors almost parallel to the tether
    %   constraint manifold: a COM error moves both bodies together, so it
    %   leaves the separation length and rate unchanged.
    %
    %   Inputs
    %     mass_2_fraction  eta = m2 / (m1 + m2), scalar in (0, 1).
    %     sigma            Struct of nonnegative scalar standard deviations in
    %                      the units of the state, with fields
    %                        com_position, com_velocity,
    %                        relative_position, relative_velocity,
    %                        payload_position, payload_velocity.
    %
    %   Outputs
    %     covariance       18x18 symmetric covariance P.
    %     factor           18x18 matrix A with P = A * A'; the columns are
    %                      ordered [c; r; p3; vc; vr; v3] (positions then
    %                      velocities).
    arguments
        mass_2_fraction (1, 1) double {mustBeFinite}
        sigma (1, 1) struct
    end

    if ~(mass_2_fraction > 0 && mass_2_fraction < 1)
        error("coa:filter:tetherInitialCovariance:InvalidMassFraction", ...
              "mass_2_fraction must lie strictly between 0 and 1.");
    end

    names = ["com_position", "com_velocity", "relative_position", ...
             "relative_velocity", "payload_position", "payload_velocity"];
    for name = names
        if ~isfield(sigma, name)
            error("coa:filter:tetherInitialCovariance:MissingField", ...
                  "sigma is missing the field '%s'.", name);
        end
        value = sigma.(name);
        if ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
             value >= 0)
            error("coa:filter:tetherInitialCovariance:InvalidSigma", ...
                  "sigma.%s must be a finite nonnegative scalar.", name);
        end
    end

    eye_3 = eye(3);
    zero_3 = zeros(3);
    % Rows are [p1; p2; p3]; columns are [c; r; p3].
    to_bodies = [eye_3, -mass_2_fraction * eye_3, zero_3; ...
                 eye_3, (1 - mass_2_fraction) * eye_3, zero_3; ...
                 zero_3, zero_3, eye_3];

    position_scale = diag(kron([sigma.com_position, ...
                                sigma.relative_position, ...
                                sigma.payload_position], ones(1, 3)));
    velocity_scale = diag(kron([sigma.com_velocity, ...
                                sigma.relative_velocity, ...
                                sigma.payload_velocity], ones(1, 3)));

    factor = blkdiag(to_bodies * position_scale, to_bodies * velocity_scale);
    covariance = factor * factor';
end
