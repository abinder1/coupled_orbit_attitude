% Copyright (c) 2026 Andrew Binder

function [z, dz_drho_1, dz_drho_2] = angularSize(rho_1, rho_2)
    % ANGULARSIZE Angle between two position vectors.
    %   z = coa.meas.angularSize(rho_1, rho_2) returns the angle
    %   z = acos(rho_1 . rho_2 / (|rho_1| |rho_2|)) in [0, pi] in radians.
    %   Implements Eq. eq:measurement-angular-size in
    %   derivations/ekf_equations.tex.
    %
    %   [z, dz_drho_1, dz_drho_2] = coa.meas.angularSize(...) also returns the
    %   1x3 Jacobians with respect to rho_1 and rho_2. The angle is defined
    %   for parallel and antiparallel vectors, but its Jacobians are singular
    %   there (sin(z) < 1e-12). Either zero input is invalid.
    % IMPL: eq:measurement-angular-size
    arguments
        rho_1 (3, 1) double {mustBeFinite}
        rho_2 (3, 1) double {mustBeFinite}
    end

    norm_1 = norm(rho_1);
    norm_2 = norm(rho_2);
    if norm_1 == 0 || norm_2 == 0
        error('coa:meas:angularSize:Singular', ...
              'Angular size is undefined for a zero vector.');
    end
    direction_1 = rho_1 / norm_1;
    direction_2 = rho_2 / norm_2;

    % Clamp against roundoff pushing the cosine just outside [-1, 1].
    cosine = min(max(dot(direction_1, direction_2), -1), 1);
    z = acos(cosine);

    if nargout < 2
        return
    end

    sine = sin(z);
    if sine < 1e-12
        error('coa:meas:angularSize:Singular', ...
              ['Angular size Jacobian is undefined for parallel or ' ...
               'antiparallel vectors.']);
    end
    dz_drho_1 = -(direction_2 - cosine * direction_1)' / (sine * norm_1);
    dz_drho_2 = -(direction_1 - cosine * direction_2)' / (sine * norm_2);
end
