% Copyright (c) 2026 Andrew Binder

function [z, dz_drho] = angles(rho)
    % ANGLES Right ascension and declination of a position vector.
    %   z = coa.meas.angles(rho) returns z = [alpha; delta] in radians for the
    %   3x1 position rho, with alpha = atan2(rho_y, rho_x) in [-pi, pi] and
    %   delta = asin(rho_z / |rho|) in [-pi/2, pi/2]. Implements h_1 and h_2 of
    %   Eq. eq:measurement-angles in derivations/ekf_equations.tex.
    %
    %   [z, dz_drho] = coa.meas.angles(rho) also returns the 2x3 Jacobian
    %   dz/drho. The model is singular on the z axis (rho_x = rho_y = 0),
    %   where the error coa:meas:angles:Singular is raised.
    % IMPL: eq:measurement-angles
    arguments
        rho (3, 1) double {mustBeFinite}
    end

    horizontal_sq = rho(1)^2 + rho(2)^2;
    if horizontal_sq == 0
        error('coa:meas:angles:Singular', ...
              'Angles are undefined for a position on the z axis.');
    end
    horizontal = sqrt(horizontal_sq);

    % atan2(rho_z, horizontal) equals asin(rho_z / |rho|) but is better
    % conditioned near the poles.
    z = [atan2(rho(2), rho(1)); atan2(rho(3), horizontal)];

    if nargout < 2
        return
    end

    rho_sq = horizontal_sq + rho(3)^2;
    dz_drho = [-rho(2) / horizontal_sq, ...
               rho(1) / horizontal_sq, ...
               0; ...
               -rho(1) * rho(3) / (rho_sq * horizontal), ...
               -rho(2) * rho(3) / (rho_sq * horizontal), ...
               horizontal / rho_sq];
end
