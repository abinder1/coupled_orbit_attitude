% Copyright (c) 2026 Andrew Binder

function [z, dz_drho, dz_drho_dot] = rangeRate(rho, rho_dot)
    % RANGERATE Range rate of a relative position and velocity.
    %   z = coa.meas.rangeRate(rho, rho_dot) returns
    %   z = rho . rho_dot / |rho| for the 3x1 position rho and velocity
    %   rho_dot. Implements Eq. eq:measurement-range-rate in
    %   derivations/ekf_equations.tex.
    %
    %   [z, dz_drho, dz_drho_dot] = coa.meas.rangeRate(...) also returns the
    %   1x3 Jacobians with respect to rho and rho_dot. The model is singular
    %   at rho = 0, where the error coa:meas:rangeRate:Singular is raised.
    % IMPL: eq:measurement-range-rate
    arguments
        rho (3, 1) double {mustBeFinite}
        rho_dot (3, 1) double {mustBeFinite}
    end

    distance = norm(rho);
    if distance == 0
        error('coa:meas:rangeRate:Singular', ...
              'Range rate is undefined at rho = 0.');
    end
    line_of_sight = rho / distance;

    z = dot(line_of_sight, rho_dot);
    if nargout < 2
        return
    end

    dz_drho = (rho_dot - z * line_of_sight)' / distance;
    dz_drho_dot = line_of_sight';
end
