% Copyright (c) 2026 Andrew Binder

function [z, dz_drho] = pnt(rho)
    % PNT Direct position measurement (lunar PN&T).
    %   z = coa.meas.pnt(rho) returns the 3x1 measurement z = rho for the 3x1
    %   position rho. Implements the lunar PN&T rows (17-19) of the sensor
    %   table in Sec. C of derivations/ekf_equations.tex. Only position is
    %   measured; no velocity measurement is available. The output uses the
    %   same coordinate frame and units as rho.
    %
    %   [z, dz_drho] = coa.meas.pnt(rho) also returns the 3x3 Jacobian
    %   dz/drho = I. The model is linear and has no singularities.
    % IMPL: eq:measurement-pnt
    arguments
        rho (3, 1) double {mustBeFinite}
    end

    z = rho;
    if nargout < 2
        return
    end
    dz_drho = eye(3);
end
