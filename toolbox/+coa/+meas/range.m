% Copyright (c) 2026 Andrew Binder

function [z, dz_drho] = range(rho)
    % RANGE Range (Euclidean norm) of a position vector.
    %   z = coa.meas.range(rho) returns z = |rho| for the 3x1 position rho.
    %   Implements Eq. eq:measurement-range in derivations/ekf_equations.tex.
    %
    %   [z, dz_drho] = coa.meas.range(rho) also returns the 1x3 Jacobian
    %   dz/drho = rho' / |rho|. The range value is defined at rho = 0, but its
    %   Jacobian is singular there and raises coa:meas:range:Singular.
    % IMPL: eq:measurement-range
    arguments
        rho (3, 1) double {mustBeFinite}
    end

    z = norm(rho);
    if z == 0
        if nargout > 1
            error('coa:meas:range:Singular', ...
                  'Range Jacobian is undefined at rho = 0.');
        end
        return
    end
    if nargout < 2
        return
    end
    dz_drho = rho' / z;
end
