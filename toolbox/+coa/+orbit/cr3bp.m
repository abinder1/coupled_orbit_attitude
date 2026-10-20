% Copyright (c) 2026 Andrew Binder

function [acceleration, da_dp, da_dv] = cr3bp(state, mu)
    % CR3BP Evaluate the CR3BP acceleration in normalized rotating-frame coordinates.
    %   a = coa.orbit.cr3bp(state, mu) returns the 3x1 acceleration for
    %   state = [x; y; z; vx; vy; vz], where mu is the secondary's mass
    %   fraction. The primaries are at x = -mu and x = 1 - mu, with unit
    %   separation and angular rate.
    %   [a, da_dp, da_dv] = coa.orbit.cr3bp(state, mu) also returns the 3x3
    %   partial derivatives of the acceleration with respect to position
    %   and velocity.
    arguments
        state (6, 1) double
        mu (1, 1) double
    end

    % Extract position and velocity
    r = state(1:3);
    v = state(4:6);

    % Compute distances to primary and secondary bodies
    re = r + [mu; 0; 0];
    rm = r - [1 - mu; 0; 0];

    ren = norm(re);
    rmn = norm(rm);

    reh = re / ren;
    rmh = rm / rmn;

    Wx = coa.util.skew([0; 0; 1]);
    Wx2 = Wx * Wx;

    % Compute acceleration
    acceleration = -(1 - mu) * reh / ren^2 - mu * rmh / rmn^2 - Wx2 * r - 2 * Wx * v;

    % Derivative outputs are only computed when requested.
    if nargout > 1
        I3 = eye(3);
        TE = (I3 - 3 * (reh * reh')) / ren^3;
        TM = (I3 - 3 * (rmh * rmh')) / rmn^3;

        da_dp = -(1 - mu) * TE - mu * TM - Wx2;
        da_dv = -2 * Wx;
    end
end
