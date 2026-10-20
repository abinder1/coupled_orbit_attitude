% Copyright (c) 2026 Andrew Binder

function [state_dot, jacobian] = cr3bp(state, mu)
    % CR3BP Evaluate the state derivative in normalized rotating-frame coordinates.
    %   state is [x; y; z; vx; vy; vz], and mu is the secondary's mass
    %   fraction. The primaries are at x = -mu and x = 1 - mu, with unit
    %   separation and angular rate. state_dot is [v; acceleration], and
    %   jacobian, when requested, is d(state_dot)/d(state).
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

    % Compute acceleration and assemble the state derivative
    acceleration = -(1 - mu) * reh / ren^2 - mu * rmh / rmn^2 - Wx2 * r - 2 * Wx * v;

    state_dot = [v; acceleration];

    % Derivative-only callers do not need the position-derivative matrices.
    if nargout > 1
        I3 = eye(3);
        TE = (I3 - 3 * (reh * reh')) / ren^3;
        TM = (I3 - 3 * (rmh * rmh')) / rmn^3;

        jacobian = [zeros(3), I3; -(1 - mu) * TE - mu * TM - Wx2, -2 * Wx];
    end
end
