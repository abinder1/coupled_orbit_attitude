% Copyright (c) 2026 Andrew Binder

function state = tetherInitialState(com_state, length_nd, mass_2_fraction, ...
                                    spin_rate_nd, length_dot, direction)
    % TETHERINITIALSTATE Initialize two tethered bodies about their center of mass.
    %   state = coa.orbit.tetherInitialState(com_state, length_nd, ...
    %       mass_2_fraction, spin_rate_nd, length_dot, direction) returns
    %   [p1; p2; v1; v2] in the nondimensional CR3BP rotating frame.
    %
    %   com_state is [position; velocity] in that frame. direction points
    %   from the center of mass toward body 1 and should be a unit 3-vector;
    %   it is used as supplied, without normalization. length_nd is the
    %   initial tether length, and length_dot is its nondimensional rate.
    %   mass_2_fraction is m2 / (m1 + m2). spin_rate_nd is the inertial spin
    %   rate about +z; the frame rotates about +z at nondimensional rate 1.
    arguments
        com_state (6, 1) double {mustBeFinite}
        length_nd (1, 1) double {mustBePositive, mustBeFinite}
        mass_2_fraction (1, 1) double {mustBeFinite, ...
                                       mustBeGreaterThanOrEqual(mass_2_fraction, 0), ...
                                       mustBeLessThanOrEqual(mass_2_fraction, 1)}
        spin_rate_nd (1, 1) double {mustBeFinite}
        length_dot (1, 1) double {mustBeFinite}
        direction (3, 1) double {mustBeFinite}
    end

    frame_rate = [0; 0; 1];
    relative_spin = spin_rate_nd * [0; 0; 1] - frame_rate;

    offset_1 = direction * (mass_2_fraction * length_nd);
    offset_2 = -direction * ((1 - mass_2_fraction) * length_nd);
    velocity_1 = cross(relative_spin, offset_1) + ...
                 mass_2_fraction * length_dot * direction;
    velocity_2 = cross(relative_spin, offset_2) - ...
                 (1 - mass_2_fraction) * length_dot * direction;

    state = [com_state(1:3) + offset_1; com_state(1:3) + offset_2; ...
             com_state(4:6) + velocity_1; com_state(4:6) + velocity_2];
end
