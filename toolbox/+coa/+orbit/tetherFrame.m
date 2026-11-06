% Copyright (c) 2026 Andrew Binder

function frame = tetherFrame(relative_position, relative_velocity)
    % TETHERFRAME Construct a tether-aligned orthonormal frame.
    %   frame = coa.orbit.tetherFrame(relative_position, relative_velocity)
    %   returns a 3x3 matrix whose rows are frame axes in the input coordinates.
    %   The first axis points along relative_position, the third along
    %   cross(relative_position, relative_velocity), and the second completes
    %   the right-handed triad. Inputs are 3x1 vectors in the same frame and
    %   compatible units; relative_velocity is the velocity in that frame.
    %   This function returns orientation only: it does not transform
    %   velocities or covariances or add rotating-frame transport terms.
    %
    %   Zero separation and a zero cross product raise
    %   coa:orbit:tetherFrame:Singular.
    arguments
        relative_position (3, 1) double {mustBeFinite}
        relative_velocity (3, 1) double {mustBeFinite}
    end

    normal = cross(relative_position, relative_velocity);
    if ~any(relative_position) || ~any(normal)
        error('coa:orbit:tetherFrame:Singular', ...
              'Tether separation and its cross product with velocity must be nonzero.');
    end
    e1 = relative_position / norm(relative_position);
    e3 = normal / norm(normal);
    frame = [e1, cross(e3, e1), e3]';
end
