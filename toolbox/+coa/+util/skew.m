% Copyright (c) 2026 Andrew Binder

function s = skew(v)
    % SKEW Skew-symmetric cross-product matrix.
    %   S = coa.util.skew(v) returns the 3x3 matrix S such that S*w = cross(v, w).
    %   Implements eq:skew (derivations/nomenclature.tex).
    arguments
        v (3, 1) double
    end
    s = [0   -v(3)  v(2)
         v(3)    0   -v(1)
         -v(2)  v(1)    0];
end
