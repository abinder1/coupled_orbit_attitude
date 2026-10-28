% Copyright (c) 2026 Andrew Binder

function mustBeStateVectorLength(packed_state)
    % Accept the 18-element state or the 342-element state and covariance,
    % as a true column (a size spec such as (:, 1) would reshape rows silently).
    if ~iscolumn(packed_state) || ~ismember(numel(packed_state), [18, 342])
        error('coa:orbit:tetheredThreeBodyEkfPropagation:InvalidPackedState', ...
              ['packed_state must be a column vector with 18 (state) or ', ...
               '342 (state and covariance) elements.']);
    end
end
