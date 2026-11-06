% Copyright (c) 2026 Andrew Binder

function [groups, angle_rows] = measurementGroups(mask)
    % MEASUREMENTGROUPS Metadata for the selected measurement-model groups.
    %   [groups, angle_rows] = coa.meas.measurementGroups(mask) returns the
    %   selected groups in ascending mask-bit order. Each struct has fields
    %   name (string), index (one-based group number), and rows (row vector
    %   of indices into the stacked measurement). angle_rows is a column
    %   vector of angular rows to wrap in coa.filter.ekfUpdate.
    %   A zero mask returns an empty struct array and a 0x1 angle vector.
    arguments
        mask (1, 1) double {mustBeInteger, mustBeNonnegative}
    end

    names = ["Ground optical angles", "Camera / star-tracker angles", ...
             "Camera angular size", "DSN range", "DSN range rate", ...
             "Crosslink range", "Crosslink range rate", "GNSS position"];
    sizes = [6, 4, 1, 3, 3, 2, 2, 9];
    angle_counts = [6, 4, 0, 0, 0, 0, 0, 0];
    if mask >= 2^numel(names)
        error('coa:meas:measurementGroups:UnknownGroup', ...
              'mask selects a group beyond the %d defined groups.', numel(names));
    end

    groups = struct('name', {}, 'index', {}, 'rows', {});
    angle_rows = zeros(0, 1);
    offset = 0;
    for group = 1:numel(names)
        if bitand(mask, 2^(group - 1)) > 0
            groups(end + 1) = struct('name', names(group), 'index', group, ...
                                     'rows', offset + (1:sizes(group))); %#ok<AGROW>
            angle_rows = [angle_rows; ...
                          (offset + (1:angle_counts(group)))']; %#ok<AGROW>
            offset = offset + sizes(group);
        end
    end
end
