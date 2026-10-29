% Copyright (c) 2026 Andrew Binder

function [h, dh_dstate, noise_cov] = measurementModels(state, mask, meas_covs)
    % MEASUREMENTMODELS Stacked measurements and Jacobian for selected sensor groups.
    %   h = coa.meas.measurementModels(state, mask) evaluates every sensor
    %   group selected by the bit mask at the 18x1 EKF state
    %   [r1; r2; r3; v1; v2; v3] and returns the measurements stacked in group
    %   order (lowest bit first, bodies in order within a group).
    %
    %   [h, dh_dstate] = coa.meas.measurementModels(state, mask) also returns
    %   the stacked measurement Jacobian dh/dstate (numel(h)-by-18).
    %   The Jacobian is computed only when requested. A zero mask returns
    %   empty measurements, a 0x18 Jacobian, and (if requested) a 0x0
    %   covariance.
    %
    %   [h, dh_dstate, noise_cov] = coa.meas.measurementModels(state, mask, meas_covs)
    %   also returns the measurement noise covariance (numel(h)-by-numel(h)),
    %   block diagonal over the selected groups and rows. Sample noise outside
    %   the function, e.g. h + mvnrnd(zeros(1, numel(h)), noise_cov)'. meas_covs
    %   is a struct with one covariance per sensor, applied to each measurement
    %   that sensor makes:
    %       ground_angles       2x2  [alpha; delta] noise, rad^2 (group 1)
    %       camera_angles       2x2  [alpha; delta] noise, rad^2 (group 2, rows 1-2)
    %       star_tracker_angles 2x2  [alpha; delta] noise, rad^2 (group 2, rows 3-4)
    %       camera_angular_size 1x1  angular size noise, rad^2 (group 3)
    %       dsn_range           1x1  range noise, state length units^2 (group 4)
    %       dsn_range_rate      1x1  range-rate noise, state speed units^2 (group 5)
    %       crosslink_range     1x1  range noise, state length units^2 (group 6)
    %       crosslink_range_rate 1x1  range-rate noise, state speed units^2 (group 7)
    %       gnss_position       3x3  position noise, state length units^2 (group 8)
    %   Only the fields needed by the selected groups must be present.
    %   Requested covariances must be finite, real, symmetric, and positive
    %   definite.
    %   Covariances are checked only when noise_cov is requested.
    %
    %   Sensor groups (bit k of mask selects group k + 1, e.g. 2^0 for group 1):
    %       1  Ground optical angles of r1, r2, r3 (6 rows)
    %       2  Onboard camera angles of r2 - r3 and star-tracker angles of
    %          r2 - r1 (4 rows)
    %       3  Onboard camera angular size of the tether (1 row)
    %       4  DSN ground radar range of r1, r2, r3 (3 rows)
    %       5  DSN ground radar range rate of r1, r2, r3 (3 rows)
    %       6  RF crosslink range of r2 - r3 and r2 - r1 (2 rows)
    %       7  RF crosslink Doppler (range rate) of r2 - r3 and r2 - r1 (2 rows)
    %       8  Lunar PN&T/GNSS position of r1, r2, r3 (9 rows)
    %
    %   Every group is a local function with the signature
    %       [h, dh_dstate, noise_cov] = ...
    %           groupName(state, meas_covs, compute_jacobian)
    %   where noise_cov is the covariance of that group's rows of h. Groups are
    %   built from the primitives in this package. To add a group, write its
    %   function and append its handle to the list in groupFunctions.
    %
    %   Only sensor rows 1-19 in Sec. C of derivations/ekf_equations.tex are
    %   implemented. Accelerometers in rows 20-22 are candidates only; they
    %   are not measurement groups in this function.
    arguments
        state (18, 1) double {mustBeFinite}
        mask (1, 1) double {mustBeInteger, mustBeNonnegative}
        meas_covs (1, 1) struct = struct()
    end

    % Covariances are only needed, and so only checked, when noise_cov is requested.
    if nargout < 3
        meas_covs = [];
    end
    compute_jacobian = nargout > 1;

    group_functions = groupFunctions();

    if mask >= 2^numel(group_functions)
        error('coa:meas:measurementModels:UnknownGroup', ...
              'mask selects a group beyond the %d defined groups.', ...
              numel(group_functions));
    end

    h = zeros(0, 1);
    dh_dstate = zeros(0, 18);
    noise_cov = zeros(0, 0);

    for group = 1:numel(group_functions)
        if bitand(mask, 2^(group - 1)) > 0
            [group_h, group_dh_dstate, group_cov] = ...
                group_functions{group}(state, meas_covs, compute_jacobian);
            h = [h; group_h]; %#ok<AGROW>
            dh_dstate = [dh_dstate; group_dh_dstate]; %#ok<AGROW>
            noise_cov = blkdiag(noise_cov, group_cov);
        end
    end
end

function group_functions = groupFunctions()
    % Position in this list is the group number, i.e. the mask bit.
    group_functions = {@groundAngles, @spaceAngles, @tetherAngularSize, ...
                       @dsnRange, @dsnRangeRate, @crosslinkRange, ...
                       @crosslinkRangeRate, @gnssPosition};
end

%% Sensor groups

function [h, dh_dstate, noise_cov] = groundAngles(state, meas_covs, compute_jacobian)
    % Group 1: right ascension and declination of each body's position
    % (Eq. eq:measurement-angles, sensors 1-3).
    h = zeros(6, 1);
    dh_dstate = zeros(6, 18);
    for body = 1:3
        rows = 2 * body - 1:2 * body;
        selector = positionSelector(body);
        if compute_jacobian
            [h(rows), dh_drho] = coa.meas.angles(selector * state);
            dh_dstate(rows, :) = dh_drho * selector;
        else
            h(rows) = coa.meas.angles(selector * state);
        end
    end
    noise_cov = kron(eye(3), requireCovariance(meas_covs, 'ground_angles', 2));
end

function [h, dh_dstate, noise_cov] = spaceAngles(state, meas_covs, compute_jacobian)
    % Group 2: relative-geometry angles (Eq. eq:measurement-angles,
    % sensors 4-5).
    %   Rows 1-2: onboard camera angles of the catch mechanism seen from the
    %             payload, r2 - r3.
    %   Rows 3-4: star-tracker angles of the tether orientation, r2 - r1.
    h = zeros(4, 1);
    dh_dstate = zeros(4, 18);

    selector = positionSelector(2) - positionSelector(3);
    if compute_jacobian
        [h(1:2), dh_drho] = coa.meas.angles(selector * state);
        dh_dstate(1:2, :) = dh_drho * selector;
    else
        h(1:2) = coa.meas.angles(selector * state);
    end

    selector = positionSelector(2) - positionSelector(1);
    if compute_jacobian
        [h(3:4), dh_drho] = coa.meas.angles(selector * state);
        dh_dstate(3:4, :) = dh_drho * selector;
    else
        h(3:4) = coa.meas.angles(selector * state);
    end
    noise_cov = blkdiag(requireCovariance(meas_covs, 'camera_angles', 2), ...
                        requireCovariance(meas_covs, 'star_tracker_angles', 2));
end

function [h, dh_dstate, noise_cov] = ...
        tetherAngularSize(state, meas_covs, compute_jacobian)
    % Group 3: onboard camera angular size of the tether in the payload's sky
    % (Eq. eq:measurement-angular-size, sensor 6), the angle between the
    % lines of sight to the tether ends, r1 - r3 and r2 - r3.
    selector_1 = positionSelector(1) - positionSelector(3);
    selector_2 = positionSelector(2) - positionSelector(3);
    dh_dstate = zeros(1, 18);
    if compute_jacobian
        [h, dh_drho_1, dh_drho_2] = ...
            coa.meas.angularSize(selector_1 * state, selector_2 * state);
        dh_dstate = dh_drho_1 * selector_1 + dh_drho_2 * selector_2;
    else
        h = coa.meas.angularSize(selector_1 * state, selector_2 * state);
    end
    noise_cov = requireCovariance(meas_covs, 'camera_angular_size', 1);
end

function [h, dh_dstate, noise_cov] = dsnRange(state, meas_covs, compute_jacobian)
    % Group 4: DSN ground radar range to each body's position
    % (Eq. eq:measurement-range, sensors 7-9).
    h = zeros(3, 1);
    dh_dstate = zeros(3, 18);
    for body = 1:3
        selector = positionSelector(body);
        if compute_jacobian
            [h(body), dh_drho] = coa.meas.range(selector * state);
            dh_dstate(body, :) = dh_drho * selector;
        else
            h(body) = coa.meas.range(selector * state);
        end
    end
    noise_cov = kron(eye(3), requireCovariance(meas_covs, 'dsn_range', 1));
end

function [h, dh_dstate, noise_cov] = ...
        dsnRangeRate(state, meas_covs, compute_jacobian)
    % Group 5: DSN ground radar range rate of each body
    % (Eq. eq:measurement-range-rate, sensors 10-12).
    h = zeros(3, 1);
    dh_dstate = zeros(3, 18);
    for body = 1:3
        position_selector = positionSelector(body);
        velocity_selector = velocitySelector(body);
        position = position_selector * state;
        velocity = velocity_selector * state;
        if compute_jacobian
            [h(body), dh_drho, dh_drho_dot] = ...
                coa.meas.rangeRate(position, velocity);
            dh_dstate(body, :) = dh_drho * position_selector + ...
                                 dh_drho_dot * velocity_selector;
        else
            h(body) = coa.meas.rangeRate(position, velocity);
        end
    end
    noise_cov = kron(eye(3), requireCovariance(meas_covs, 'dsn_range_rate', 1));
end

function [h, dh_dstate, noise_cov] = ...
        crosslinkRange(state, meas_covs, compute_jacobian)
    % Group 6: RF crosslink range (Eq. eq:measurement-range, sensors 13 and 15).
    %   Row 1: payload to catch mechanism, |r2 - r3|.
    %   Row 2: catch mechanism to tether end, |r2 - r1|.
    h = zeros(2, 1);
    dh_dstate = zeros(2, 18);

    selector = positionSelector(2) - positionSelector(3);
    if compute_jacobian
        [h(1), dh_drho] = coa.meas.range(selector * state);
        dh_dstate(1, :) = dh_drho * selector;
    else
        h(1) = coa.meas.range(selector * state);
    end

    selector = positionSelector(2) - positionSelector(1);
    if compute_jacobian
        [h(2), dh_drho] = coa.meas.range(selector * state);
        dh_dstate(2, :) = dh_drho * selector;
    else
        h(2) = coa.meas.range(selector * state);
    end
    noise_cov = kron(eye(2), requireCovariance(meas_covs, 'crosslink_range', 1));
end

function [h, dh_dstate, noise_cov] = ...
        crosslinkRangeRate(state, meas_covs, compute_jacobian)
    % Group 7: RF crosslink Doppler, i.e. range rate
    % (Eq. eq:measurement-range-rate, sensors 14 and 16).
    %   Row 1: payload to catch mechanism, using r2 - r3 and v2 - v3.
    %   Row 2: catch mechanism to tether end, using r2 - r1 and v2 - v1.
    h = zeros(2, 1);
    dh_dstate = zeros(2, 18);

    position_selector = positionSelector(2) - positionSelector(3);
    velocity_selector = velocitySelector(2) - velocitySelector(3);
    position = position_selector * state;
    velocity = velocity_selector * state;
    if compute_jacobian
        [h(1), dh_drho, dh_drho_dot] = ...
            coa.meas.rangeRate(position, velocity);
        dh_dstate(1, :) = dh_drho * position_selector + ...
                          dh_drho_dot * velocity_selector;
    else
        h(1) = coa.meas.rangeRate(position, velocity);
    end

    position_selector = positionSelector(2) - positionSelector(1);
    velocity_selector = velocitySelector(2) - velocitySelector(1);
    position = position_selector * state;
    velocity = velocity_selector * state;
    if compute_jacobian
        [h(2), dh_drho, dh_drho_dot] = ...
            coa.meas.rangeRate(position, velocity);
        dh_dstate(2, :) = dh_drho * position_selector + ...
                          dh_drho_dot * velocity_selector;
    else
        h(2) = coa.meas.rangeRate(position, velocity);
    end
    noise_cov = kron(eye(2), requireCovariance(meas_covs, 'crosslink_range_rate', 1));
end

function [h, dh_dstate, noise_cov] = ...
        gnssPosition(state, meas_covs, compute_jacobian)
    % Group 8: lunar PN&T/GNSS position fix of each body (sensors 17-19).
    h = zeros(9, 1);
    dh_dstate = zeros(9, 18);
    for body = 1:3
        rows = 3 * body - 2:3 * body;
        selector = positionSelector(body);
        if compute_jacobian
            [h(rows), dh_drho] = coa.meas.pnt(selector * state);
            dh_dstate(rows, :) = dh_drho * selector;
        else
            h(rows) = coa.meas.pnt(selector * state);
        end
    end
    noise_cov = kron(eye(3), requireCovariance(meas_covs, 'gnss_position', 3));
end

%% State selectors
% A primitive's input is a linear function of the state, input = S * state,
% so the measurement Jacobian follows from the chain rule, dh/dstate =
% dh/dinput * S. Differences of bodies are differences of selectors, for
% example positionSelector(2) - positionSelector(3) for r2 - r3.

function selector = positionSelector(body)
    % 3x18 matrix picking the position of body 1, 2, or 3.
    identity = eye(18);
    selector = identity(3 * body - 2:3 * body, :);
end

function selector = velocitySelector(body)
    % 3x18 matrix picking the velocity of body 1, 2, or 3.
    identity = eye(18);
    selector = identity(9 + 3 * body - 2:9 + 3 * body, :);
end

%% Noise helpers

function covariance = requireCovariance(meas_covs, field, dimension)
    % Fetch a finite, real, symmetric positive-definite covariance.
    % An empty meas_covs means noise was not requested, so a placeholder is returned.
    if isempty(meas_covs)
        covariance = zeros(dimension);
        return
    end
    if ~isfield(meas_covs, field)
        error('coa:meas:measurementModels:MissingCovariance', ...
              'meas_covs.%s is required by the selected groups.', field);
    end
    covariance = meas_covs.(field);
    valid = isequal(size(covariance), [dimension, dimension]) && ...
            isfloat(covariance) && isreal(covariance) && ...
            all(isfinite(covariance), 'all');
    if valid
        if issymmetric(covariance)
            [~, chol_flag] = chol(covariance);
            valid = chol_flag == 0;
        else
            valid = false;
        end
    end
    if ~valid
        error('coa:meas:measurementModels:InvalidCovariance', ...
              ['meas_covs.%s must be a finite, real, symmetric positive ' ...
               'definite %dx%d matrix.'], ...
              field, dimension, dimension);
    end
end
