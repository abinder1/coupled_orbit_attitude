% Copyright (c) 2026 Andrew Binder

function sp = fitWindowedSpline(t, position, velocity, opts)
    % FITWINDOWEDSPLINE Windowed B-spline fit to position and velocity.
    %   sp = coa.interp.fitWindowedSpline(t, position, velocity) fits an
    %   N-by-D position and velocity history with a quintic B-spline, one
    %   independent spline per column. t is a strictly increasing N-element
    %   time vector; velocity is in position units per unit of t.
    %
    %   sp = coa.interp.fitWindowedSpline(t, position, velocity, opts)
    %   overrides any of the following fields in the optional scalar struct:
    %     order  B-spline order (default 6); interior continuity is C^(order-2).
    %     window_size      Window length (default 40, at least order).
    %     stride           Window stride (default 20, at least 1 and less
    %                      than window_size).
    %     velocity_weight  Positive relative velocity weight (default 1).
    %                      Velocity residuals are multiplied by this
    %                      weight times median(diff(t)).
    %
    %   The fit minimizes the weighted position and velocity residuals in
    %   each window while holding previously committed coefficients fixed.
    %   It is a least-squares approximation, not exact interpolation or a
    %   globally optimal fit. All node times are needed up front. Lookahead
    %   samples are used and there is no online append interface. The
    %   default quintic has continuous acceleration and jerk at all knots.
    %
    %   sp is a scalar struct containing the shifted nodes and knot vector,
    %   coefficients (one column per component), spline order and degree,
    %   and the indices of window starts. Treat its contents as read-only.
    %   Use coa.interp.evaluateSpline to evaluate the result. The interval
    %   [t(1), t(end)] is the valid evaluation domain.
    %
    %   Example:
    %     sp = coa.interp.fitWindowedSpline(t, position, velocity);
    %     derivatives = coa.interp.evaluateSpline(sp, query_times, 3);
    %     jerk = derivatives(:, :, 4);

    if nargin < 4
        opts = struct();
    end

    validateattributes(t, {'double'}, ...
                       {'real', 'finite', 'vector', 'nonempty'}, mfilename, 't');
    validateattributes(position, {'double'}, ...
                       {'real', 'finite', '2d', 'nonempty'}, mfilename, 'position');
    validateattributes(velocity, {'double'}, ...
                       {'real', 'finite', '2d', 'nonempty'}, mfilename, 'velocity');
    defaults = struct('order', 6, 'window_size', 40, 'stride', 20, ...
                      'velocity_weight', 1);
    if ~isstruct(opts) || ~isscalar(opts) || ...
            ~all(ismember(fieldnames(opts), fieldnames(defaults)))
        error('coa:interp:fitWindowedSpline:InvalidOptions', ...
              ['opts must be a scalar struct with only order, window_size, ' ...
               'stride, and velocity_weight.']);
    end
    names = fieldnames(defaults);
    for i = 1:numel(names)
        if ~isfield(opts, names{i})
            opts.(names{i}) = defaults.(names{i});
        end
    end
    validateattributes(opts.order, {'double'}, ...
                       {'real', 'finite', 'scalar', 'integer', '>=', 2});
    validateattributes(opts.window_size, {'double'}, ...
                       {'real', 'finite', 'scalar', 'integer', '>=', opts.order});
    validateattributes(opts.stride, {'double'}, ...
                       {'real', 'finite', 'scalar', 'integer', ...
                        '>=', 1, '<', opts.window_size});
    validateattributes(opts.velocity_weight, {'double'}, ...
                       {'real', 'finite', 'scalar', 'positive'});

    t = t(:);
    n = numel(t);
    d = size(position, 2);
    k = opts.order;
    p = k - 1;
    if size(position, 1) ~= n || ~isequal(size(position), size(velocity))
        error('coa:interp:fitWindowedSpline:SizeMismatch', ...
              'Position and velocity must both be N-by-D arrays.');
    end
    if n < k || any(diff(t) <= 0)
        error('coa:interp:fitWindowedSpline:InvalidTimes', ...
              'Times must be strictly increasing with at least order samples.');
    end

    t0 = t(1);
    nodes = t - t0;
    knots = [repmat(nodes(1), 1, k), nodes(2:n - 1).', ...
             repmat(nodes(n), 1, k)];
    n_coefficients = n + k - 2;
    velocity_weight = opts.velocity_weight * median(diff(nodes));
    coefficients = zeros(n_coefficients, d);
    last_committed = 0;
    window_starts = zeros(ceil(n / opts.stride), 1);
    n_windows = 0;
    a = 1;

    while true
        b = min(a + opts.window_size - 1, n);
        is_last = b == n;
        cols = a:b + k - 2;
        n_window = b - a + 1;
        values = zeros(n_window, numel(cols));
        derivatives = zeros(n_window, numel(cols));
        for r = 1:n_window
            i = a + r - 1;
            j = min(i, n - 1);
            basis = splineBasisDerivatives(p + j - 1, nodes(i), p, 1, knots);
            local_cols = (j:j + p) - cols(1) + 1;
            keep = local_cols <= numel(cols);
            values(r, local_cols(keep)) = basis(1, keep);
            derivatives(r, local_cols(keep)) = basis(2, keep);
        end

        matrix = [values; velocity_weight * derivatives];
        rhs = [position(a:b, :); velocity_weight * velocity(a:b, :)];
        frozen = cols <= last_committed;
        rhs = rhs - matrix(:, frozen) * coefficients(cols(frozen), :);
        coefficients(cols(~frozen), :) = matrix(:, ~frozen) \ rhs;

        if is_last
            last_committed = n_coefficients;
        else
            last_committed = a + opts.stride + k - 2;
        end
        n_windows = n_windows + 1;
        window_starts(n_windows) = a;
        if is_last
            break
        end
        a = a + opts.stride;
    end

    sp = struct('t0', t0, 'nodes', nodes, 'knots', knots, 'order', k, ...
                'degree', p, 'coefficients', coefficients, ...
                'window_starts', window_starts(1:n_windows));
end
