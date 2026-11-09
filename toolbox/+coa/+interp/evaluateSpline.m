% Copyright (c) 2026 Andrew Binder

function values = evaluateSpline(sp, tq, nd, side)
    % EVALUATESPLINE Evaluate a fitted spline and its time derivatives.
    %   values = coa.interp.evaluateSpline(sp, tq) evaluates position at
    %   query times tq for the struct returned by coa.interp.fitWindowedSpline.
    %
    %   values = coa.interp.evaluateSpline(sp, tq, nd) returns a
    %   numel(tq)-by-D-by-(nd+1) array. Slice 1 holds position and slice
    %   m+1 holds the m-th time derivative in the original time units.
    %   nd defaults to 0 and must be an integer no greater than the degree.
    %   tq must be a real finite vector within the fitted interval, or empty.
    %
    %   values = coa.interp.evaluateSpline(sp, tq, nd, side) selects the
    %   'right' (default) or 'left' polynomial piece at interior knots. The
    %   first and last endpoints use the only available adjacent piece.

    if nargin < 3
        nd = 0;
    end
    if nargin < 4
        side = 'right';
    end
    validateattributes(tq, {'double'}, {'real', 'finite', '2d'}, mfilename, 'tq');
    if ~isempty(tq) && ~isvector(tq)
        error('coa:interp:evaluateSpline:InvalidQuerySize', ...
              'Query times must be a vector.');
    end
    validateattributes(nd, {'double'}, ...
                       {'real', 'finite', 'scalar', 'integer', ...
                        '>=', 0, '<=', sp.degree}, mfilename, 'nd');
    if ~(ischar(side) && isrow(side) && ...
         (strcmp(side, 'left') || strcmp(side, 'right')))
        error('coa:interp:evaluateSpline:InvalidSide', ...
              'side must be ''left'' or ''right''.');
    end

    u = tq(:) - sp.t0;
    if any(u < sp.nodes(1) | u > sp.nodes(end))
        error('coa:interp:evaluateSpline:OutOfRange', ...
              'Query times must lie in the fitted interval.');
    end
    d = size(sp.coefficients, 2);
    values = zeros(numel(u), d, nd + 1);
    if isempty(u)
        return
    end
    spans = discretize(u, sp.nodes);
    for q = 1:numel(u)
        j = spans(q);
        if strcmp(side, 'left') && j > 1 && u(q) == sp.nodes(j)
            j = j - 1;
        end
        basis = splineBasisDerivatives(sp.degree + j - 1, u(q), ...
                                       sp.degree, nd, sp.knots);
        values(q, :, :) = permute(basis * sp.coefficients(j:j + sp.degree, :), ...
                                  [3, 2, 1]);
    end
end
