% Copyright (c) 2026 Andrew Binder

function nees = computeNees(errors, covariances, blocks)
    % COMPUTENEES Normalized estimation error squared for trials and blocks.
    %   nees = coa.filter.computeNees(errors, covariances, blocks) computes
    %   e' * (P \ e) for each trial, update, and caller-defined state block.
    %
    %   errors has size [run, update, state]. covariances has size
    %   [state, state, update, run], with each trial's own covariance.
    %   blocks is a nonempty struct array with a nonempty vector of state
    %   indices in each rows field. The output has size [run, update, block].
    %   Each selected covariance block is expected to be nonsingular; no
    %   symmetrization or other covariance repair is performed.
    arguments
        errors double {mustBeFinite}
        covariances double {mustBeFinite}
        blocks struct
    end

    [run_count, step_count, state_count] = size(errors);
    if ndims(errors) > 3 || run_count == 0 || step_count == 0 || state_count == 0
        error('coa:filter:computeNees:InvalidErrors', ...
              'errors must be a nonempty [run, update, state] array.');
    end
    if ndims(covariances) > 4 || ...
            ~isequal(size(covariances, 1), state_count) || ...
            ~isequal(size(covariances, 2), state_count) || ...
            ~isequal(size(covariances, 3), step_count) || ...
            ~isequal(size(covariances, 4), run_count)
        error('coa:filter:computeNees:SizeMismatch', ...
              'covariances must have size [state, state, update, run].');
    end
    if isempty(blocks) || ~isfield(blocks, 'rows')
        error('coa:filter:computeNees:InvalidBlocks', ...
              'blocks must be a nonempty struct array with a rows field.');
    end
    for b = 1:numel(blocks)
        rows = blocks(b).rows;
        if ~isnumeric(rows) || ~isreal(rows) || ~isvector(rows) || ...
                isempty(rows) || any(~isfinite(rows)) || ...
                any(rows ~= fix(rows)) || any(rows < 1 | rows > state_count)
            error('coa:filter:computeNees:InvalidBlocks', ...
                  'Each block must contain valid state indices in rows.');
        end
    end

    nees = zeros(run_count, step_count, numel(blocks));
    for run = 1:run_count
        for k = 1:step_count
            for b = 1:numel(blocks)
                rows = blocks(b).rows;
                e = squeeze(errors(run, k, rows));
                nees(run, k, b) = e' * (covariances(rows, rows, k, run) \ e);
            end
        end
    end
end
