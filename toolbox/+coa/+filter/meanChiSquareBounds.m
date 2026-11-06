% Copyright (c) 2026 Andrew Binder

function [lower, upper] = meanChiSquareBounds(run_count, dof, confidence)
    % MEANCHISQUAREBOUNDS Two-sided bounds for a mean chi-square statistic.
    %   [lower, upper] = coa.filter.meanChiSquareBounds(run_count, dof, ...
    %       confidence) returns the two-sided acceptance bounds for the mean
    %   of run_count independent chi-square statistics, each with dof degrees
    %   of freedom. dof may be a scalar or matrix; both outputs have its size.
    %   confidence is strictly between 0 and 1.
    arguments
        run_count (1, 1) double {mustBeFinite, mustBeInteger, mustBePositive}
        dof (:, :) double {mustBeNonempty, mustBeFinite, ...
                           mustBeInteger, mustBePositive}
        confidence (1, 1) double {mustBeFinite, ...
                                  mustBeGreaterThan(confidence, 0), ...
                                  mustBeLessThan(confidence, 1)}
    end

    tail = (1 - confidence) / 2;
    total_dof = run_count * dof;
    lower = 2 * gammaincinv(tail, total_dof / 2) / run_count;
    upper = 2 * gammaincinv(1 - tail, total_dof / 2) / run_count;
end
