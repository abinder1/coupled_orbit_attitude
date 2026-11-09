% Copyright (c) 2026 Andrew Binder
function setupProjectPaths(root)
    % Add the project and its pinned dependency's production code to the path.
    % Run git submodule update --init --recursive before calling this function.
    if nargin == 0
        root = fileparts(mfilename('fullpath'));
    end

    dependency = fullfile(root, 'external', 'ephem_integrator', 'toolbox');
    fit_file = fullfile(dependency, '+emframe', '+interp', 'fitWindowedSpline.m');
    eval_file = fullfile(dependency, '+emframe', '+interp', 'evaluateSpline.m');
    if ~isfile(fit_file) || ~isfile(eval_file)
        error('coa:setupProjectPaths:MissingDependency', ...
              ['The ephem_integrator submodule is missing or incomplete. ', ...
               'From the repository root, run ', ...
               'git submodule update --init --recursive.']);
    end

    addpath(fullfile(root, 'toolbox'), fullfile(root, 'examples'), dependency);
end
