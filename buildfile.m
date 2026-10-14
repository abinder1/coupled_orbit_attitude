% Copyright (c) 2026 Andrew Binder
function plan = buildfile
    % Build plan for coupled_orbit_attitude.
    %   buildtool          -> setup, check, test (default)
    %   buildtool test     -> same as above
    %   buildtool clean    -> delete build outputs (results/)
    % Requires MATLAB R2024a or later.
    import matlab.buildtool.tasks.*

    plan = buildplan(localfunctions);

    plan("clean") = CleanTask;
    plan("check") = CodeIssuesTask("toolbox");
    plan("test")  = TestTask("tests", ...
                             SourceFiles         = "toolbox", ...
                             TestResults         = "test_results/test-results.xml", ...
                             CodeCoverageResults = "test_results/coverage.xml");

    plan("check").Dependencies = "setup";
    plan("test").Dependencies  = ["setup", "check"];
    plan.DefaultTasks = "test";
end

function setupTask(~)
    % Put the project root (and all subfolders) on the MATLAB path.
    root = fileparts(mfilename("fullpath")); % folder holding buildfile (project root)
    addpath(genpath(root));
end
