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
    TestResults         = "results/test-results.xml", ...
    CodeCoverageResults = "results/coverage.xml");

plan("check").Dependencies = "setup";
plan("test").Dependencies  = ["setup", "check"];
plan.DefaultTasks = "test";
end

function setupTask(~)
% Put the toolbox on the MATLAB path.
addpath(fullfile(fileparts(mfilename("fullpath")), "toolbox"));
end
