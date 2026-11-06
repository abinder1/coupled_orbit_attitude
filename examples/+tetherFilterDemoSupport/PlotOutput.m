% Copyright (c) 2026 Andrew Binder

classdef PlotOutput
    % Validate plot modes and manage displayed or saved figures.
    methods (Static)

        function keys = plotKeys()
            % Every plot that can be switched between "off", "save", and "show".
            keys = ["prediction_history", "trajectory", "anees", "rmse", ...
                    "innovations", ...
                    "residuals", "corner", "relative_corner"];
        end

        function preparePlotOutput(config)
            % Fail before the long run on a mistyped plot setting; create the output
            % folder if any plot is saved.
            import tetherFilterDemoSupport.*
            given = string(fieldnames(config.plot_modes))';
            unknown = setdiff(given, PlotOutput.plotKeys());
            missing = setdiff(PlotOutput.plotKeys(), given);
            if ~isempty(unknown) || ~isempty(missing)
                error('tetherFilterDemo:PlotModeKeys', ...
                      ['config.plot_modes has unknown keys [%s] and lacks ' ...
                       'keys [%s].'], ...
                      strjoin(unknown, ', '), strjoin(missing, ', '));
            end

            modes = string(struct2cell(config.plot_modes));
            invalid = ~ismember(modes, ["off", "save", "show"]);
            if any(invalid)
                error('tetherFilterDemo:PlotModeValue', ...
                      'Plot modes must be "off", "save", or "show"; invalid: %s.', ...
                      strjoin(given(invalid'), ', '));
            end
            if any(modes == "save") && ~isfolder(config.plot_output_folder)
                mkdir(config.plot_output_folder);
            end
        end

        function mode = plotMode(config, key)
            % Output mode ("off", "save", or "show") configured for one plot.
            mode = string(config.plot_modes.(key));
        end

        function fig = newFigure(config, key, varargin)
            % Create a figure, hidden when it will only be saved. Rendering offscreen
            % avoids the cost of drawing a large figure to the display.
            import tetherFilterDemoSupport.*
            visibility = 'on';
            if PlotOutput.plotMode(config, key) == "save"
                visibility = 'off';
            end
            fig = figure(varargin{:}, 'Visible', visibility);
        end

        function finishFigure(fig, config, key)
            % Write a saved figure to file and close it; leave shown figures open.
            import tetherFilterDemoSupport.*
            if PlotOutput.plotMode(config, key) ~= "save"
                return
            end
            stem = strip(regexprep(string(fig.Name), '[^A-Za-z0-9]+', '_'), '_');
            file = fullfile(config.plot_output_folder, ...
                            key + "_" + stem + "." + config.plot_file_format);
            exportgraphics(fig, file, 'Resolution', config.plot_resolution);
            close(fig);
            fprintf('  Saved %s\n', file);
        end

    end
end
