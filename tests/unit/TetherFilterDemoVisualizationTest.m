% Copyright (c) 2026 Andrew Binder

classdef TetherFilterDemoVisualizationTest < matlab.unittest.TestCase
    methods (TestMethodSetup)

        function addExamplePath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            original_path = path;
            addpath(fullfile(root, 'examples'));
            tc.addTeardown(@path, original_path);
        end

    end

    methods (Test)

        function plotsComponentMagnitudesAndBounds(tc)
            import tetherFilterDemoSupport.*
            history.time = [0, 1, 1, 2];
            history.position_error = (1:9)' * [-1e4, -1, 0, 1e-5];
            history.velocity_error = (1:9)' * [100, -0.1, 0, -1e-7];
            history.position_sigma = (1:9)' * [100, 1, 0, 1e-4];
            history.velocity_sigma = (1:9)' * [1, 0.1, 0, 1e-6];
            config.time_scale = 60;
            config.measurement_times = 1;
            config.body_names = ["Test body 1", "Test body 2", "Test body 3"];
            config.plot_modes.prediction_history = "show";

            existing = findall(groot, 'Type', 'figure');
            tc.addTeardown(@closeNewFigures, existing);
            Visualization.plotPredictionHistory(history, config);
            figures = setdiff(findall(groot, 'Type', 'figure'), existing);
            tc.assertNumElements(figures, 3);

            for body = 1:3
                fig = figures(strcmp(string({figures.Name}), ...
                                     config.body_names(body) + " error magnitudes"));
                tc.assertNumElements(fig, 1);
                axes_list = findall(fig, 'Type', 'axes');
                tc.assertNumElements(axes_list, 6);
                for ax = reshape(axes_list, 1, [])
                    tile = ax.Layout.Tile;
                    row = 3 * (body - 1) + ceil(tile / 2);
                    if mod(tile, 2) == 1
                        errors = history.position_error(row, :);
                        sigma = history.position_sigma(row, :);
                    else
                        errors = history.velocity_error(row, :);
                        sigma = history.velocity_sigma(row, :);
                    end

                    tc.verifyEqual(ax.XScale, 'linear');
                    tc.verifyEqual(ax.YScale, 'log');
                    lines = findall(ax, 'Type', 'line');
                    tc.assertNumElements(lines, 2);
                    error_line = lines([lines.LineWidth] == 1.1);
                    bound_line = lines([lines.LineWidth] == 0.5);
                    tc.assertNumElements(error_line, 1);
                    tc.assertNumElements(bound_line, 1);
                    tc.verifyEqual(error_line.XData, history.time);
                    tc.verifyEqual(error_line.YData, abs(errors));
                    tc.verifyEqual(bound_line.YData, 3 * sigma);
                    tc.verifyEqual(bound_line.LineStyle, '-');
                    tc.verifyEqual(bound_line.Color, [0.45, 0.45, 0.45]);
                    shade = findall(ax, 'Type', 'patch');
                    tc.assertNumElements(shade, 1);
                    tc.verifyEqual(shade.XData(:)', ...
                                   [history.time, fliplr(history.time)]);
                    tc.verifyEqual(shade.YData(:)', ...
                                   [max(3 * sigma, ax.YLim(1)), ...
                                    ax.YLim(1) * ones(size(history.time))]);
                    tc.verifyEqual(shade.EdgeColor, 'none');
                    tc.verifyLessThanOrEqual(ax.YLim(1), min(abs(errors(errors ~= 0))));
                    tc.verifyGreaterThanOrEqual(ax.YLim(2), max(abs(errors)));
                end
            end
        end

        function disabledHistoryDoesNotCreateFigures(tc)
            import tetherFilterDemoSupport.*
            existing = findall(groot, 'Type', 'figure');
            config.plot_modes.prediction_history = "off";
            Visualization.plotPredictionHistory(struct(), config);
            tc.verifyEqual(findall(groot, 'Type', 'figure'), existing);
        end

        function savesAndClosesHistoryFigures(tc)
            import tetherFilterDemoSupport.*
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            history.time = [0, 1];
            history.position_error = repmat([-2, 1], 9, 1);
            history.velocity_error = repmat([2, -1], 9, 1);
            history.position_sigma = ones(9, 2);
            history.velocity_sigma = ones(9, 2);
            config.time_scale = 60;
            config.measurement_times = 1;
            config.body_names = ["Test body 1", "Test body 2", "Test body 3"];
            config.plot_modes.prediction_history = "save";
            config.plot_output_folder = folder.Folder;
            config.plot_file_format = "png";
            config.plot_resolution = 50;

            existing = findall(groot, 'Type', 'figure');
            tc.addTeardown(@closeNewFigures, existing);
            Visualization.plotPredictionHistory(history, config);
            tc.verifyEqual(findall(groot, 'Type', 'figure'), existing);
            images = dir(fullfile(folder.Folder, 'prediction_history_*.png'));
            tc.assertNumElements(images, 3);
            tc.verifyGreaterThan([images.bytes], 0);
        end

        function marksMeasurementTimesOnRangePlot(tc)
            import tetherFilterDemoSupport.*
            state = (1:18)';
            scenario.truth_solution = ode113(@(~, x) ones(size(x)), [0, 2], state);
            config.plot_modes.trajectory = "show";
            config.duration = 2;
            config.measurement_count = 2;
            config.measurement_times = [0.5, 1.5];
            config.plot_samples_per_interval = 3;
            config.time_scale = 120;
            config.length_scale = 1;
            existing = findall(groot, 'Type', 'figure');
            tc.addTeardown(@closeNewFigures, existing);

            Visualization.plotTruthTrajectory(scenario, config);
            fig = setdiff(findall(groot, 'Type', 'figure'), existing);
            tc.assertNumElements(fig, 1);
            axes_list = findall(fig, 'Type', 'axes');
            range_ax = axes_list(arrayfun(@(ax) ax.Layout.Tile == 2, axes_list));
            tc.assertNumElements(range_ax, 1);
            markers = findall(range_ax, 'Type', 'constantline');
            tc.assertNumElements(markers, 2);
            tc.verifyEqual(sort([markers.Value]), [1, 3]);
            tc.verifyEqual({markers.LineStyle}, {'--', '--'});
            tc.verifyEqual([markers.LineWidth], [0.5, 0.5]);
            tc.verifyEqual(numel(findall(range_ax, 'Type', 'line')), 1);
        end

        function rmseUsesLogarithmicYAxes(tc)
            import tetherFilterDemoSupport.*
            monte_carlo.errors = reshape(1:108, [2, 3, 18]) / 100;
            monte_carlo.covariances = repmat(eye(18), 1, 1, 3);
            monte_carlo.time_minutes = [1, 2, 3];
            monte_carlo.run_count = 2;
            config.length_scale = 2;
            config.time_scale = 4;
            config.body_names = ["A", "B", "C"];
            config.plot_modes.rmse = "show";
            existing = findall(groot, 'Type', 'figure');
            tc.addTeardown(@closeNewFigures, existing);

            Analysis.plotRmse(monte_carlo, config);
            fig = setdiff(findall(groot, 'Type', 'figure'), existing);
            tc.assertNumElements(fig, 1);
            axes_list = findall(fig, 'Type', 'axes');
            tc.assertNumElements(axes_list, 18);
            for ax = reshape(axes_list, 1, [])
                tile = ax.Layout.Tile;
                row = ceil(tile / 3);
                body = mod(tile - 1, 3) + 1;
                is_velocity = row > 3;
                state = 9 * is_velocity + 3 * (body - 1) + row - 3 * is_velocity;
                scale = config.length_scale;
                if is_velocity
                    scale = scale / config.time_scale * 1000;
                end
                expected_rmse = sqrt(mean(monte_carlo.errors(:, :, state).^2, 1));
                tc.verifyEqual(ax.XScale, 'linear');
                tc.verifyEqual(ax.YScale, 'log');
                lines = findall(ax, 'Type', 'line');
                tc.assertNumElements(lines, 2);
                rmse_line = lines(strcmp({lines.Marker}, 'o'));
                sigma_line = lines(strcmp({lines.LineStyle}, '--'));
                tc.assertNumElements(rmse_line, 1);
                tc.assertNumElements(sigma_line, 1);
                tc.verifyEqual(rmse_line.YData, expected_rmse * scale);
                tc.verifyEqual(sigma_line.YData, ones(1, 3) * scale);
            end
        end

        function relativeCornerUsesMetersPerSecond(tc)
            import tetherFilterDemoSupport.*
            state = zeros(18, 1);
            state(4) = 1;
            state(14) = 1;
            endpoint.reference_state = state;
            endpoint.covariance = eye(18);
            endpoint.errors = reshape(1:144, [8, 18]) / 100;
            endpoint.label = 'first';
            second = endpoint;
            second.label = 'second';
            monte_carlo.endpoints = [endpoint, second];
            monte_carlo.run_count = 8;
            config.length_scale = 2;
            config.time_scale = 4;
            config.corner_ellipse_sigma = 3;
            config.plot_modes.relative_corner = "show";
            relative = Analysis.relativeStateDistribution(endpoint, config);
            expected_sigma = sqrt(diag(relative.covariance)) / 1000;
            existing = findall(groot, 'Type', 'figure');
            tc.addTeardown(@closeNewFigures, existing);

            Analysis.plotRelativeCorner(monte_carlo, config);
            figures = setdiff(findall(groot, 'Type', 'figure'), existing);
            tc.assertNumElements(figures, 2);
            for fig = reshape(figures, 1, [])
                axes_list = findall(fig, 'Type', 'axes');
                texts = strings(0);
                for ax = reshape(axes_list, 1, [])
                    texts = [texts; string(ax.XLabel.String(:)); ...
                             string(ax.YLabel.String(:))]; %#ok<AGROW>
                end
                for index = 4:6
                    label = sprintf('\\sigma=%.2g m/s', expected_sigma(index));
                    tc.verifyTrue(any(texts == label));
                end
                tc.verifyFalse(any(contains(texts, 'mm/s')));
            end
        end

    end
end

function closeNewFigures(existing)
    close(setdiff(findall(groot, 'Type', 'figure'), existing));
end
