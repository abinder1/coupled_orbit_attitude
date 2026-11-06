% Copyright (c) 2026 Andrew Binder

classdef Visualization
    % Reference-filter summaries and truth-trajectory visualizations.
    methods (Static)

        function reportFilterRun(filter_run, config)
            % Report prior/posterior errors and uncertainty for every measurement.
            import tetherFilterDemoSupport.*
            updates = filter_run.updates;
            fprintf(['EKF demonstration: %.0f-second propagation with %d evenly ' ...
                     'spaced measurement updates\n'], config.duration_s, ...
                    config.measurement_count);
            fprintf('Measurements: all implemented coa.meas measurement groups\n');
            fprintf(['Update  Time [min]  Innovation norm  Pos. error [km]  ' ...
                     'Vel. error [m/s]\n']);

            for k = 1:numel(updates)
                prior_error = updates(k).prior_state - updates(k).truth_state;
                posterior_error = updates(k).posterior_state - updates(k).truth_state;
                pos_rows = config.position_rows;
                vel_rows = config.velocity_rows;
                fprintf('%6d  %10.2f  %15.3g  %14.3g -> %-8.3g  %14.3g -> %-8.3g\n', ...
                        k, config.measurement_times(k) * config.time_scale / 60, ...
                        norm(updates(k).innovation), ...
                        Visualization.stateNorm(prior_error, pos_rows, config, ...
                                                'position'), ...
                        Visualization.stateNorm(posterior_error, pos_rows, config, ...
                                                'position'), ...
                        Visualization.stateNorm(prior_error, vel_rows, config, ...
                                                'velocity'), ...
                        Visualization.stateNorm(posterior_error, vel_rows, config, ...
                                                'velocity'));
            end
        end

        function plotPredictionHistory(history, config)
            % Plot signed XYZ errors and shaded componentwise +/-3-sigma bounds.
            import tetherFilterDemoSupport.*
            if PlotOutput.plotMode(config, 'prediction_history') == "off"
                return
            end
            time_minutes = history.time * config.time_scale / 60;
            component_colors = lines(3);
            component_names = ["x", "y", "z"];
            position_limit = ...
                Visualization.symmetricPlotLimit(history.position_error, ...
                                                 history.position_sigma);
            velocity_limit = ...
                Visualization.symmetricPlotLimit(history.velocity_error, ...
                                                 history.velocity_sigma);

            for body = 1:3
                body_rows = 3 * (body - 1) + (1:3);
                figure_name = config.body_names(body) + " signed errors";
                fig = PlotOutput.newFigure(config, 'prediction_history', 'Name', ...
                                           figure_name, 'Color', 'w');
                layout = tiledlayout(3, 2, 'TileSpacing', 'compact', ...
                                     'Padding', 'compact');

                for component = 1:3
                    row = body_rows(component);
                    position_axes = nexttile(layout, 2 * component - 1);
                    position_error = history.position_error(row, :);
                    position_sigma = history.position_sigma(row, :);
                    Visualization.plotSignedComponent( ...
                                                      position_axes, time_minutes, ...
                                                      position_error, ...
                                                      position_sigma, config, ...
                                                      component_colors(component, :));
                    ylim(position_axes, position_limit);
                    title(position_axes, component_names(component) + " position");
                    ylabel(position_axes, 'Position error [km]');

                    velocity_axes = nexttile(layout, 2 * component);
                    velocity_error = history.velocity_error(row, :);
                    velocity_sigma = history.velocity_sigma(row, :);
                    Visualization.plotSignedComponent( ...
                                                      velocity_axes, time_minutes, ...
                                                      velocity_error, ...
                                                      velocity_sigma, config, ...
                                                      component_colors(component, :));
                    ylim(velocity_axes, velocity_limit);
                    title(velocity_axes, component_names(component) + " velocity");
                    ylabel(velocity_axes, 'Velocity error [m/s]');
                end

                xlabel(layout, 'Time [min]');
                title(layout, config.body_names(body) + ...
                      " signed component errors with shaded \pm3\sigma bounds");
                PlotOutput.finishFigure(fig, config, 'prediction_history');
            end
        end

        function limits = symmetricPlotLimit(errors, sigmas)
            % Include all signed errors and three-sigma envelopes in a shared range.
            extent = max([abs(errors(:)); 3 * sigmas(:)]);
            extent = max(extent, eps(extent));
            limits = 1.05 * [-extent, extent];
        end

        function plotSignedComponent(ax, time, error, sigma, config, color)
            % Draw one signed component error and its +/-3-sigma envelope.
            hold(ax, 'on');

            bound = 3 * sigma;
            fill(ax, [time, fliplr(time)], [bound, -fliplr(bound)], ...
                 color, 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');
            boundary_color = [0.45, 0.45, 0.45];
            plot(ax, time, bound, '-', 'Color', boundary_color, 'LineWidth', 0.5, ...
                 'HandleVisibility', 'off');
            plot(ax, time, -bound, '-', 'Color', boundary_color, 'LineWidth', 0.5, ...
                 'HandleVisibility', 'off');

            update_times = config.measurement_times * config.time_scale / 60;
            xline(ax, update_times, ':', 'Color', [0.65, 0.65, 0.65], ...
                  'HandleVisibility', 'off');

            plot(ax, time, error, 'Color', color, 'LineWidth', 1.1);

            yline(ax, 0, 'k-', 'HandleVisibility', 'off');
            grid(ax, 'on');
        end

        function plotTruthTrajectory(scenario, config)
            % Plot the nominal truth paths, payload-to-catch range, and inbound speed.
            import tetherFilterDemoSupport.*
            if PlotOutput.plotMode(config, 'trajectory') == "off"
                return
            end

            sample_count = config.measurement_count * ...
                           config.plot_samples_per_interval + 1;
            sample_time = linspace(0, config.duration, sample_count);
            truth = deval(scenario.truth_solution, sample_time);
            position = truth(1:9, :) * config.length_scale;
            velocity = truth(10:18, :) * config.length_scale / config.time_scale;
            time_minutes = sample_time * config.time_scale / 60;

            relative_position = position(7:9, :) - position(4:6, :);
            relative_velocity = velocity(7:9, :) - velocity(4:6, :);
            payload_to_catch_range = vecnorm(relative_position);
            inbound_speed = -sum(relative_position .* relative_velocity, 1) ./ ...
                            payload_to_catch_range * 1000;

            fig = PlotOutput.newFigure(config, 'trajectory', 'Name', ...
                                       'Truth motion and payload range', 'Color', ...
                                       'w', ...
                                       'Position', [60, 60, 1200, 720]);
            layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', ...
                                 'Padding', 'compact');

            trajectory_axes = nexttile(layout, 1, [2, 1]);
            hold(trajectory_axes, 'on');
            colors = lines(3);
            bodies = {1:3, 4:6, 7:9};
            names = ["Counterweight (m1)", "Catch mechanism (m2)", ...
                     "Payload (m3)"];
            for body = 1:numel(bodies)
                rows = bodies{body};
                plot3(trajectory_axes, position(rows(1), :), position(rows(2), :), ...
                      position(rows(3), :), 'Color', colors(body, :), ...
                      'LineWidth', 1.2);
            end
            grid(trajectory_axes, 'on');
            axis(trajectory_axes, 'equal');
            view(trajectory_axes, 3);
            xlabel(trajectory_axes, 'X [km]');
            ylabel(trajectory_axes, 'Y [km]');
            zlabel(trajectory_axes, 'Z [km]');
            title(trajectory_axes, 'Truth trajectories in the rotating frame');
            legend(trajectory_axes, names, 'Location', 'best');

            range_axes = nexttile(layout, 2);
            plot(range_axes, time_minutes, payload_to_catch_range, 'LineWidth', 1.2);
            grid(range_axes, 'on');
            ylabel(range_axes, 'Range |r_3 - r_2| [km]');
            title(range_axes, 'Payload (m3) to catch mechanism (m2)');

            speed_axes = nexttile(layout, 4);
            plot(speed_axes, time_minutes, inbound_speed, 'LineWidth', 1.2);
            yline(speed_axes, 0, 'k-', 'HandleVisibility', 'off');
            grid(speed_axes, 'on');
            xlabel(speed_axes, 'Time [min]');
            ylabel(speed_axes, 'Inbound speed [m/s]');
            title(speed_axes, 'Positive when closing');

            sgtitle(layout, 'Truth motion and payload-to-catch range history');
            PlotOutput.finishFigure(fig, config, 'trajectory');
        end

        function value = stateNorm(state_error, rows, config, quantity)
            % Convert a state-error norm to km or m/s.
            value = norm(state_error(rows));
            if strcmp(quantity, 'position')
                value = value * config.length_scale;
            else
                value = value * config.length_scale / config.time_scale * 1000;
            end
        end

    end
end
