% Copyright (c) 2026 Andrew Binder

classdef Analysis
    % Reports and visualizations derived from Monte Carlo filter results.
    methods (Static)

        function analyses = monteCarloAnalyses()
            % Registry of Monte Carlo analyses, run in order. To add one, write a
            % function taking (monte_carlo, config) -- a report, a plot, or both -- and
            % append its handle(s) here.
            import tetherFilterDemoSupport.*
            analyses = {@Analysis.reportAnees, @Analysis.reportNis, ...
                        @Analysis.plotAnees, @Analysis.plotCorner, ...
                        @Analysis.plotInnovations, @Analysis.plotResiduals, ...
                        @Analysis.plotRmse, @Analysis.reportRelativeState, ...
                        @Analysis.plotRelativeCorner};
        end

        function runMonteCarloAnalyses(analyses, monte_carlo, config)
            % Run every registered analysis on the shared Monte Carlo results.
            for analysis = analyses
                analysis{1}(monte_carlo, config);
            end
        end

        function reportAnees(monte_carlo, config)
            % Summarize ANEES against its acceptance region for each block.
            fprintf('\nANEES over %d full-EKF Monte Carlo runs (%.0f%% region)\n', ...
                    monte_carlo.run_count, 100 * config.nees_confidence);
            fprintf('%-20s %5s  %9s  %-19s  %s\n', 'Block', 'dof', 'Mean ANEES', ...
                    'Acceptance region', 'Steps inside');

            for b = 1:numel(monte_carlo.blocks)
                anees = monte_carlo.anees(:, b);
                lower = monte_carlo.bounds.lower(b);
                upper = monte_carlo.bounds.upper(b);
                inside = mean(anees >= lower & anees <= upper);
                fprintf('%-20s %5d  %9.3g  [%7.3g, %7.3g]  %5.1f%%\n', ...
                        monte_carlo.blocks(b).name, monte_carlo.blocks(b).dof, ...
                        mean(anees), lower, upper, 100 * inside);
            end
        end

        function reportNis(monte_carlo, config)
            % Average NIS per update against its chi-square acceptance region.
            [bounds.lower, bounds.upper] = ...
                coa.filter.meanChiSquareBounds(monte_carlo.run_count, ...
                                               monte_carlo.nis_dof, ...
                                               config.nees_confidence);
            fprintf('\nANIS per update (dof %d, region [%.1f, %.1f])\n', ...
                    monte_carlo.nis_dof, bounds.lower, bounds.upper);
            for k = 1:numel(monte_carlo.anis)
                fprintf('  t = %7.1f min: ANIS = %9.3g\n', ...
                        monte_carlo.time_minutes(k), ...
                        monte_carlo.anis(k));
            end
        end

        function plotAnees(monte_carlo, config)
            % Plot per-body ANEES versus time with the chi-square acceptance region.
            import tetherFilterDemoSupport.*
            if PlotOutput.plotMode(config, 'anees') == "off"
                return
            end
            fig = PlotOutput.newFigure(config, 'anees', 'Name', ...
                                       'ANEES consistency', 'Color', 'w');
            body_blocks = 2:numel(monte_carlo.blocks);
            layout = tiledlayout(1, numel(body_blocks), ...
                                 'TileSpacing', 'compact', 'Padding', 'compact');
            time = monte_carlo.time_minutes;
            colors = lines(numel(body_blocks));

            for panel = 1:numel(body_blocks)
                b = body_blocks(panel);
                ax = nexttile(layout);
                hold(ax, 'on');
                % Normalize by the block dimension so the expected value is 1.
                dof = monte_carlo.blocks(b).dof;
                lower = monte_carlo.bounds.lower(b) / dof;
                upper = monte_carlo.bounds.upper(b) / dof;
                fill(ax, [time, fliplr(time)], ...
                     [lower * ones(size(time)), upper * ones(size(time))], ...
                     [0.5, 0.5, 0.5], 'FaceAlpha', 0.2, ...
                     'EdgeColor', [0.45, 0.45, 0.45], ...
                     'LineWidth', 0.5);
                yline(ax, 1, 'k--');
                plot(ax, time, monte_carlo.anees(:, b) / dof, '-o', ...
                     'Color', colors(panel, :), 'LineWidth', 1.1, 'MarkerSize', 3);
                set(ax, 'YScale', 'log');
                grid(ax, 'on');
                title(ax, sprintf('%s (n = %d)', monte_carlo.blocks(b).name, ...
                                  monte_carlo.blocks(b).dof));
                ylabel(ax, 'ANEES / n');
            end

            xlabel(layout, 'Time [min]');
            title(layout, ...
                  sprintf(['Normalized ANEES over %d runs with %.0f%% acceptance ' ...
                           'region (dashed: expected value 1)'],  ...
                          monte_carlo.run_count, 100 * config.nees_confidence));
            PlotOutput.finishFigure(fig, config, 'anees');
        end

        function plotRmse(monte_carlo, config)
            % Dimensional RMSE over the Monte Carlo runs vs. the filter's predicted
            % sigma at each update. Columns are bodies; rows are position (km) and
            % velocity (m/s) axes. RMSE above the predicted sigma means overconfident.
            import tetherFilterDemoSupport.*
            if PlotOutput.plotMode(config, 'rmse') == "off"
                return
            end
            rmse = squeeze(sqrt(mean(monte_carlo.errors.^2, 1)));   % [update, state]
            sigma = sqrt(max(reshape(Analysis.diagBlocks(monte_carlo.covariances), ...
                                     18, []), 0))';
            position_scale = config.length_scale;                      % km
            velocity_scale = config.length_scale / config.time_scale * 1000;  % m/s
            scale = [position_scale * ones(1, 9), velocity_scale * ones(1, 9)];
            rmse = rmse .* scale;
            sigma = sigma .* scale;

            time = monte_carlo.time_minutes;
            axis_names = ["x", "y", "z"];
            fig = PlotOutput.newFigure(config, 'rmse', 'Name', ...
                                       'RMSE vs predicted sigma', ...
                                       'Color', 'w', 'Position', [80, 40, 1000, 900]);
            layout = tiledlayout(6, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
            for row = 1:6
                for body = 1:3
                    is_velocity = row > 3;
                    axis_index = row - 3 * is_velocity;
                    state = 9 * is_velocity + 3 * (body - 1) + axis_index;
                    ax = nexttile(layout);
                    hold(ax, 'on');
                    plot(ax, time, sigma(:, state), 'k--', 'LineWidth', 1.1);
                    plot(ax, time, rmse(:, state), '-o', 'LineWidth', 1.1, ...
                         'MarkerSize', 3, 'Color', [0, 0.447, 0.741]);
                    set(ax, 'FontSize', 7, 'YScale', 'log', 'XScale', 'linear');
                    grid(ax, 'on');
                    if body == 1
                        prefix = ["p", "v"];
                        unit = ["km", "m/s"];
                        ylabel(ax, sprintf('%s%s [%s]', prefix(1 + is_velocity), ...
                                           axis_names(axis_index), ...
                                           unit(1 + is_velocity)));
                    end
                    if row == 1
                        title(ax, config.body_names(body), 'FontSize', 9);
                    end
                    if row < 6
                        set(ax, 'XTickLabel', []);
                    end
                end
            end
            lg = legend(nexttile(1), {'predicted \sigma', 'RMSE'});
            lg.Layout.Tile = 'east';
            lg.FontSize = 8;
            xlabel(layout, 'Time [min]');
            title(layout, sprintf(['RMSE over %d runs vs. predicted \\sigma ' ...
                                   '(RMSE above \\sigma: overconfident)'], ...
                                  monte_carlo.run_count));
            PlotOutput.finishFigure(fig, config, 'rmse');
        end

        function d = diagBlocks(covariances)
            % Diagonals of every 18x18 page, returned as [18, 1, pages].
            pages = size(covariances, 3);
            d = zeros(18, 1, pages);
            for k = 1:pages
                d(:, 1, k) = diag(covariances(:, :, k));
            end
        end

        function plotInnovations(monte_carlo, config)
            % Pre-fit innovations, normalized by sqrt(diag(S)).
            import tetherFilterDemoSupport.*
            Analysis.plotMeasurementHistory(monte_carlo, config, 'innovation', ...
                                            'innovation_covariance', ...
                                            "Pre-fit innovations");
        end

        function plotResiduals(monte_carlo, config)
            % Post-fit residuals, normalized by sqrt(diag(R - H P+ H')), floored at
            % config.residual_sigma_floor times the sensor sigma.
            import tetherFilterDemoSupport.*
            Analysis.plotMeasurementHistory(monte_carlo, config, 'residual', ...
                                            'residual_covariance', ...
                                            "Post-fit residuals");
        end

        function plotMeasurementHistory(monte_carlo, config, value_field, ...
                                        covariance_field, ...
                                        title_text)
            % One tile per active sensor group. Values are normalized by their predicted
            % sigma, so every group shares the same axis and the +/-3 bounds mean the
            % same thing everywhere. Gray dots are all Monte Carlo runs (reused from
            % the consistency check); colored lines are the single reference filter run.
            % The plot key is the plural of value_field ('innovations' or 'residuals').
            import tetherFilterDemoSupport.*
            key = Analysis.plural(value_field);
            if PlotOutput.plotMode(config, key) == "off"
                return
            end
            updates = monte_carlo.reference_updates;
            sigma = sqrt(max(cell2mat(arrayfun(@(u) diag(u.(covariance_field)), ...
                                               updates', 'UniformOutput', false)), 0));
            % Post-fit variances R - H P+ H' can vanish (or round negative) for
            % near-exact measurements; floor sigma at a fraction of the sensor sigma so
            % the normalized values stay finite. Innovations are never affected, since
            % S >= R.
            noise_sigma = sqrt(cell2mat(arrayfun(@(u) sum(u.noise_chol.^2, 2), ...
                                                 updates', 'UniformOutput', false)));
            sigma = max(sigma, config.residual_sigma_floor * noise_sigma);
            reference = cell2mat(arrayfun(@(u) u.(value_field), updates', ...
                                          'UniformOutput', false)) ./ sigma;
            samples = monte_carlo.(Analysis.plural(value_field));
            time = monte_carlo.time_minutes;
            groups = coa.meas.measurementGroups(config.measurement_mask);

            fig = PlotOutput.newFigure(config, key, 'Name', title_text, 'Color', ...
                                       'w', ...
                                       'Position', [60, 60, 1100, 800]);
            layout = tiledlayout(ceil(numel(groups) / 2), 2, ...
                                 'TileSpacing', 'compact', 'Padding', 'compact');
            for g = 1:numel(groups)
                ax = nexttile(layout);
                hold(ax, 'on');
                rows = groups(g).rows;
                colors = lines(numel(rows));
                for r = 1:numel(rows)
                    normalized = squeeze(samples(:, :, rows(r))) ./ sigma(rows(r), :);
                    plot(ax, repmat(time, size(normalized, 1), 1), normalized, '.', ...
                         'Color', [0.7, 0.7, 0.7], 'MarkerSize', 4);
                end
                yline(ax, [-3, 3], 'k--');
                for r = 1:numel(rows)
                    plot(ax, time, reference(rows(r), :), 'o', ...
                         'Color', colors(r, :), 'LineWidth', 1, 'MarkerSize', 3);
                end
                grid(ax, 'on');
                title(ax, sprintf('%s (%d rows)', groups(g).name, numel(rows)), ...
                      'FontSize', 9);
                ylim(ax, [-5, 5]);
            end
            xlabel(layout, 'Time [min]');
            ylabel(layout, 'Value / predicted \sigma');
            title(layout, ...
                  sprintf(['%s: reference run (colored, one per row) over %d ' ...
                           'Monte Carlo runs (gray); dashed: \\pm3\\sigma'], ...
                          title_text, monte_carlo.run_count));
            PlotOutput.finishFigure(fig, config, key);
        end

        function name = plural(field)
            % Name of the Monte Carlo array holding a per-update field.
            name = [field, 's'];
        end

        function plotCorner(monte_carlo, config)
            % 18-D corner plot of the full-state errors at each final-time endpoint.
            import tetherFilterDemoSupport.*
            if PlotOutput.plotMode(config, 'corner') == "off"
                return
            end
            colors = lines(3);
            colors = colors(arrayfun(@Analysis.bodyOf, 1:18), :);
            for endpoint = monte_carlo.endpoints
                Analysis.drawCornerPlot(config, 'corner', ...
                                        "Error corner plot, " + endpoint.label, ...
                                        endpoint.errors, ...
                                        endpoint.covariance, ...
                                        Analysis.stateLabels(), colors, ...
                                        sprintf('Error corner plot at %s, %d runs', ...
                                                endpoint.label, monte_carlo.run_count));
            end
        end

        function drawCornerPlot(config, key, name, errors, covariance, labels, ...
                                colors, ...
                                title_text)
            % Generic corner plot of [runs x n] errors against an n x n covariance,
            % standardized by predicted sigma. Standardizing gives every axis the same
            % +/-limit range, so cells can be packed edge to edge with tick labels only
            % on the outer border. Labels and colors are per state. The figure is shown
            % or saved according to config.plot_modes.(key).
            %   Diagonal: marginal histogram vs. the predicted N(0,1).
            %   Lower:    samples with the predicted covariance ellipse.
            %   Upper:    sample minus predicted correlation (blue/red = under/over).
            import tetherFilterDemoSupport.*
            radius = config.corner_ellipse_sigma;
            sigma = sqrt(diag(covariance));
            z = errors ./ sigma';
            predicted_corr = covariance ./ (sigma * sigma');
            sample_corr = corrcoef(z);
            n = size(z, 2);
            font_size = max(6, min(10, 60 / n));
            marker_size = max(3, min(10, 54 / n));

            limit = 4;
            corr_limit = 0.5;
            cmap = [linspace(0, 1, 32)', linspace(0, 1, 32)', ones(32, 1); ...
                    ones(32, 1), linspace(1, 0, 32)', linspace(1, 0, 32)'];
            angle = linspace(0, 2 * pi, 100);
            unit_circle = [cos(angle); sin(angle)];
            edges = linspace(-limit, limit, 17);
            grid_points = linspace(-limit, limit, 200);
            bar_x = [edges(1:end - 1); edges(2:end); edges(2:end); edges(1:end - 1)];

            screen = get(0, 'ScreenSize');
            side = 0.9 * min(screen(3:4));
            fig = PlotOutput.newFigure(config, key, 'Name', name, 'Color', 'w', ...
                                       'Position', [screen(1) + 20, ...
                                                    screen(2) + 40, side, side]);
            extra = 0.1 * (size(labels, 1) > 1);
            margin = struct('left', 0.06 + extra, 'right', 0.10, ...
                            'bottom', 0.06 + 0.5 * extra, 'top', 0.09);
            cell_w = (1 - margin.left - margin.right) / n;
            cell_h = (1 - margin.bottom - margin.top) / n;

            for i = 1:n
                for j = 1:n
                    position = [margin.left + (j - 1) * cell_w, ...
                                margin.bottom + (n - i) * cell_h, cell_w, cell_h];
                    ax = axes(fig, 'Position', position, 'Box', 'on', ...
                              'FontSize', font_size, 'TickLength', [0.03, 0.03], ...
                              'XLim', [-limit, limit], 'YLim', [-limit, limit], ...
                              'XTick', [-2, 0, 2], 'YTick', [-2, 0, 2], ...
                              'NextPlot', 'add'); %#ok<LAXES>

                    if i == j
                        height = histcounts(z(:, j), edges, 'Normalization', 'pdf');
                        patch(ax, bar_x, [0; 0; 1; 1] * height, colors(j, :), ...
                              'EdgeColor', 'none', 'FaceAlpha', 0.8);
                        plot(ax, grid_points, ...
                             exp(-grid_points.^2 / 2) / sqrt(2 * pi), ...
                             'k', 'LineWidth', 0.8);
                        set(ax, 'YLim', [0, max(0.6, 1.05 * max(height))], 'YTick', []);
                    elseif i > j
                        scatter(ax, z(:, j), z(:, i), marker_size, colors(i, :), ...
                                'filled', 'MarkerFaceAlpha', 0.5);
                        pair = [j, i];
                        ellipse = radius * chol(predicted_corr(pair, pair), ...
                                                'lower') * ...
                                  unit_circle;
                        plot(ax, ellipse(1, :), ellipse(2, :), 'k-', 'LineWidth', 0.7);
                    else
                        difference = sample_corr(i, j) - predicted_corr(i, j);
                        level = (min(max(difference, -corr_limit), corr_limit) + ...
                                 corr_limit) / (2 * corr_limit);
                        color_index = 1 + round(level * (size(cmap, 1) - 1));
                        set(ax, 'Color', cmap(color_index, :), ...
                            'XTick', [], 'YTick', []);
                        text(ax, 0, 0, sprintf('%+.2f', difference), ...
                             'FontSize', font_size, 'HorizontalAlignment', 'center');
                    end

                    if i < n
                        set(ax, 'XTickLabel', []);
                    end
                    if j > 1 || i == 1
                        set(ax, 'YTickLabel', []);
                    end
                    if i == n
                        xlabel(ax, labels(:, j), 'FontSize', font_size + 1, ...
                               'Color', colors(j, :));
                    end
                    if j == 1
                        ylabel(ax, labels(:, i), 'FontSize', font_size + 1, ...
                               'Rotation', 0, ...
                               'HorizontalAlignment', 'right', 'Color', colors(i, :));
                    end
                end
            end

            colormap(fig, cmap);
            colorbar_ax = axes(fig, 'Position', [0, 0, 1e-3, 1e-3], 'Visible', 'off');
            caxis(colorbar_ax, [-corr_limit, corr_limit]);
            bar = colorbar(colorbar_ax, 'Position', ...
                           [1 - margin.right + 0.03, margin.bottom + 0.3, 0.012, 0.3]);
            bar.Label.String = 'Sample - predicted correlation (upper triangle)';
            bar.FontSize = 7;

            sgtitle(fig, {sprintf('%s, axes in predicted \\sigma', title_text), ...
                          sprintf(['lower: samples with predicted %g\\sigma ' ...
                                   'ellipse | diagonal: histogram vs. N(0,1) | ' ...
                                   'upper: correlation mismatch'], radius)}, ...
                    'FontSize', 10);
            PlotOutput.finishFigure(fig, config, key);
        end

        function relative = relativeStateDistribution(endpoint, config)
            % Errors and covariance of d = [r3 - r2; v3 - v2], expressed in the tether
            % frame of the reference run's estimate at the endpoint, in meters and mm/s.
            % d is linear in the state, so errors and covariance map exactly:
            % e_d = M e, P_d = M P M'. No new samples are needed.
            state = endpoint.reference_state;
            frame = coa.orbit.tetherFrame(state(4:6) - state(1:3), ...
                                          state(13:15) - state(10:12));

            difference = zeros(6, 18);
            difference(1:3, 4:6) = -eye(3);
            difference(1:3, 7:9) = eye(3);
            difference(4:6, 13:15) = -eye(3);
            difference(4:6, 16:18) = eye(3);

            units = [config.length_scale * 1e3 * ones(3, 1); ...
                     config.length_scale / config.time_scale * 1e6 * ones(3, 1)];
            map = diag(units) * blkdiag(frame, frame) * difference;
            covariance = map * endpoint.covariance * map';

            relative.frame = frame;
            relative.errors = endpoint.errors * map';
            relative.covariance = (covariance + covariance') / 2;
            relative.names = ["\Delta r_{e1}", "\Delta r_{e2}", "\Delta r_{e3}", ...
                              "\Delta v_{e1}", "\Delta v_{e2}", "\Delta v_{e3}"];
            relative.units = ["m", "m", "m", "mm/s", "mm/s", "mm/s"];
        end

        function reportRelativeState(monte_carlo, config)
            % Predicted vs. sample spread of the relative state in the tether frame.
            import tetherFilterDemoSupport.*
            for endpoint = monte_carlo.endpoints
                relative = Analysis.relativeStateDistribution(endpoint, config);
                predicted = sqrt(diag(relative.covariance));
                sample = sqrt(mean(relative.errors.^2, 1))';
                fprintf('\nRelative state r3-r2, v3-v2 in tether frame at %s\n', ...
                        endpoint.label);
                fprintf('%-8s %-6s  %12s  %12s  %s\n', 'Comp.', 'Unit', ...
                        'Predicted sd', ...
                        'Sample RMS', 'Ratio');
                for k = 1:6
                    fprintf('%-8s %-6s  %12.4g  %12.4g  %5.2f\n', ...
                            char(erase(relative.names(k), {'\', '{', '}'})), ...
                            relative.units(k), predicted(k), sample(k), ...
                            sample(k) / predicted(k));
                end
            end
        end

        function plotRelativeCorner(monte_carlo, config)
            % 6-D corner plot of the relative state in the tether frame, per endpoint.
            import tetherFilterDemoSupport.*
            if PlotOutput.plotMode(config, 'relative_corner') == "off"
                return
            end
            colors = [repmat([0, 0.447, 0.741], 3, 1); ...
                      repmat([0.85, 0.325, 0.098], 3, 1)];
            for endpoint = monte_carlo.endpoints
                relative = Analysis.relativeStateDistribution(endpoint, config);
                units = diag([1, 1, 1, 1e-3, 1e-3, 1e-3]);
                relative.errors = relative.errors * units;
                relative.covariance = units * relative.covariance * units;
                relative.units(4:6) = "m/s";
                sigma = sqrt(diag(relative.covariance));
                labels = strings(2, 6);
                for k = 1:6
                    labels(1, k) = relative.names(k);
                    labels(2, k) = sprintf('\\sigma=%.2g %s', sigma(k), ...
                                           relative.units(k));
                end
                name = "Relative state corner plot (tether frame), " + endpoint.label;
                plot_title = sprintf(['Error of [r_3-r_2; v_3-v_2] in the ' ...
                                      'tether frame at %s, %d runs'], ...
                                     endpoint.label, monte_carlo.run_count);
                Analysis.drawCornerPlot(config, 'relative_corner', name, ...
                                        relative.errors, relative.covariance, ...
                                        labels, colors, plot_title);
            end
        end

        function body = bodyOf(index)
            % Body (1-3) owning a state index; positions then velocities, 3 per body.
            body = mod(floor((index - 1) / 3), 3) + 1;
        end

        function labels = stateLabels()
            % Short labels p<body><axis> and v<body><axis> in state-vector order.
            import tetherFilterDemoSupport.*
            axes_names = ["x", "y", "z"];
            labels = strings(1, 18);
            for index = 1:18
                kind = "p";
                if index > 9
                    kind = "v";
                end
                labels(index) = kind + Analysis.bodyOf(index) + ...
                                axes_names(mod(index - 1, 3) + 1);
            end
        end

    end
end
