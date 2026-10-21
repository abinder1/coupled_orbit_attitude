% Copyright (c) 2026 Andrew Binder

% tetherHaloDemo
%   Demonstrates coa.orbit.tetherAccel by propagating a spinning two-mass tether
%   that is carried along a CR3BP halo orbit, with its length prescribed as
%   L(t) = L0 + Ldot * t. Three cases are compared: reeling in (Ldot = -0.5 m/s),
%   a fixed-length tether (Ldot = 0), and reeling out (Ldot = +0.5 m/s).
%
%   Each mass moves under the unconstrained CR3BP acceleration plus the
%   constraint force that keeps their separation on the prescribed schedule.
%   tetherAccel returns those constrained accelerations directly, given a handle
%   to the single-particle CR3BP physics.
%
%   Figures:
%     1. Tether length versus time, with the constraint residual.
%     2. Tether geometry relative to its center of mass, colored by time.
%     3. Context: the tether's center of mass along the halo orbit.
%
%   State convention (nondimensional, CR3BP rotating frame):
%     sv = [p_cw; p_cm; v_cw; v_cm], with CW = counterweight, CM = catch mechanism.

%% Script setup
clear;
clc;
close('all');
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'toolbox'));

save_figures = false;  % Set true to export each figure as a vector PDF

ode_opts = odeset('AbsTol', 1e-10, 'RelTol', 1e-10);

%% System parameters
% Earth-Moon CR3BP characteristic scales
L_star = 384400;            % Length unit [km]
T_star = 375190.258892627;  % Time unit [s]
mu = 0.012150585;           % Moon's mass fraction

% Tether parameterization
L = 200;                    % Initial tether length [km]
spin_rate = 2 * pi / (55 * 60);  % Inertial spin rate about +z [rad/s]
eta = 0.5;                  % Catch mechanism's mass fraction, m_cm / (m_cw + m_cm)
lhat_0 = [1; 0; 0];         % Initial unit vector from the CM to the CW

% Halo orbit that the tether's center of mass starts on (rotating frame)
halo_IC = [0.831987145310294; 0; 0.12575552414079; 0; 0.239887839725451; 0];
halo_PD = 2.78261469021996;  % Halo period (nondimensional)

% Propagation span and the reel-rate cases to compare. Over two days, +-0.5 m/s
% changes the length by +-86 km, which keeps the 200 km tether well away from
% zero length (where the separation constraint is singular).
prop_days = 2;
rates_ms = [-0.5, 0, 0.5];  % Prescribed Ldot [m/s]
n_samples = 3000;

%% Case definitions
% Colorblind-safe palette (Okabe-Ito)
colors = [0.835, 0.369, 0.000; ...  % vermillion: reeling in
          0.000, 0.447, 0.698; ...  % blue: fixed length
          0.000, 0.619, 0.451];     % green: reeling out
labels = ["Reel in, $\dot{L} = -0.5$ m/s", ...
          "Fixed length, $\dot{L} = 0$", ...
          "Reel out, $\dot{L} = +0.5$ m/s"];

% The physics of a single unconstrained particle. tetherAccel calls this once
% per mass, so any other force model with the same signature can be swapped in.
accel_fn = @(x) coa.orbit.cr3bp(x, mu);

t_eval = linspace(0, prop_days * 86400 / T_star, n_samples);

cases = struct('rate_ms', num2cell(rates_ms));

%% Propagate each case
for k = 1:numel(cases)
    length_dot = cases(k).rate_ms / 1000 * T_star / L_star;  % nondimensional

    sv_0 = tetherInitialState(halo_IC, L / L_star, eta, ...
                              spin_rate * T_star, length_dot, lhat_0);

    ode_func = @(~, sv) tetherDynamics(sv, eta, length_dot, accel_fn);
    sol = ode113(ode_func, [t_eval(1), t_eval(end)], sv_0, ode_opts);
    sv = deval(sol, t_eval);

    cw = sv(1:3, :) * L_star;
    cm = sv(4:6, :) * L_star;

    cases(k).label = labels(k);
    cases(k).color = colors(k, :);
    cases(k).t_days = t_eval * T_star / 86400;
    cases(k).cw = cw;
    cases(k).cm = cm;
    cases(k).com = (1 - eta) * cw + eta * cm;
    cases(k).length = vecnorm(cm - cw);

    % The constraint should hold: |p_cm - p_cw| = L0 + Ldot * t
    prescribed = L + cases(k).rate_ms / 1000 * t_eval * T_star;
    cases(k).prescribed = prescribed;
    cases(k).residual_m = abs(cases(k).length - prescribed) * 1000;
end

%% Reference halo orbit for context (single particle, no tether)
halo_func = @(~, x) [x(4:6); coa.orbit.cr3bp(x, mu)];
halo_sol = ode113(halo_func, [0, halo_PD], halo_IC, ode_opts);
halo_xyz = deval(halo_sol, linspace(0, halo_PD, 600)) * L_star;
halo_xyz = halo_xyz(1:3, :);

%% Figures
figs = [plotLengthHistory(cases), ...
        plotTetherGeometry(cases, eta, prop_days), ...
        plotHaloContext(cases, halo_xyz, mu, L_star), ...
        plotEndpointMotion(cases(rates_ms == 0), 1)];

if save_figures
    out_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');
    exportFigures(figs, out_dir, ["tetherLength", "tetherGeometry", ...
                                  "tetherHalo", "tetherEndpoints"]);
end

%% Local functions: dynamics and initial conditions
function sv_0 = tetherInitialState(com_state, length_nd, eta, spin_rate_nd, ...
                                   length_dot, lhat_0)
    % Initial state [p_cw; p_cm; v_cw; v_cm] for a tether whose center of mass
    % starts at com_state = [p; v] in the rotating frame.
    omega_frame = [0; 0; 1];            % Rotating-frame angular velocity
    spin = spin_rate_nd * [0; 0; 1];    % Tether's inertial angular velocity

    % The CW is a distance eta * L from the center of mass on the +lhat side and
    % the CM is (1 - eta) * L on the other side.
    rel_cw = lhat_0 * (eta * length_nd);
    rel_cm = -lhat_0 * ((1 - eta) * length_nd);

    % Rigid spin as seen in the rotating frame, plus the extension velocity of
    % each mass. The extension terms leave the center of mass unaffected.
    spin_rel = spin - omega_frame;
    vel_cw = cross(spin_rel, rel_cw) + eta * length_dot * lhat_0;
    vel_cm = cross(spin_rel, rel_cm) - (1 - eta) * length_dot * lhat_0;

    com_pos = com_state(1:3);
    com_vel = com_state(4:6);

    sv_0 = [com_pos + rel_cw; com_pos + rel_cm; com_vel + vel_cw; com_vel + vel_cm];
end

function sv_dot = tetherDynamics(sv, eta, length_dot, accel_fn)
    % Time derivative of sv = [p_cw; p_cm; v_cw; v_cm]. tetherAccel returns the
    % stacked constrained accelerations [a_cw; a_cm].
    x_cw = sv([1:3, 7:9]);
    x_cm = sv([4:6, 10:12]);

    a = coa.orbit.tetherAccel(x_cw, x_cm, eta, length_dot, accel_fn);

    sv_dot = [sv(7:12); a];
end

%% Local functions: plotting
function fig = plotLengthHistory(cases)
    % Tether length versus time, and how well the constraint is satisfied.
    fig = figure('Units', 'inches', 'Position', [1, 1, 6.5, 5.2], 'Color', 'w');
    tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax1 = nexttile(tl);
    hold(ax1, 'on');
    handles = gobjects(1, numel(cases) + 1);
    for k = 1:numel(cases)
        h_ref = plot(ax1, cases(k).t_days, cases(k).prescribed, '-', ...
                     'Color', [0.8, 0.8, 0.8], 'LineWidth', 4);
        if k == 1
            handles(end) = h_ref;
        end
    end
    for k = 1:numel(cases)
        handles(k) = plot(ax1, cases(k).t_days, cases(k).length, '-', ...
                          'Color', cases(k).color, 'LineWidth', 1.8);
    end
    styleAxes(ax1);
    ylabel(ax1, 'Tether length [km]', 'Interpreter', 'latex');
    title(ax1, 'Tether Length Follows the Prescribed Schedule', ...
          'Interpreter', 'latex');
    legend(ax1, handles, [cases.label, "Prescribed $L_0 + \dot{L}\,t$"], ...
           'Interpreter', 'latex', 'Location', 'northwest', 'Box', 'off');
    ax1.XTickLabel = [];

    ax2 = nexttile(tl);
    hold(ax2, 'on');
    for k = 1:numel(cases)
        % Skip t = 0, where the residual is exactly zero
        semilogy(ax2, cases(k).t_days(2:end), cases(k).residual_m(2:end), '-', ...
                 'Color', cases(k).color, 'LineWidth', 1.4);
    end
    ax2.YScale = 'log';
    styleAxes(ax2);
    xlabel(ax2, 'Elapsed time [days]', 'Interpreter', 'latex');
    ylabel(ax2, 'Constraint residual [m]', 'Interpreter', 'latex');
    title(ax2, ['Constraint Residual ', ...
                '$\left| \, \|\mathbf{d}\|^2 - (L_0 + \dot{L}t)^2 \right|$'], ...
          'Interpreter', 'latex');
end

function fig = plotTetherGeometry(cases, eta, prop_days)
    % Tether endpoints relative to the tether's center of mass, in the rotating
    % frame's x-y plane. Snapshots of the tether are colored by elapsed time.
    n_snapshots = 24;
    cmap = parula(256);

    fig = figure('Units', 'inches', 'Position', [1, 1, 7.0, 3.3], 'Color', 'w');
    tl = tiledlayout(fig, 1, numel(cases), 'TileSpacing', 'compact', ...
                     'Padding', 'compact');

    % Common axis limits so the three panels can be compared directly
    lim = 0;
    for k = 1:numel(cases)
        rel = [cases(k).cw - cases(k).com, cases(k).cm - cases(k).com];
        lim = max(lim, max(abs(rel(1:2, :)), [], 'all'));
    end
    lim = 1.08 * lim;

    handles = gobjects(1, 4);
    for k = 1:numel(cases)
        ax = nexttile(tl);
        hold(ax, 'on');

        rel_cw = cases(k).cw - cases(k).com;
        rel_cm = cases(k).cm - cases(k).com;
        n = size(rel_cw, 2);

        % Tether snapshots, colored by time
        for idx = round(linspace(1, n, n_snapshots))
            color = cmap(max(1, ceil(cases(k).t_days(idx) / prop_days * 256)), :);
            plot(ax, [rel_cw(1, idx), rel_cm(1, idx)], ...
                 [rel_cw(2, idx), rel_cm(2, idx)], '-', ...
                 'Color', [color, 0.8], 'LineWidth', 0.6);
        end

        % Endpoint paths and the t = 0 positions
        h_cw = plot(ax, rel_cw(1, :), rel_cw(2, :), '-', ...
                    'Color', [0.15, 0.15, 0.15], 'LineWidth', 0.8);
        h_cm = plot(ax, rel_cm(1, :), rel_cm(2, :), '-', ...
                    'Color', [0.55, 0.55, 0.55], 'LineWidth', 0.8);
        h_cw0 = plot(ax, rel_cw(1, 1), rel_cw(2, 1), 's', 'MarkerSize', 5, ...
                     'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
        h_cm0 = plot(ax, rel_cm(1, 1), rel_cm(2, 1), 'o', 'MarkerSize', 5, ...
                     'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
        handles = [h_cw, h_cm, h_cw0, h_cm0];

        axis(ax, 'equal');
        xlim(ax, [-lim, lim]);
        ylim(ax, [-lim, lim]);
        colormap(ax, cmap);
        clim(ax, [0, prop_days]);
        styleAxes(ax);
        xlabel(ax, '$x - x_{\mathrm{COM}}$ [km]', 'Interpreter', 'latex');
        if k == 1
            ylabel(ax, '$y - y_{\mathrm{COM}}$ [km]', 'Interpreter', 'latex');
        else
            ax.YTickLabel = [];
        end
        title(ax, cases(k).label, 'Interpreter', 'latex', 'FontSize', 9);
    end

    title(tl, sprintf(['Tether Motion About Its Center of Mass ', ...
                       '($\\eta = %.1f$, rotating frame)'], eta), ...
          'Interpreter', 'latex');

    cb = colorbar(ax);
    cb.Layout.Tile = 'east';
    cb.Label.String = 'Elapsed time [days]';
    cb.Label.Interpreter = 'latex';
    cb.TickLabelInterpreter = 'latex';

    lg = legend(ax, handles, {'CW path', 'CM path', 'CW at $t=0$', 'CM at $t=0$'}, ...
                'Interpreter', 'latex', 'Orientation', 'horizontal', 'Box', 'off');
    lg.Layout.Tile = 'south';
end

function fig = plotHaloContext(cases, halo_xyz, mu, l_star)
    % Where the tether is, relative to the halo orbit and the Moon.
    scale = 1e3;  % Plot in thousands of km

    fig = figure('Units', 'inches', 'Position', [1, 1, 5.5, 4.6], 'Color', 'w');
    ax = axes(fig);
    hold(ax, 'on');

    plot3(ax, halo_xyz(1, :) / scale, halo_xyz(2, :) / scale, ...
          halo_xyz(3, :) / scale, '-', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.2, ...
          'DisplayName', 'Halo orbit');

    % The center-of-mass paths nearly coincide; thinner lines drawn on top keep
    % all three visible.
    widths = linspace(4, 1.2, numel(cases));
    for k = 1:numel(cases)
        com = cases(k).com / scale;
        plot3(ax, com(1, :), com(2, :), com(3, :), '-', ...
              'Color', cases(k).color, 'LineWidth', widths(k), ...
              'DisplayName', cases(k).label);
    end

    com_0 = cases(1).com(:, 1) / scale;
    scatter3(ax, com_0(1), com_0(2), com_0(3), 40, 'k', 'filled', 'o', ...
             'DisplayName', 'Start');
    scatter3(ax, (1 - mu) * l_star / scale, 0, 0, 160, [0.45, 0.45, 0.45], ...
             'filled', 'DisplayName', 'Moon');

    axis(ax, 'equal');
    view(ax, 3);
    styleAxes(ax);
    xlabel(ax, '$x$ [$10^3$ km]', 'Interpreter', 'latex');
    ylabel(ax, '$y$ [$10^3$ km]', 'Interpreter', 'latex');
    zlabel(ax, '$z$ [$10^3$ km]', 'Interpreter', 'latex');
    title(ax, 'Tether Center of Mass Along the Halo Orbit', 'Interpreter', 'latex');
    legend(ax, 'Interpreter', 'latex', 'Location', 'northeast', 'Box', 'off');
end

function fig = plotEndpointMotion(tether, span_days)
    % Both endpoints of the fixed-length tether in the CR3BP rotating frame (x-y
    % projection) over the first span_days. At full scale the 200 km tether is
    % far smaller than the path, so the left panel shows the whole arc and the
    % right panel zooms in on the highlighted window, where the spin of the two
    % masses about their center of mass becomes visible.
    scale = 1e3;               % Left panel in thousands of km
    zoom_hours = 2;            % Width of the zoomed time window
    n_snapshots = 7;
    cmap = jet(256);

    keep = tether.t_days <= span_days;
    t = tether.t_days(keep);
    cw = tether.cw(:, keep);
    cm = tether.cm(:, keep);
    in_zoom = t <= zoom_hours / 24;

    fig = figure('Units', 'inches', 'Position', [1, 1, 7.0, 3.5], 'Color', 'w');
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Full-arc view
    ax1 = nexttile(tl);
    hold(ax1, 'on');
    h_cw = plot(ax1, cw(1, :) / scale, cw(2, :) / scale, '-', ...
                'Color', [0.835, 0.369, 0.000], 'LineWidth', 1.6);
    h_cm = plot(ax1, cm(1, :) / scale, cm(2, :) / scale, '--', ...
                'Color', [0.000, 0.447, 0.698], 'LineWidth', 1.6);
    h_0 = plot(ax1, cw(1, 1) / scale, cw(2, 1) / scale, 'o', 'MarkerSize', 5, ...
               'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    h_f = plot(ax1, cw(1, end) / scale, cw(2, end) / scale, 's', ...
               'MarkerSize', 5, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');

    % Translucent highlight marking the zoomed region; the zoom panel's frame
    % uses the same accent color so the two read as linked.
    accent = [0.000, 0.350, 0.150];  % Dark green
    both = [cw(1:2, in_zoom), cm(1:2, in_zoom)] / scale;
    lo = min(both, [], 2);
    hi = max(both, [], 2);
    pad = 0.25 * max(hi - lo);
    box_x = [lo(1) - pad, hi(1) + pad, hi(1) + pad, lo(1) - pad];
    box_y = [lo(2) - pad, lo(2) - pad, hi(2) + pad, hi(2) + pad];
    h_box = patch(ax1, box_x, box_y, accent, 'FaceAlpha', 0.25, ...
                  'EdgeColor', accent, 'LineWidth', 1.8);
    axis(ax1, 'equal');
    styleAxes(ax1);
    xlabel(ax1, '$x$ [$10^3$ km]', 'Interpreter', 'latex');
    ylabel(ax1, '$y$ [$10^3$ km]', 'Interpreter', 'latex');
    title(ax1, sprintf('Full %g-Day Arc', span_days), 'Interpreter', 'latex');
    legend(ax1, [h_cw, h_cm, h_0, h_f, h_box], ...
           {'CW', 'CM', '$t = 0$', sprintf('$t = %g$ d', span_days)}, ...
           'Interpreter', 'latex', 'Location', 'northwest', 'Box', 'off');

    % Zoomed view with tether snapshots colored by time
    ax2 = nexttile(tl);
    hold(ax2, 'on');
    idx_zoom = find(in_zoom);
    origin = cw(1:2, 1);  % Offsets keep the tick labels short
    cw_z = cw(1:2, :) - origin;
    cm_z = cm(1:2, :) - origin;
    for idx = idx_zoom(round(linspace(1, numel(idx_zoom), n_snapshots)))
        color = cmap(max(1, ceil(t(idx) / (zoom_hours / 24) * 256)), :);
        plot(ax2, [cw_z(1, idx), cm_z(1, idx)], [cw_z(2, idx), cm_z(2, idx)], ...
             '-', 'Color', [color, 0.8], 'LineWidth', 1.5);
    end
    plot(ax2, cw_z(1, in_zoom), cw_z(2, in_zoom), '-', ...
         'Color', [0.835, 0.369, 0.000], 'LineWidth', 1.6);
    plot(ax2, cm_z(1, in_zoom), cm_z(2, in_zoom), '-', ...
         'Color', [0.000, 0.447, 0.698], 'LineWidth', 1.6);
    plot(ax2, 0, 0, 'o', 'MarkerSize', 5, ...
         'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    axis(ax2, 'equal');
    % Explicit limits leave a small margin of whitespace around the motion
    zoom_pts = [cw_z(:, in_zoom), cm_z(:, in_zoom)];
    zoom_lo = min(zoom_pts, [], 2);
    zoom_hi = max(zoom_pts, [], 2);
    margin = 0.08 * max(zoom_hi - zoom_lo);
    xlim(ax2, [zoom_lo(1) - margin, zoom_hi(1) + margin]);
    ylim(ax2, [zoom_lo(2) - margin, zoom_hi(2) + margin]);
    colormap(ax2, cmap);
    clim(ax2, [0, zoom_hours]);
    styleAxes(ax2);
    % The window is tall and narrow, so thin the x ticks to keep labels legible
    x_step = 200;
    x_lim = xlim(ax2);
    xticks(ax2, x_step * (ceil(x_lim(1) / x_step):floor(x_lim(2) / x_step)));
    ax2.XMinorTick = 'off';
    xtickangle(ax2, 0);
    set(ax2, 'XColor', accent, 'YColor', accent, 'LineWidth', 1.8);
    xlabel(ax2, '$x - x_0$ [km]', 'Interpreter', 'latex');
    ylabel(ax2, '$y - y_0$ [km]', 'Interpreter', 'latex');
    title(ax2, sprintf('Zoom: First %g Hours', zoom_hours), ...
          'Interpreter', 'latex');
    cb = colorbar(ax2);
    cb.Label.String = 'Elapsed time [hours]';
    cb.Label.Interpreter = 'latex';
    cb.TickLabelInterpreter = 'latex';

    title(tl, 'Tether Endpoints in the CR3BP Rotating Frame', ...
          'Interpreter', 'latex');

    % Connector lines from the highlighted box to the zoom panel's corners. They
    % are placed in figure coordinates, so they track the layout at the size the
    % figure is drawn (re-run if the window is resized).
    drawnow;
    p1 = dataToFigure(ax1, [box_x(2), box_x(3)], [box_y(2), box_y(3)]);
    p2 = dataToFigure(ax2, [1, 1] * min(xlim(ax2)), ylim(ax2));
    for k = 1:2
        annotation(fig, 'line', [p1(1, k), p2(1, k)], [p1(2, k), p2(2, k)], ...
                   'Color', accent, 'LineWidth', 1.2);
    end
end

function pts = dataToFigure(ax, x, y)
    % Map data points to normalized figure coordinates, accounting for the
    % smaller plot box that "axis equal" leaves inside the axes position.
    fig = ancestor(ax, 'figure');
    old_units = fig.Units;
    fig.Units = 'pixels';
    fig_size = fig.Position(3:4);
    fig.Units = old_units;

    inner = ax.InnerPosition .* [fig_size, fig_size];
    x_lim = xlim(ax);
    y_lim = ylim(ax);
    data_aspect = diff(y_lim) / diff(x_lim);
    if inner(4) / inner(3) > data_aspect
        box_w = inner(3);
        box_h = box_w * data_aspect;
    else
        box_h = inner(4);
        box_w = box_h / data_aspect;
    end
    box_x0 = inner(1) + (inner(3) - box_w) / 2;
    box_y0 = inner(2) + (inner(4) - box_h) / 2;

    pts = [(box_x0 + (x - x_lim(1)) / diff(x_lim) * box_w) / fig_size(1); ...
           (box_y0 + (y - y_lim(1)) / diff(y_lim) * box_h) / fig_size(2)];
end

function styleAxes(ax)
    % Shared publication-style formatting.
    set(ax, 'FontSize', 10, 'LineWidth', 0.75, 'Box', 'on', 'TickDir', 'in', ...
        'TickLabelInterpreter', 'latex', 'Layer', 'top', ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'GridAlpha', 0.15, 'MinorGridAlpha', 0.08);
    grid(ax, 'on');
end

function exportFigures(figs, out_dir, names)
    % Save each figure as a vector PDF.
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end
    for k = 1:numel(figs)
        exportgraphics(figs(k), fullfile(out_dir, names(k) + ".pdf"), ...
                       'ContentType', 'vector');
    end
end
