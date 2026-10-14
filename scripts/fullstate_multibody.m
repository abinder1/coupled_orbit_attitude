%% Script setup
clear; clc; close('all');
addpath(genpath('./utilities'))

ode_opts = odeset('AbsTol', 1e-6, 'RelTol', 1e-6);

L_star = 384400;
T_star = 375190.258892627;

%% Start defining parameters of the system
% Tether parameterization
spin_rate = 2 * pi / (55 * 60);
L = 200;             % Length of the tether, kilometers
eta = 0.5;           % Tether mass balance, eta = catch mechanism / total
lhat_0 = [1; 0; 0];  % Initial pointing direction, CW -> CM
spin_0 = [0; 0; 1];  % Initial tether spin direction

% Halo orbit parameterization
mu = 0.012150585;
halo_IC = [0.831987145310294; 0; 0.12575552414079; 0; 0.239887839725451; 0];
halo_IC(4:6) = halo_IC(4:6) + cross([0; 0; 1], halo_IC(1:3));  % Inertial frame IC

halo_PD = 2.78261469021996;

%% Derive the payload and tether's initial condition dimensionally
xs = eta * L;

counterweight_relative_position = -lhat_0 * xs;
catchmech_relative_position = lhat_0 * (L - xs);

spin_av = spin_rate * spin_0;

counterweight_delta_state = [counterweight_relative_position / L_star; ...
                             cross(spin_av, counterweight_relative_position) * T_star / L_star];

catchmech_delta_state = [catchmech_relative_position / L_star; ...
                         cross(spin_av, catchmech_relative_position) * T_star / L_star];

cw_IC = halo_IC - counterweight_delta_state; 
cm_IC = halo_IC - catchmech_delta_state;
pl_IC = halo_IC; % cm_IC;

% Reorder IC vector to be [r1; r2; r3; v1; v2; v3]
total_IC = [cw_IC(1:3); cm_IC(1:3); pl_IC(1:3); cw_IC(4:6); cm_IC(4:6); pl_IC(4:6)];

Ldot = 0.00 * T_star / L_star;

total_IC(10:12) = total_IC(10:12) + eta * Ldot * lhat_0;
total_IC(13:15) = total_IC(13:15) - (1 - eta) * Ldot * lhat_0;

%% Check Jacobian of the dynamics with finite differencing
% Augment the initial conditions with an 18 x 18 state transition matrix
sv_0 = [total_IC; reshape(eye(18), [324, 1])];
prop_time = 0.1;
ode_func = @(tau, sv) constrained_gravity_STM(tau, sv, eta, mu);

% Propagate backwards first by five minutes to get the right IC at that time
sol_struct = ode89(ode_func, [0, prop_time], sv_0, ode_opts);

%% Propagate state and state transition matrix nondimensionally
STM_f = reshape(sol_struct.y(19:342, end), [18, 18]);

cw_pos = sol_struct.y(1:3, :); cm_pos = sol_struct.y(4:6, :);
pl_pos = sol_struct.y(7:9, :);

figure; hold on; grid on; axis equal;
plot3(cw_pos(1, :) * L_star, cw_pos(2, :) * L_star, cw_pos(3, :) * L_star, 'k')
plot3(cm_pos(1, :) * L_star, cm_pos(2, :) * L_star, cm_pos(3, :) * L_star, 'b')
plot3(pl_pos(1, :) * L_star, pl_pos(2, :) * L_star, pl_pos(3, :) * L_star, 'r')
scatter3(cw_pos(1, 1) * L_star, cw_pos(2, 1) * L_star, cw_pos(3, 1) * L_star, 'ks', 'filled')
scatter3(cm_pos(1, 1) * L_star, cm_pos(2, 1) * L_star, cm_pos(3, 1) * L_star, 'bs', 'filled')
scatter3(pl_pos(1, 1) * L_star, pl_pos(2, 1) * L_star, pl_pos(3, 1) * L_star, 'rs', 'filled')

title_text = compose("The Last %d Seconds Before Catch", prop_time * T_star);

title(title_text)

xlabel("X-position [km]")
ylabel("Y-position [km]")
zlabel("Z-position [km]")

legend("Counterweight", "Catch Mechanism", "Payload", ...
       "CW IC", "CM IC", "PL IC", "Location", "best")

figure; hold on; grid on;
plot(sol_struct.x * T_star, vecnorm(cw_pos - cm_pos) * L_star)

%% Visualize the payload trajectory in the rotating frame
cw_rot = zeros(size(pl_pos));
cm_rot = zeros(size(pl_pos));
pl_rot = zeros(size(pl_pos));

for k = 1:1:length(sol_struct.x)
    tau = sol_struct.x(k);

    DCM = [cos(tau), sin(tau), 0; -sin(tau), cos(tau), 0; 0, 0, 1];

    cw_rot(:, k) = DCM * cw_pos(:, k);
    cm_rot(:, k) = DCM * cm_pos(:, k);
    pl_rot(:, k) = DCM * pl_pos(:, k);
end

figure; hold on; grid on; axis equal;
plot3(cw_rot(1, :) * L_star, cw_rot(2, :) * L_star, cw_rot(3, :) * L_star, 'k')
plot3(cm_rot(1, :) * L_star, cm_rot(2, :) * L_star, cm_rot(3, :) * L_star, 'b')
plot3(pl_rot(1, :) * L_star, pl_rot(2, :) * L_star, pl_rot(3, :) * L_star, 'r')
scatter3(cw_rot(1, 1) * L_star, cw_rot(2, 1) * L_star, cw_rot(3, 1) * L_star, 'ks', 'filled')
scatter3(cm_rot(1, 1) * L_star, cm_rot(2, 1) * L_star, cm_rot(3, 1) * L_star, 'bs', 'filled')
scatter3(pl_rot(1, 1) * L_star, pl_rot(2, 1) * L_star, pl_rot(3, 1) * L_star, 'rs', 'filled')

title_text = compose("The Last %d Seconds Before Catch", prop_time * T_star);

title(title_text)

xlabel("X-position [km]")
ylabel("Y-position [km]")
zlabel("Z-position [km]")

legend("Counterweight", "Catch Mechanism", "Payload", ...
       "CW IC", "CM IC", "PL IC", "Location", "best")

%% Check the STM against finite differentiating
ode_func = @(tau, sv) constrained_gravity(tau, sv, eta, mu);
fd_stepsize = 1e-7;

STM_fd = zeros(size(STM_f));

for k = 1:1:18
    k
    d_IC = zeros(size(total_IC));  d_IC(k) = fd_stepsize;

    sol_struct = ode89(ode_func, [0, prop_time], total_IC + d_IC, ode_opts);
    xf_plus = sol_struct.y(:, end);

    sol_struct = ode89(ode_func, [0, prop_time], total_IC - d_IC, ode_opts);
    xf_minus = sol_struct.y(:, end);

    STM_fd(:, k) = (xf_plus - xf_minus) / 2 / fd_stepsize;
end

%% Local function definition
function sv_dot = constrained_gravity_STM(tau, sv, eta, mu)
    cw_pos = sv(1:3); cm_pos = sv(4:6); pl_pos = sv(7:9); 
    cw_vel = sv(10:12); cm_vel = sv(13:15); pl_vel = sv(16:18);

    stm = reshape(sv(19:342), [18, 18]);

    % Calculate the virtual accelerations of the two point masses
    [cw_virt_accel, cw_A] = virtual_acceleration(tau, cw_pos, mu);
    [cm_virt_accel, cm_A] = virtual_acceleration(tau, cm_pos, mu);
    [pl_virt_accel, pl_A] = virtual_acceleration(tau, pl_pos, mu);

    % Construct relative quantities
    dv = cm_pos - cw_pos;
    dvd = cm_vel - cw_vel;

    dh = dv / norm(dv);

    D = dh * dh';
    qv = dot(dvd, dvd) / dot(dv, dv) * dv;

    % Construct constrained accelerations
    O3 = zeros(3);  I3 = eye(3);

    accel_scaling = [I3 - eta * D, eta * D, O3; ...
                     (1 - eta) * D, I3 - (1 - eta)*D, O3; ...
                     O3, O3, I3];

    combined_qv = [eta * qv; -(1 - eta) * qv; 0; 0; 0];

    true_accels = accel_scaling * [cw_virt_accel; cm_virt_accel; pl_virt_accel] + combined_qv;
    
    % Construct the system Jacobian
    Z1 = dot(dh, cw_virt_accel) * I3 + (dh * cw_virt_accel');
    Z2 = dot(dh, cm_virt_accel) * I3 + (dh * cm_virt_accel');
    Z3 = (I3 - D) / norm(dv);
    Z4 = (dot(dvd, dvd) / dot(dv, dv)) * (I3 - 2 * D);
    Z5 = -2 * (dv * dvd') / dot(dv, dv);

    Z_combo = (Z1 - Z2) * Z3 - Z4;

    da1_dr1 = (I3 - eta * D) * cw_A + eta * Z_combo;
    da1_dv1 = -eta * Z5;
    da1_dr2 = eta * D * cm_A - eta * Z_combo;

    da2_dr1 = (1 - eta) * D * cw_A - (1 - eta) * Z_combo;
    da2_dv1 = -(1 - eta) * Z5;
    da2_dr2 = (I3 - (1 - eta) * D) * cm_A + (1 - eta) * Z_combo;

    acceleration_jacobian = [da1_dr1, da1_dr2, O3, -da1_dv1, da1_dv1, O3; ...
                             da2_dr1, da2_dr2, O3, da2_dv1, -da2_dv1, O3; ...
                             O3, O3, pl_A, zeros([3, 9])];

    full_jacobian = [zeros(9), eye(9); acceleration_jacobian];

    % Return sv_dot for ODE
    stm_dot = reshape(full_jacobian * stm, [324, 1]);
    sv_dot = [cw_vel; cm_vel; pl_vel; true_accels; stm_dot];
end


function sv_dot = constrained_gravity(tau, sv, eta, mu)
    cw_pos = sv(1:3); cm_pos = sv(4:6); pl_pos = sv(7:9); 
    cw_vel = sv(10:12); cm_vel = sv(13:15); pl_vel = sv(16:18);
        
    % Calculate the virtual accelerations of the two point masses
    [cw_virt_accel, ~] = virtual_acceleration(tau, cw_pos, mu);
    [cm_virt_accel, ~] = virtual_acceleration(tau, cm_pos, mu);
    [pl_virt_accel, ~] = virtual_acceleration(tau, pl_pos, mu);
    
    % Construct relative quantities
    dv = cm_pos - cw_pos;
    dvd = cm_vel - cw_vel;
    
    dh = dv / norm(dv);
    
    D = dh * dh';
    qv = dot(dvd, dvd) / dot(dv, dv) * dv;
    
    % Construct constrained accelerations
    O3 = zeros(3);  I3 = eye(3);
    
    accel_scaling = [I3 - eta * D, eta * D, O3; ...
        (1 - eta) * D, I3 - (1 - eta)*D, O3; ...
        O3, O3, I3];
    
    combined_qv = [eta * qv; -(1 - eta) * qv; 0; 0; 0];
    
    true_accels = accel_scaling * [cw_virt_accel; cm_virt_accel; pl_virt_accel] + combined_qv;
    
    sv_dot = [cw_vel; cm_vel; pl_vel; true_accels];
end

function [acc, Jac] = virtual_acceleration(tau, pos, mu)
    xhat = [cos(tau); sin(tau); 0];

    ep = -mu * xhat;  mp = (1-mu) * xhat;
    epr = pos - ep;  mpr = pos - mp;

    epn = norm(epr);  mpn = norm(mpr);
    eph = epr / epn;  mph = mpr / mpn;

    earth_acc = -(1 - mu) * eph / epn^2;
    moon_acc = -mu * mph / mpn^2;

    earth_Jac = -(1 - mu) * (eye(3) - 3 * (eph * eph')) / epn^3;
    moon_Jac = -mu * (eye(3) - 3 * (mph * mph')) / mpn^3;

    acc = earth_acc + moon_acc;
    Jac = earth_Jac + moon_Jac;
end