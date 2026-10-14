%% Script setup
clear; clc; close('all');
addpath(genpath('./utilities'))

ode_opts = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);

%% Start defining parameters of the system
mu_e = 398600.4415;
omega = 2 * pi / (22.5 * 60);
L = 300;
m1 = 500;
m2 = 500;

R = 42124;           % Radius of tether C.o.M.'s orbit
V = sqrt(mu_e / R);  % Speed of tether C.o.M.'s orbit

%% Derive the payload and tether's initial condition dimensionally
xs = m2 * L / (m1 + m2);

lhat_0 = [1; 0; 0];

counterweight_relative_position = lhat_0 * (-xs);
catchmech_relative_position = lhat_0 * (L - xs);

omv = [0; 0; omega];

counterweight_delta_state = [counterweight_relative_position; cross(omv, counterweight_relative_position)];
catchmech_delta_state = [catchmech_relative_position; cross(omv, catchmech_relative_position)];

tether_com_IC = [R; 0; 0; 0; V; 0];

cw_IC = tether_com_IC - counterweight_delta_state; 
cm_IC = tether_com_IC - catchmech_delta_state;
pl_IC = cm_IC;

% Reorder IC vector to be [r1; r2; r3; v1; v2; v3]
total_IC = [cw_IC(1:3); cm_IC(1:3); pl_IC(1:3); cw_IC(4:6); cm_IC(4:6); pl_IC(4:6)];
total_IC(7) = total_IC(7) - 4.5 / 1000;

Ldot = 0.001;

total_IC(10:12) = total_IC(10:12) + Ldot / 2 * lhat_0;
total_IC(13:15) = total_IC(13:15) - Ldot / 2 * lhat_0;

%% Nondimensionalize all quantities
L_star = R;
T_star = sqrt(L_star^3 / mu_e);
V_star = L_star / T_star;

mu_t = m2 / (m1 + m2);

total_IC_nd = total_IC;

total_IC_nd(1:9) = total_IC_nd(1:9) / L_star;
total_IC_nd(10:18) = total_IC_nd(10:18) / V_star;

% Augment the initial conditions with an 18 x 18 state transition matrix
sv_0 = [total_IC_nd; reshape(eye(18), [324, 1])];

%% Check Jacobian of the dynamics with finite differencing
ode_func = @(tau, sv) constrained_gravity(tau, sv, mu_t);

ref_sv_dot = ode_func(0, sv_0);
DF_ref = reshape(ref_sv_dot(19:342), [18, 18]);

DF_fd = zeros(18);
fd_stepsize = 1e-6;

for k = 1:1:18
    d_sv_0 = zeros(size(sv_0));
    d_sv_0(k) = fd_stepsize;

    svd_plus = ode_func(0, sv_0 + d_sv_0);
    svd_minus = ode_func(0, sv_0 - d_sv_0);

    DF_fd(:, k) = (svd_plus(1:18) - svd_minus(1:18)) / 2 / fd_stepsize;
end

DF_err = abs(DF_fd - DF_ref);

%% Propagate state and state transition matrix nondimensionally
prop_time = 1000 / T_star;

% Propagate backwards first by five minutes to get the right IC at that time
sol_struct = ode89(ode_func, [0, prop_time], sv_0, ode_opts);

STM_f = reshape(sol_struct.y(19:342, end), [18, 18]);

cw_pos = sol_struct.y(1:3, :); cm_pos = sol_struct.y(4:6, :); pl_pos = sol_struct.y(7:9, :);

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

dxf_dldot = STM_f * [zeros([9, 1]); lhat_0 / 2; -lhat_0 / 2; zeros([3, 1])];

%% Check the Ldot partial derivative

dsv0 = zeros(size(sv_0));

fd_stepsize = 1e-6;

dsv0(10:12) =  fd_stepsize / 2 * lhat_0;
dsv0(13:15) = -fd_stepsize / 2 * lhat_0;

sv0_pert = sv_0 + dsv0;
sol_struct = ode89(ode_func, [0, prop_time], sv0_pert, ode_opts);
final_state_plus = sol_struct.y(1:18, end);

sv0_pert = sv_0 - dsv0;
sol_struct = ode89(ode_func, [0, prop_time], sv0_pert, ode_opts);
final_state_minus = sol_struct.y(1:18, end);

dxf_dldot_fd = (final_state_plus - final_state_minus) / 2 / fd_stepsize;

abs(dxf_dldot - dxf_dldot_fd)


%% Local function definition
function sv_dot = constrained_gravity(~, sv, eta)
    cw_pos = sv(1:3); cm_pos = sv(4:6); pl_pos = sv(7:9); 
    cw_vel = sv(10:12); cm_vel = sv(13:15); pl_vel = sv(16:18);

    stm = reshape(sv(19:342), [18, 18]);

    % Calculate the virtual accelerations of the two point masses
    [cw_virt_accel, cw_A] = virtual_acceleration(cw_pos);
    [cm_virt_accel, cm_A] = virtual_acceleration(cm_pos);
    [pl_virt_accel, pl_A] = virtual_acceleration(pl_pos);

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

function [acc, Jac] = virtual_acceleration(pos)
    rhat = pos / norm(pos);

    acc = -rhat / dot(pos, pos);
    Jac = (3 * (rhat * rhat') - eye(3)) / dot(pos, pos)^(3/2);
end