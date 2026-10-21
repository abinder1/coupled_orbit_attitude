% Copyright (c) 2026 Andrew Binder

% checkLagrangianDerivation
%   Symbolic check of lagrangian_based.tex. Each check builds a quantity two
%   independent ways and asserts that the difference is exactly zero.
%     Sections 2.1-2.4: constrained equations of motion (exact symbolic
%                       simplification).
%     Section 3:        Jacobian of the dynamics (exact evaluation at a
%                       rational test point, with generic a_i).
%
%   Requires the Symbolic Math Toolbox. Run from any folder.

clear;
clc;

%% Symbols
% p_i: position, v_i: velocity (pdot), pdd_i: acceleration (pddot), all in
% the rotating CR3BP frame. The two masses are weighted (1-eta) and eta.
p1  = sym('p1',  [3 1], 'real');
p2  = sym('p2',  [3 1], 'real');
v1  = sym('v1',  [3 1], 'real');
v2  = sym('v2',  [3 1], 'real');
pdd1 = sym('pdd1', [3 1], 'real');
pdd2 = sym('pdd2', [3 1], 'real');
syms mu eta Ld L0 t lambda real;

w  = [0; 0; 1];                 % CR3BP angular velocity (nondimensional)
Wx = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];  % Wx * x = w x x
rE = [-mu; 0; 0];               % primary (Earth) location
rM = [1 - mu; 0; 0];            % secondary (Moon) location
I3 = eye(3);

%% Section 2.1: Euler-Lagrange equations with a multiplier
% Specific Lagrangian (per unit total mass), as in the "specific energies"
% equation: kinetic terms use the inertial velocity v + w x p.
vin1 = v1 + Wx * p1;
vin2 = v2 + Wx * p2;
T = (1 - eta) / 2 * (vin1.' * vin1) + eta / 2 * (vin2.' * vin2);
V1 = -(1 - eta) * ((1 - mu) / norm(p1 - rE) + mu / norm(p1 - rM));
V2 = -eta * ((1 - mu) / norm(p2 - rE) + mu / norm(p2 - rM));
V = V1 + V2;
Lag = T - V;

% Constraint g and its gradients (row vectors, as written in the .tex).
d = p2 - p1;
g = d.' * d - (Ld * t + L0)^2;
dgdp1 = jacobian(g, p1);
dgdp2 = jacobian(g, p2);

% Check: the gradients are -2 d' and +2 d' (stated after the d^2g/dt^2 steps).
assertZero(dgdp1 + 2 * d.', 'dg/dp1 = -2 d''');
assertZero(dgdp2 - 2 * d.', 'dg/dp2 = +2 d''');

% d/dt (dL/dpdot) by the chain rule: only p and pdot carry time dependence.
dLdv1 = jacobian(Lag, v1).';
dLdv2 = jacobian(Lag, v2).';
ddt = @(f, p, v, pdd) jacobian(f, p) * v + jacobian(f, v) * pdd;
momRate1 = ddt(dLdv1, p1, v1, pdd1);
momRate2 = ddt(dLdv2, p2, v2, pdd2);

dLdp1 = jacobian(Lag, p1).';
dLdp2 = jacobian(Lag, p2).';

%% Section 2.2: momentum rate and unconstrained accelerations
% Check: d/dt dL/dpdot_i = mbar_i (pddot_i + w x pdot_i), with mbar_1 = 1-eta
% and mbar_2 = eta.
assertZero(momRate1 - (1 - eta) * (pdd1 + Wx * v1), 'momentum rate, mass 1');
assertZero(momRate2 - eta * (pdd2 + Wx * v2), 'momentum rate, mass 2');

% Check: the EL equations rearrange to
%   mbar_i pddot_i - lambda dg/dp_i' = dL/dp_i - mbar_i (w x pdot_i).
% Both sides are built independently: the left from the EL residual
% (momRate - dL/dp - lambda dg/dp'), the right from the rearranged form.
el1 = momRate1 - dLdp1 - lambda * dgdp1.';
el2 = momRate2 - dLdp2 - lambda * dgdp2.';
rearr1 = (1 - eta) * pdd1 - lambda * dgdp1.' - (dLdp1 - (1 - eta) * (Wx * v1));
rearr2 = eta * pdd2 - lambda * dgdp2.' - (dLdp2 - eta * (Wx * v2));
assertZero(el1 - rearr1, 'rearranged EL, mass 1');
assertZero(el2 - rearr2, 'rearranged EL, mass 2');

% Unconstrained accelerations (lambda = 0): mbar_i a_i = dL/dp_i - mbar_i w x pdot_i.
a1 = (dLdp1 - (1 - eta) * (Wx * v1)) / (1 - eta);
a2 = (dLdp2 - eta * (Wx * v2)) / eta;

% Check: a_i equals the standard CR3BP acceleration of a single particle.
cr3bpAcc = @(p, v) -(1 - mu) * (p - rE) / norm(p - rE)^3 - ...
                   mu * (p - rM) / norm(p - rM)^3 - Wx * (Wx * p) - 2 * Wx * v;
assertZero(a1 - cr3bpAcc(p1, v1), 'a1 = CR3BP acceleration');
assertZero(a2 - cr3bpAcc(p2, v2), 'a2 = CR3BP acceleration');

%% Section 2.3: second derivative of the constraint and the linear system
% Check d^2 g/dt^2 using genuine time-dependent positions, so that diff()
% differentiates through p(t) independently of the formula in the .tex.
p1t = [str2sym('p1x(t)'); str2sym('p1y(t)'); str2sym('p1z(t)')];
p2t = [str2sym('p2x(t)'); str2sym('p2y(t)'); str2sym('p2z(t)')];
dt_ = p2t - p1t;
gt = dt_.' * dt_ - (Ld * t + L0)^2;
ddot = diff(dt_, t);
ddd_acc = diff(p2t, t, 2) - diff(p1t, t, 2);
claimed = 2 * (ddot.' * ddot) + 2 * dt_.' * ddd_acc - 2 * Ld^2;
assertZero(diff(gt, t, 2) - claimed, 'd^2g/dt^2 formula');

% Linear system A * [pdd1; pdd2; lambda] = rhs (the 7x7 matrix equation).
% Generic symbols are used for d, ddot (= v2 - v1), and a_i so the result
% does not depend on the specific forms of the accelerations.
ds   = sym('d',  [3 1], 'real');
dds  = sym('dd', [3 1], 'real');
as1  = sym('a1_', [3 1], 'real');
as2  = sym('a2_', [3 1], 'real');
A = [(1 - eta) * I3, zeros(3),  2 * ds
     zeros(3),       eta * I3, -2 * ds
     -2 * ds.',      2 * ds.',  0];
rhs = [(1 - eta) * as1; eta * as2; 2 * Ld^2 - 2 * (dds.' * dds)];

% Check: rows 1-6 are the rearranged EL equations (with dg/dp = -/+ 2 d')
% and row 7 is the rearranged acceleration-level constraint.
x = [pdd1; pdd2; lambda];
rows16 = A(1:6, :) * x - rhs(1:6);
expect16 = [(1 - eta) * (pdd1 - as1) + lambda * 2 * ds
            eta * (pdd2 - as2) - lambda * 2 * ds];
assertZero(rows16 - expect16, 'rows 1-6 of the linear system');
row7 = A(7, :) * x - rhs(7);
expect7 = 2 * ds.' * (pdd2 - pdd1) - (2 * Ld^2 - 2 * (dds.' * dds));
assertZero(row7 - expect7, 'row 7 of the linear system');

%% Section 2.4: closed-form constrained accelerations
% Solve the 7x7 system symbolically and compare against the closed form
% with D = d d'/(d'd) and q = (ddot'ddot - Ld^2)/(d'd) d.
sol = simplify(A \ rhs);
D = (ds * ds.') / (ds.' * ds);
q = ((dds.' * dds - Ld^2) / (ds.' * ds)) * ds;
closed = [(I3 - eta * D) * as1 + eta * D * as2 + eta * q
          (1 - eta) * D * as1 + (I3 - (1 - eta) * D) * as2 - (1 - eta) * q];
assertZero(sol(1:6) - closed, 'closed-form constrained accelerations');

% Extra check: the constrained accelerations satisfy the constraint, i.e.
% d'(pdd2 - pdd1) = Ld^2 - ddot'ddot.
assertZero(ds.' * (closed(4:6) - closed(1:3)) - (Ld^2 - dds.' * dds), ...
           'closed form satisfies the constraint');

disp('All Section 2.1-2.4 checks passed.');

%% Section 3: Jacobian of the constrained dynamics
% Ground truth: build the constrained accelerations explicitly from
% Section 3.1 and let jacobian() differentiate them. Each formula in the
% .tex is then compared against this (or against a direct derivative of a
% sub-expression). Generic nonlinear accelerations a1(p1,v1) and a2(p2,v2)
% are used so the identities are shown to hold for arbitrary a_i, not just
% the CR3BP form. a_i depends only on its own mass's state, as in the .tex.
ag1 = [p1(1) * v1(2)^2 + sin(p1(2)); ...
       p1(3) * p1(1) + v1(1)^3; ...
       exp(p1(1)) * v1(3) + p1(2)^2 * v1(1)];
ag2 = [p2(2) * v2(1) * v2(3) + p2(3)^2; ...
       cos(p2(1)) + v2(2)^2 * p2(3); ...
       p2(1) * p2(2) * v2(2) + v2(3)^3];

dv  = p2 - p1;                          % d
ddv = v2 - v1;                          % ddot
dd_ = dv.' * dv;                        % d . d
Dm  = (dv * dv.') / dd_;                % projector D
kap = ddv.' * ddv - Ld^2;               % kappa
qv  = kap / dd_ * dv;                   % q

% Section 3.1: the expanded form equals the block form of Section 2.4, and
% D(a2 - a1) = [d'(a2 - a1) / d'd] d.
pdd1 = ag1 + eta * (Dm * (ag2 - ag1) + qv);
pdd2 = ag2 - (1 - eta) * (Dm * (ag2 - ag1) + qv);
blk1 = (I3 - eta * Dm) * ag1 + eta * Dm * ag2 + eta * qv;
blk2 = (1 - eta) * Dm * ag1 + (I3 - (1 - eta) * Dm) * ag2 - (1 - eta) * qv;
X = [p1; p2; v1; v2; eta; Ld];          % every free symbol in Section 3
% Exact rational test point
Xv = [sym([3 -2 5 7 1 4 -3 2 6 -1 8 3].') / 7; 2 / sym(7); 3 / sym(5)];
check = @(expr, name) assertZeroAt(expr, name, X, Xv);
check([pdd1; pdd2] - [blk1; blk2], '3.1 expanded form = block form');
check(Dm * (ag2 - ag1) - (dv.' * (ag2 - ag1)) / dd_ * dv, '3.1 D(a2-a1) scalar form');

% Section 3.2: derivatives of D*a_i. Partial derivatives of d are -I (p1)
% and +I (p2); the general formula is checked for every (i, j) pair, then
% the four specialized position blocks and the velocity blocks.
P = {p1, p2};
V = {v1, v2};
A = {ag1, ag2};
dDp = {-I3, I3};
for i = 1:2
    for j = 1:2
        % d(d.a_i)/dp_j and d((d.d)^-1)/dp_j (stated before the general result)
        d_da_dot = A{i}.' * dDp{j} + dv.' * jacobian(A{i}, P{j});
        check(jacobian(dv.' * A{i}, P{j}) - d_da_dot, ...
              sprintf('3.2 d(d.a%d)/dp%d', i, j));
        check(jacobian(1 / dd_, P{j}) + 2 / dd_^2 * dv.' * dDp{j}, ...
              sprintf('3.2 d((d.d)^-1)/dp%d', j));
        % General formula for d(D a_i)/dp_j
        gen = Dm * jacobian(A{i}, P{j}) + ...
              (dv * A{i}.' + (dv.' * A{i}) * (I3 - 2 * Dm)) / dd_ * dDp{j};
        check(jacobian(Dm * A{i}, P{j}) - gen, ...
              sprintf('3.2 d(D a%d)/dp%d general', i, j));
        % Velocity derivatives: D depends only on positions, and a_i only
        % on its own velocity.
        check(jacobian(Dm * A{i}, V{j}) - Dm * jacobian(A{i}, V{j}), ...
              sprintf('3.2 d(D a%d)/dv%d', i, j));
    end
end
% Specialized position blocks written with T1..T4 (Section 3.4 definitions).
T1 = (dv * ag1.') / dd_;
T2 = (dv.' * ag1) / dd_ * (I3 - 2 * Dm);
T3 = (dv * ag2.') / dd_;
T4 = (dv.' * ag2) / dd_ * (I3 - 2 * Dm);
check(jacobian(Dm * ag1, p1) - (Dm * jacobian(ag1, p1) - T1 - T2), '3.2 d(D a1)/dp1');
check(jacobian(Dm * ag1, p2) - (T1 + T2),                           '3.2 d(D a1)/dp2');
check(jacobian(Dm * ag2, p1) + (T3 + T4),                           '3.2 d(D a2)/dp1');
check(jacobian(Dm * ag2, p2) - (Dm * jacobian(ag2, p2) + T3 + T4),  '3.2 d(D a2)/dp2');

% Section 3.3: derivatives of q with respect to position, velocity, and Ld.
K = kap / dd_ * (I3 - 2 * Dm);          % K = kappa/(d.d) (I - 2D)
W = 2 / dd_ * (dv * ddv.');             % W = 2/(d.d) d ddot'
check(jacobian(qv, p1) + K, '3.3 dq/dp1 = -K');
check(jacobian(qv, p2) - K, '3.3 dq/dp2 = +K');
check(jacobian(qv, v1) + W, '3.3 dq/dv1 = -W');
check(jacobian(qv, v2) - W, '3.3 dq/dv2 = +W');
check(diff(qv, Ld) + 2 * Ld / dd_ * dv, '3.3 dq/dLd');

% Section 3.4: the full 6x12 Jacobian of [pdd1; pdd2] with respect to
% [p1; p2; v1; v2], against both tables in the .tex.
Jtrue = jacobian([pdd1; pdd2], [p1; p2; v1; v2]);

% Table 1: assembled from d(D a_i) and dq blocks.
Da1p1 = Dm * jacobian(ag1, p1) - T1 - T2;
Da1p2 = T1 + T2;
Da2p1 = -(T3 + T4);
Da2p2 = Dm * jacobian(ag2, p2) + T3 + T4;
q_p1 = -K;
q_p2 = K;
q_v1 = -W;
q_v2 = W;
Z3 = zeros(3);
J1 = [jacobian(ag1, p1) + eta * (Da2p1 - Da1p1 + q_p1), ...
      eta * (Da2p2 - Da1p2 + q_p2), ...
      jacobian(ag1, v1) + eta * (q_v1 - Dm * jacobian(ag1, v1)), ...
      eta * (q_v2 + Dm * jacobian(ag2, v2))];
J2 = [-(1 - eta) * (Da2p1 - Da1p1 + q_p1), ...
      jacobian(ag2, p2) - (1 - eta) * (Da2p2 - Da1p2 + q_p2), ...
      -(1 - eta) * (q_v1 - Dm * jacobian(ag1, v1)), ...
      jacobian(ag2, v2) - (1 - eta) * (q_v2 + Dm * jacobian(ag2, v2))];
check(Jtrue - [J1; J2], '3.4 Table 1 (all eight blocks)');

% Table 2: grouped by external Jacobian, using S = T3 + T4 - T1 - T2.
S = (dv * (ag2 - ag1).' + (dv.' * (ag2 - ag1)) * (I3 - 2 * Dm)) / dd_;
check(S - (T3 + T4 - T1 - T2), '3.4 S definition');
Ap1 = jacobian(ag1, p1);
Ap2 = jacobian(ag2, p2);
Av1 = jacobian(ag1, v1);
Av2 = jacobian(ag2, v2);
G1 = [(I3 - eta * Dm) * Ap1 - eta * (S + K), ...
      eta * Dm * Ap2 + eta * (S + K), ...
      (I3 - eta * Dm) * Av1 - eta * W, ...
      eta * Dm * Av2 + eta * W];
G2 = [(1 - eta) * Dm * Ap1 + (1 - eta) * (S + K), ...
      (I3 - (1 - eta) * Dm) * Ap2 - (1 - eta) * (S + K), ...
      (1 - eta) * Dm * Av1 + (1 - eta) * W, ...
      (I3 - (1 - eta) * Dm) * Av2 - (1 - eta) * W];
check(Jtrue - [G1; G2], '3.4 Table 2 (grouped form, all eight blocks)');

disp('All Section 3 checks passed.');

%% Local helpers
function assertZero(expr, name)
    % Simplify and require an identically zero symbolic result.
    r = simplify(expr);
    assert(all(r(:) == 0), 'Check failed: %s (residual: %s)', name, char(r(1)));
    fprintf('  ok: %s\n', name);
end

function assertZeroAt(expr, name, vars, vals)
    % Section 3 expressions are large rational functions, so instead of
    % simplifying symbolically we substitute an exact rational test point
    % and require an exactly zero result. A nonzero rational function
    % vanishes at a generic point with probability zero, so a pass here is
    % strong evidence the identity holds (and a fail is a definite error).
    r = subs(expr, vars, vals);
    assert(all(r(:) == 0), 'Check failed: %s (residual: %s)', name, char(r(1)));
    fprintf('  ok: %s\n', name);
end
