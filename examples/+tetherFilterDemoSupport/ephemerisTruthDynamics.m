% Copyright (c) 2026 Andrew Binder

function derivative = ephemerisTruthDynamics(tau, state, config, parameters)
    % EPHEMERISTRUTHDYNAMICS Evaluate tethered ephemeris truth dynamics.
    %   derivative = tetherFilterDemoSupport.ephemerisTruthDynamics( ...
    %       tau, state, config, parameters) returns the 18x1 derivative for
    %   [p1; p2; p3; v1; v2; v3] in instantaneous EMB rotating axes.
    %   tau and the state use the fixed nondimensional scales supplied to
    %   emframe.prepareScaledDynamics. Velocities are coordinate derivatives.
    %
    %   config supplies mass_2_fraction, tether_length_rate, and
    %   spring_damper. parameters is the prepared ephemeris dynamics struct.
    %   The existing COA tether model couples bodies 1 and 2; body 3 is free.
    %   This truth-only function does not propagate covariance.
    arguments
        tau (1, 1) double {mustBeReal, mustBeFinite}
        state (18, 1) double {mustBeReal, mustBeFinite}
        config (1, 1) struct
        parameters (1, 1) struct
    end

    accel_fn = @(t, body_state) ephemerisAcceleration(t, body_state, parameters);
    propagate = @coa.orbit.tetheredThreeBodyEkfPropagation;
    derivative = propagate(tau, state, config.mass_2_fraction, ...
                           config.tether_length_rate, accel_fn, struct(), ...
                           config.spring_damper);
end

function acceleration = ephemerisAcceleration(tau, state, parameters)
    derivative = emframe.embRotatingDerivative(tau, state, parameters);
    acceleration = derivative(4:6);
end
