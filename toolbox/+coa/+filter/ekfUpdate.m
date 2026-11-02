% Copyright (c) 2026 Andrew Binder

function [state, covariance, info] = ekfUpdate(prior_state, prior_covariance, ...
                                               observation, measurement_model, ...
                                               angle_rows)
    % EKFUPDATE Joseph-form extended Kalman filter measurement update.
    %   [state, covariance] = coa.filter.ekfUpdate(prior_state, ...
    %       prior_covariance, observation, measurement_model)
    %   corrects a prior estimate with one observation. The function is
    %   deterministic and uses neither truth nor random numbers, so it applies
    %   equally to simulated and real measurements.
    %
    %   [___, info] = coa.filter.ekfUpdate(___) also returns the quantities
    %   used for filter-health diagnostics (see Outputs).
    %
    %   [___] = coa.filter.ekfUpdate(prior_state, prior_covariance, ...
    %       observation, measurement_model, angle_rows) wraps the listed
    %   innovation and residual rows to [-pi, pi). Without this, an angle
    %   near the atan2 branch cut produces a spurious 2*pi innovation.
    %
    %   With H the measurement Jacobian at the prior and R the measurement
    %   noise covariance, the update is
    %       nu  = z - h(x^-)                           (innovation)
    %       S   = H*P^-*H' + R                         (innovation covariance)
    %       K   = P^-*H'*S^-1                          (gain)
    %       x^+ = x^- + K*nu
    %       P^+ = (I - K*H)*P^-*(I - K*H)' + K*R*K'    (Joseph form)
    %   The Joseph form stays symmetric and positive semidefinite under
    %   roundoff and for suboptimal gains. P^+ is symmetrized before return.
    %   See derivations/ekf_health_diagnostics.tex for the diagnostics built
    %   on info.
    %
    %   Inputs:
    %     prior_state         n-by-1 prior estimate x^-.
    %     prior_covariance    n-by-n prior covariance P^-, symmetric.
    %     observation         m-by-1 measurement z, in the same row order,
    %                         frame, and units as the model output.
    %     measurement_model   Function handle,
    %                             [h, H, R] = measurement_model(state)
    %                         returning the m-by-1 predicted measurement, the
    %                         m-by-n Jacobian dh/dstate, and the m-by-m noise
    %                         covariance. It is called at the prior state
    %                         (all three outputs) and at the posterior state
    %                         (predicted measurement only, for the post-fit
    %                         residual), so it must be free of side effects.
    %                         For coa.meas.measurementModels, bind the mask
    %                         and covariances:
    %                             @(x) coa.meas.measurementModels(x, mask, covs)
    %     angle_rows          (optional) Row indices of angular measurements
    %                         to wrap. Default: none.
    %
    %   Outputs:
    %     state       n-by-1 corrected (posterior) estimate x^+.
    %     covariance  n-by-n corrected (posterior) covariance P^+.
    %     info        Scalar struct:
    %       innovation              m-by-1 pre-fit innovation nu (wrapped).
    %       innovation_covariance   m-by-m S. The normalized innovation
    %                               squared nu'*(S\nu) is chi-square with m
    %                               degrees of freedom for a consistent filter.
    %       gain                    n-by-m Kalman gain K.
    %       residual                m-by-1 post-fit residual z - h(x^+)
    %                               (wrapped).
    %       residual_covariance     m-by-m R - H*P^+*H', the post-fit residual
    %                               covariance for the optimal gain, with H
    %                               the prior-point Jacobian.
    %
    %   Notes:
    %     - S must be positive definite; otherwise the function raises
    %       coa:filter:ekfUpdate:SingularInnovationCovariance.
    %     - H is evaluated once, at the prior state; the update is not
    %       iterated.
    %     - A zero-row measurement returns the prior unchanged.
    %
    %   Example:
    %       model = @(x) coa.meas.measurementModels(x, mask, covs);
    %       [x_post, P_post, info] = coa.filter.ekfUpdate(x_prior, P_prior, ...
    %                                                     z, model, 1:6);
    %       nis = info.innovation' * (info.innovation_covariance \ info.innovation);
    arguments
        prior_state (:, 1) double {mustBeFinite}
        prior_covariance (:, :) double {mustBeFinite}
        observation (:, 1) double {mustBeFinite}
        measurement_model (1, 1) function_handle
        angle_rows (:, 1) double {mustBeInteger, mustBePositive} = zeros(0, 1)
    end

    state_count = numel(prior_state);
    if ~isequal(size(prior_covariance), [state_count, state_count])
        error('coa:filter:ekfUpdate:SizeMismatch', ...
              'prior_covariance must be %d-by-%d to match prior_state.', ...
              state_count, state_count);
    end

    [predicted_measurement, jacobian, noise_covariance] = ...
        measurement_model(prior_state);

    measurement_count = numel(predicted_measurement);
    if numel(observation) ~= measurement_count || ...
            ~isequal(size(jacobian), [measurement_count, state_count]) || ...
            ~isequal(size(noise_covariance), ...
                     [measurement_count, measurement_count])
        error('coa:filter:ekfUpdate:SizeMismatch', ...
              ['observation (%d rows), Jacobian, and noise covariance must ' ...
               'match the %d-row measurement model output and the %d-state ' ...
               'estimate.'], numel(observation), measurement_count, state_count);
    end
    if any(angle_rows > measurement_count)
        error('coa:filter:ekfUpdate:SizeMismatch', ...
              'angle_rows must not exceed the %d measurement rows.', ...
              measurement_count);
    end

    innovation = wrapRows(observation - predicted_measurement, angle_rows);
    cross_covariance = jacobian * prior_covariance;
    innovation_covariance = cross_covariance * jacobian' + noise_covariance;

    [~, not_positive_definite] = chol(innovation_covariance);
    if not_positive_definite
        error('coa:filter:ekfUpdate:SingularInnovationCovariance', ...
              'The innovation covariance is not positive definite.');
    end

    gain = (innovation_covariance \ cross_covariance)';
    state = prior_state + gain * innovation;

    correction = eye(state_count) - gain * jacobian;
    covariance = correction * prior_covariance * correction' + ...
                 gain * noise_covariance * gain';
    covariance = (covariance + covariance') / 2;

    postfit = observation - measurement_model(state);

    info.innovation = innovation;
    info.innovation_covariance = innovation_covariance;
    info.gain = gain;
    info.residual = wrapRows(postfit, angle_rows);
    info.residual_covariance = noise_covariance - ...
                               jacobian * covariance * jacobian';
end

function values = wrapRows(values, rows)
    % Wrap the selected rows to [-pi, pi).
    values(rows) = mod(values(rows) + pi, 2 * pi) - pi;
end
