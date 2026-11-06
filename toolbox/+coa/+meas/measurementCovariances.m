% Copyright (c) 2026 Andrew Binder

function meas_covs = measurementCovariances(length_scale, speed_scale)
    % MEASUREMENTCOVARIANCES Nominal sensor noise covariances for measurementModels.
    %   meas_covs = coa.meas.measurementCovariances(length_scale, speed_scale)
    %   returns the struct of noise covariances consumed by
    %   coa.meas.measurementModels, built from nominal 1-sigma sensor accuracies.
    %   length_scale (km) and speed_scale (km/s) are the characteristic scales
    %   that nondimensionalize the EKF state. The returned fields are variances
    %   or covariance matrices: angular fields are in rad^2, range fields in
    %   squared nondimensional length units, and range-rate fields in squared
    %   nondimensional speed units. These are measurement-noise covariances;
    %   they do not model sensor bias or bias drift. Only selected groups'
    %   fields need to be passed to coa.meas.measurementModels.
    %
    %   Fields and nominal 1-sigma values:
    %       ground_angles         2x2  10 arcsec, ground optical telescope
    %       camera_angles         2x2  30 arcsec, onboard camera
    %       star_tracker_angles   2x2  30 arcsec, star tracker cross-boresight
    %       camera_angular_size  1x1  sqrt(2) * 30 arcsec, two independent
    %                                      camera-angle errors
    %       dsn_range             1x1  10 m, ground DSN-class ranging
    %       dsn_range_rate        1x1  10 mm/s, ground DSN-class Doppler
    %       crosslink_range       1x1  5 m, RF crosslink code ranging
    %       crosslink_range_rate  1x1  100 mm/s, RF crosslink Doppler
    %       gnss_position         3x3  10 m per axis, lunar PN&T/GNSS fix
    arguments
        length_scale (1, 1) double {mustBeFinite, mustBePositive}
        speed_scale (1, 1) double {mustBeFinite, mustBePositive}
    end

    % Unit conversions into radians and nondimensional state units
    arcsec = pi / 180 / 3600;
    meter = 1e-3 / length_scale;
    meter_per_sec = 1e-3 / speed_scale;

    sigma_ground_angles = 10 * arcsec;
    sigma_camera_angles = 30 * arcsec;
    sigma_star_tracker_angles = 30 * arcsec;
    sigma_dsn_range = 10 * meter;
    sigma_dsn_range_rate = 1e-2 * meter_per_sec;
    sigma_crosslink_range = 5 * meter;
    sigma_crosslink_range_rate = 0.1 * meter_per_sec;
    sigma_gnss_position = 10 * meter;

    meas_covs.ground_angles = sigma_ground_angles^2 * eye(2);
    meas_covs.camera_angles = sigma_camera_angles^2 * eye(2);
    meas_covs.star_tracker_angles = sigma_star_tracker_angles^2 * eye(2);
    meas_covs.camera_angular_size = 2 * sigma_camera_angles^2;
    meas_covs.dsn_range = sigma_dsn_range^2;
    meas_covs.dsn_range_rate = sigma_dsn_range_rate^2;
    meas_covs.crosslink_range = sigma_crosslink_range^2;
    meas_covs.crosslink_range_rate = sigma_crosslink_range_rate^2;
    meas_covs.gnss_position = sigma_gnss_position^2 * eye(3);
end
