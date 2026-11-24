% Copyright (c) 2026 Andrew Binder

classdef Ephemeris
    % Prepare the filter demo's fixed-scale ephemeris truth model.
    methods (Static)

        function parameters = loadParameters(config, config_path)
            % Read the local DE440 paths and prepare the selected scenario window.
            if ~isfile(config_path)
                error('tetherFilterDemoSupport:Ephemeris:MissingConfig', ...
                      'Create the local ephemeris configuration file: %s.', ...
                      config_path);
            end
            settings = jsondecode(fileread(config_path));
            fields = {'EMFRAME_DE440_HEADER', 'EMFRAME_DE440_ASCII'};
            if ~isstruct(settings) || ~isscalar(settings) || ...
                    ~all(isfield(settings, fields))
                error('tetherFilterDemoSupport:Ephemeris:InvalidConfig', ...
                      'Ephemeris configuration requires DE440 header and ASCII paths.');
            end
            for k = 1:numel(fields)
                value = settings.(fields{k});
                if ~ischar(value) || ~isrow(value) || isempty(strtrim(value))
                    error('tetherFilterDemoSupport:Ephemeris:InvalidConfig', ...
                          'Ephemeris configuration field %s must be a path.', ...
                          fields{k});
                end
            end
            data = emframe.interp.readDe440Ascii(settings.EMFRAME_DE440_HEADER, ...
                                                 settings.EMFRAME_DE440_ASCII);
            parameters = ...
                tetherFilterDemoSupport.Ephemeris.prepareParameters(config, data);
        end

        function parameters = prepareParameters(config, data)
            % Convert DE constants to a force model and prepare scaled dynamics.
            names = {'GMS', 'GM1', 'GM2', 'GMB', 'GM4', 'GM5', 'GM6', ...
                     'GM7', 'GM8', 'GM9', 'AU', 'EMRAT'};
            [present, indices] = ismember(names, data.constantNames);
            if ~all(present)
                error('tetherFilterDemoSupport:Ephemeris:MissingConstant', ...
                      'DE header lacks required constant(s): %s.', ...
                      strjoin(names(~present), ', '));
            end
            values = data.constantValues(indices);
            gm = values(1:10) * values(11)^3 / 86400^2;
            ratio = values(12);
            validateattributes(ratio, {'double'}, ...
                               {'real', 'finite', 'positive', 'scalar'});
            earth_gm = gm(4) * ratio / (1 + ratio);
            moon_gm = gm(4) / (1 + ratio);
            constants = struct('ids', [10, 199, 299, 399, 301, 4:9], ...
                               'gm', [gm(1:3), earth_gm, moon_gm, gm(5:10)]);
            force = emframe.internal.makeForceModel(constants);
            window = [config.ephemeris_start_et, ...
                      config.ephemeris_start_et + config.duration_s];
            scales = struct('length', config.length_scale, 'time', config.time_scale);
            parameters = emframe.prepareScaledDynamics(data, window, force, scales);
        end

    end
end
