% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function set_h5_tseries_attributes(input_file, varargin)
%SET_H5_TSERIES_ATTRIBUTES function that sets the attributes in a begonia
%tseries h5 file
%data is initialized with an optional mdf file or default values and
%completed by the supplied parameters
%parameter names are identical with the Begonia parameter names

    p = inputParser;
    addParameter(p, 'mdf_reference', [], @(x) isa(x, 'MDF_file_reader') || ...
        (@(x) isstring(x) && ...
        @(x) isfile(x)));
    addParameter(p, 'name', []);
    addParameter(p, 'frame_count', []);
    addParameter(p, 'slices', []);
    addParameter(p, 'channel_names', []);
    addParameter(p, 'channels', []);
    addParameter(p, 'img_dim', []);
    addParameter(p, 'dt', []);
    addParameter(p, 'dx', []);
    addParameter(p, 'dy', []);
    addParameter(p, 'zoom', []);
    addParameter(p, 'start_time', []);
    
    parse(p, varargin{:})
    
    [~, ~, ext] = fileparts(input_file);
    if ~strcmp(ext, '.h5')
        warning('Begonia only accepts files with the .h5 extension. To open the file in Begonia manually change the extension.');
    end
    
    info = h5info(input_file);
    idx = strcmp({info.Datasets.Name}, 'recording');
    assert(~isempty(idx), ...
        'Please pass a h5 file that contains already a dataset recording, e.g. with convert_to_tseriesh5 or TSERIESH5_file_writer!');

    tmp = num2cell(info.Datasets(idx).Dataspace.Size);
    [m, n, n_channels, n_frames] = tmp{:};
    channel_names = cell(1, n_channels);
    for i = 1:n_channels
        channel_names{i} = ['ch' int2str(i)];
    end
    
    if ~isempty(p.Results.mdf_reference)
        % MDF mode:
        if isa(p.Results.mdf_reference, 'MDF_file_reader')
            mdf_file = p.Results.mdf_reference;
        else
            mdf_file = MDF_file_reader(p.Results.mdf_reference);
        end
        metadata = mdf_file.get_tseries_metadata();
    else
        % Default values:
        metadata.name = 'recording';
        metadata.source = 'FlowReg';
        metadata.frame_count = n_frames;
        metadata.slices = n_frames ;
        metadata.channel_names = channel_names;
        metadata.channels = n_channels;
        metadata.img_dim = [m, n];
        % choose default fps of 30.9:
        metadata.dt = 1 / 30.9;
        metadata.dx = 1;
        metadata.dy = 1;
        metadata.zoom = 1;
        tmp = datetime();
        tmp.Format = 'uuuu/MM/dd HH:mm:ss';
        metadata.start_time = tmp;
    end
    
    for i = 1:length(p.Parameters)
        pm_id =  p.Parameters{i};
        param = p.Results.(pm_id);
        if (strcmp({'mdf_reference'}, pm_id) || isempty(param))
            continue;
        end
        metadata.(pm_id) = param;
    end
    metadata.duration = metadata.dt * n_frames;
    if isdatetime(metadata.start_time)
        metadata.start_time.Format = 'uuuu/MM/dd HH:mm:ss';
    end
    
    h5writeatt(input_file, '/recording', ...
        'json_metadata', jsonencode(metadata));
end

