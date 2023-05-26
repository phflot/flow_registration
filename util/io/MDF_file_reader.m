% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

% class for reading blocks from an mdf file with binning
classdef MDF_file_reader < Video_file_reader
    
    properties (Access = private)
        mfile;
        channel_idx = [1, 2];
        channels;
    end
    
    properties(GetAccess = public, SetAccess = private)
        file_name;
    end
    
    properties(GetAccess = public, Constant)
        datatype = 'MDF';
    end
    
    properties
        current_frame = 0;
    end
    
    methods
        function obj = MDF_file_reader(input_file, buffer_size, ...
                bin_size, varargin)
            
            if nargin > 1 && ~isempty(buffer_size)
                obj.buffer_size = buffer_size;
            end
            if nargin > 2 && ~isempty(bin_size)
                obj.bin_size = bin_size;
            end
            
            obj.file_name = input_file;

            [~, f_name, ~] = fileparts(input_file);
            obj.input_file_name = f_name;
                        
            obj.buffer_size = buffer_size;
            obj.bin_size = bin_size;
            
            % previously inside of the open function:
            obj.mfile = actxserver('MCSX.Data',[0 0 0 0]);
            if obj.mfile.invoke('OpenMCSFile', input_file)
                error('Only one MDF file instance can be opened at once! E.g. close the MDF Viewer and clear the Matlab workspace.');
            end

            obj.frame_count = str2num(invoke(obj.mfile, 'ReadParameter', 'Frame Count'));
            obj.m = str2num(invoke(obj.mfile, 'ReadParameter', 'Frame Height'));
            obj.n = str2num(invoke(obj.mfile, 'ReadParameter', 'Frame Width'));
            tmpbit = strsplit(invoke(obj.mfile, 'ReadParameter', 'Frame Bit Depth'), '-');
            obj.bitdepth = str2num(tmpbit{1});
            
            obj.channels = [];
            for i = 0:2
                if ~isempty(invoke(obj.mfile, 'ReadParameter', ['Scanning Ch ' num2str(i) ' Name']))
                    obj.channels(end + 1) = i + 1;
                end
            end
            obj.channel_idx = obj.channels;
            
            p = inputParser;
            addParameter(p, 'channel_idx', obj.channel_idx, ...
                @(x) isvector(x) && isnumeric(x));
            parse(p, varargin{:});
            obj.channel_idx = p.Results.channel_idx;
            
            obj.n_channels = min(length(obj.channels), ...
                length(obj.channel_idx));
            obj.mat_data_type = obj.get_mat_data_type(obj.bitdepth);
        end
        
        function reset(obj)
            obj.current_frame = 0;
        end
        
        function buffer = read_batch(obj)
            if obj.current_frame > obj.frame_count
                buffer = [];
                return;
            end
            
            n_elem_left = min(obj.buffer_size * obj.bin_size, ...
                obj.frame_count - obj.current_frame);
            buffer = zeros(obj.m, obj.n, obj.n_channels, n_elem_left, ...
                obj.mat_data_type);
            
            for i = 1:n_elem_left
                obj.current_frame = obj.current_frame + 1;
                
                for j = 1:obj.n_channels
                    buffer(:, :, j, i) = cast(...
                        obj.mfile.ReadFrame(obj.channel_idx(j), ...
                        obj.current_frame), ...
                        obj.mat_data_type)';
                end
            end

            buffer = obj.bin_buffer(buffer);
        end
        
        function buffer = read_frames(obj, idx)
            
            assert(sum(idx > obj.frame_count) == 0);
            
            n_elements = length(idx);
            
            buffer = zeros(obj.m, obj.n, obj.n_channels, n_elements, obj.mat_data_type);
            
            for i = 1:n_elements             
                for j = 1:obj.n_channels
                    buffer(:, :, j, i) = cast(...
                        obj.mfile.ReadFrame(obj.channel_idx(j), idx(i)), ...
                        obj.mat_data_type)';
                end
            end
            
            buffer = obj.bin_buffer(buffer);
        end
        
        function result = has_batch(obj)
            result = obj.current_frame < obj.frame_count;
        end
        
        function left = batches_left(obj)
            left = ceil((obj.frame_count - obj.current_frame) / ...
                (obj.buffer_size * obj.bin_size));
        end
        
        function delete(obj)
            delete(obj.mfile);
        end
        
        function metadata = get_tseries_metadata(obj)
        
            microns_per_pixel = invoke(...
                obj.mfile, 'ReadParameter', 'Microns per Pixel');
            frame_duration = invoke(...
                obj.mfile, 'ReadParameter', 'Frame Duration (s)');
            frame_interval = invoke(...
                obj.mfile, 'ReadParameter', 'Frame Interval (ms)');
            magnification = invoke(...
                obj.mfile, 'ReadParameter', 'Magnification');
            
            magnification = str2double(...
                replace(replace(magnification, ',', '.'), 'x', ''));            
            microns_per_pixel = str2double(...
                replace(replace(microns_per_pixel, ',', '.'), ' µm', ''));
            frame_duration = str2double(...
                replace(replace(frame_duration, ',', '.'), ' s', ''));
            frame_interval = str2double(...
                replace(replace(frame_interval, ',', '.'), ' ms', ''));
            
            if isempty(microns_per_pixel)
                microns_per_pixel = 1;
                warning('microns per pixel could not be read from mdf');
            end
            if isempty(magnification)
                magnification = 1;
                warning('magnification could not be read from mdf');
            end
            if isempty(frame_duration) || isempty(frame_interval)
                dt = 1 / 30.91;
                warning('frame duration and frame interval could not be read from mdf');
            else                
                dt = frame_duration + frame_interval;
            end
            
            channel_names = cell(1, obj.n_channels);
            for i = 1:obj.n_channels
                channel_names{i} = invoke(...
                obj.mfile, 'ReadParameter', ...
                ['Scanning Ch ' int2str(i - 1) ' Name']);
            end
            
            [~, recording_name, ~] = fileparts(obj.file_name);
            
            metadata.name = recording_name;
            metadata.frame_count = obj.frame_count;
            metadata.slices = obj.frame_count;
            metadata.channel_names = channel_names;
            metadata.channels = obj.n_channels;
            metadata.img_dim = [obj.m, obj.n];
            metadata.dt = dt * obj.bin_size;
            metadata.dx = microns_per_pixel;
            metadata.dy = microns_per_pixel;
            metadata.zoom = magnification;
            metadata.source = 'FlowReg';
            tmp = datetime(...
                invoke(obj.mfile, 'ReadParameter', 'Created On'), ...
                'InputFormat', 'eeee, MMMM d, yyyy h:mm:ss a');
            tmp.Format = 'uuuu/MM/dd HH:mm:ss';
            metadata.start_time = tmp;
        end
    end
end

