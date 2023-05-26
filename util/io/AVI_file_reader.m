% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

% class for reading blocks from an mdf file with binning
classdef AVI_file_reader < Video_file_reader
    
    properties (Access = private)
        video_reader;
    end
    
    properties (GetAccess = public, Constant)
        datatype = 'VIDEO';
    end
    
    properties
        current_frame = 0;
    end
    
    methods
        function obj = AVI_file_reader(input_file, buffer_size, ...
                bin_size, varargin)
                                   
            if nargin > 1 && ~isempty(buffer_size)
                obj.buffer_size = buffer_size;
            end
            if nargin > 2 && ~isempty(bin_size)
                obj.bin_size = bin_size;
            end
            
            [~, f_name, ~] = fileparts(input_file);
            obj.input_file_name = f_name;

            obj.video_reader = VideoReader(input_file);
            obj.frame_count = obj.video_reader.NumberOfFrames;
            
            obj.m = obj.video_reader.Height;
            obj.n = obj.video_reader.Width;
            
            tmp_format = split(obj.video_reader.VideoFormat, ' ');
            switch tmp_format{1}
                case 'RGB24'
                    obj.n_channels = 3;
                    obj.bitdepth = 8;
                case 'RGB48'
                    obj.n_channels = 3;
                    obj.bitdepth = 16;
                case 'Indexed'
                    error('Indexed video files are currently not supported, export data as hdf5 of mat files first instead!');
                case 'Grayscale'
                    obj.n_channels = 1;
                    obj.bitdepth = obj.video_reader.BitsPerPixel;
                otherwise
                    obj.n_channels = 1;
                    obj.bitdepth = obj.video_reader.BitsPerPixel;
            end
            
            obj.mat_data_type = obj.get_mat_data_type(obj.bitdepth, ...
                length(tmp_format) > 1 && strcmp(tmp_format{2}, 'Signed'));
        end
        
        function reset(obj)
            obj.video_reader.CurrentTime = 0;
            obj.current_frame = 0;
        end
        
        function buffer = read_batch(obj)
            if obj.current_frame > obj.frame_count
                buffer = [];
                return;
            end
            
            n_elem_left = min(obj.buffer_size * obj.bin_size, ...
                obj.frame_count - obj.current_frame);
            buffer = zeros(obj.m, obj.n, obj.n_channels, n_elem_left, obj.mat_data_type);
            
            for i = 1:n_elem_left
                obj.current_frame = obj.current_frame + 1;
                
                buffer(:, :, :, i) =...
                    cast(obj.video_reader.readFrame, obj.mat_data_type);
            end
            
            if obj.bin_size > 1
                buffer = cast(convn(buffer, obj.downsampling_kernel, 'same'), obj.mat_data_type);
                buffer = buffer(:, :, :, ceil(obj.bin_size / 2):obj.bin_size:end);
            end
        end
        
        function buffer = read_frames(obj, idx)
            
            assert(sum(idx > obj.frame_count) == 0);
            
            current_time = obj.video_reader.CurrentTime;
            
            n_elements = length(idx);
            
            buffer = zeros(obj.m, obj.n, obj.n_channels, n_elements, obj.mat_data_type);
            
            for i = 1:n_elements                
                obj.video_reader.CurrentTime = ...
                    (idx(i) - 1) / obj.video_reader.FrameRate;
                
                buffer(:, :, :, i) = ...
                    cast(obj.video_reader.readFrame(), obj.mat_data_type);
            end
            
            if obj.bin_size > 1
                buffer = cast(convn(buffer, obj.downsampling_kernel, 'same'), obj.mat_data_type);
                if (obj.bin_size >= size(buffer, 4))
                    buffer = buffer(:, :, :, round(size(buffer, 4) / 2));
                else
                    buffer = buffer(:, :, :, ceil(obj.bin_size / 2):obj.bin_size:end);
                end
            end
            
            obj.video_reader.CurrentTime = current_time;
        end
        
        function result = has_batch(obj)
            result = obj.current_frame < obj.frame_count;
        end
        
        function left = batches_left(obj)
            left = ceil((obj.frame_count - obj.current_frame) / ...
                (obj.buffer_size * obj.bin_size));
        end
        
        function close(obj)
            close(obj.video_reader);
        end
    end
end

