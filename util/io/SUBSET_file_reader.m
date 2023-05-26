% Author   : Philipp Flotho
% Copyright 2022 by Philipp Flotho, All rights reserved.


classdef SUBSET_file_reader < Video_file_reader
    %SUBSET_FILE_READER class that outputs a subset of a video reader
    
    properties(SetAccess = protected, GetAccess = public)
        video_file_reader;
        idx;
    end

    properties
        current_frame = 0;
    end
    
    properties(GetAccess = public, Constant)
        datatype = 'SUBSET'
    end

    methods
        function obj = SUBSET_file_reader(video_file_reader, idx)
            obj.video_file_reader = video_file_reader;
            obj.buffer_size = video_file_reader.buffer_size;
            obj.bin_size = video_file_reader.bin_size;
            obj.downsampling_kernel = video_file_reader.downsampling_kernel;
            obj.idx = idx;
            obj.m = video_file_reader.get_height();
            obj.n = video_file_reader.get_width();
            obj.n_channels = video_file_reader.n_channels;
            obj.bitdepth = video_file_reader.bitdepth;
            obj.mat_data_type = video_file_reader.mat_data_type;
            obj.frame_count = length(idx);
            
            obj.input_file_name = video_file_reader.input_file_name;
        end
        
        function buffer = read_batch(obj)
            obj.video_file_reader.bin_size = 1;

            if obj.current_frame > obj.frame_count
                buffer = [];
                return;
            end
            
            n_elem_left = min(obj.buffer_size * obj.bin_size, ...
                obj.frame_count - obj.current_frame);
            buffer = zeros(obj.m, obj.n, obj.n_channels, ...
                n_elem_left, obj.mat_data_type);

            for i = 1:n_elem_left
                buffer(:, :, :, i) = obj.video_file_reader.read_frames(obj.idx(i));
            end

            obj.video_file_reader.bin_size = obj.bin_size;

            obj.current_frame = obj.current_frame + n_elem_left;
            
            buffer = obj.bin_buffer(buffer);
        end

        function buffer = read_frames(obj, idx)
            obj.video_file_reader.bin_size = 1;
            buffer = obj.video_file_reader.read_frames(obj.idx(idx));
            obj.video_file_reader.bin_size = obj.bin_size;
        end

        function left = batches_left(obj)
            left = ceil((obj.frame_count - obj.current_frame) / ...
                (obj.buffer_size * obj.bin_size));
        end

        function result = has_batch(obj)
            result = obj.current_frame < obj.frame_count;
        end

        function reset(obj)
            obj.current_frame = 0;
        end
    end
end

