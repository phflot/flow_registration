% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

% class for reading blocks from a raw avi file with binning and debayering
classdef RAWAVI_file_reader < Video_file_reader
    
    properties (Access = private)
        video_reader;
        options;
    end
    
    properties (SetAccess = private, GetAccess = public)
        frame_count;
        buffer_size;
        bin_size;
        n_channels = 3;
        m;
        n;
    end
    
    properties
        current_frame = 0;
    end
    
    methods
        function obj = AVI_file_reader(options)
            assert(isa(options, 'OF_options'));
            
            obj.buffer_size = options.buffer_size;
            obj.bin_size = options.bin_size;
            obj.options = options;
        end
        
        function open(obj)
            obj.video_reader = VideoReader(obj.options.input_file);
            obj.frame_count = obj.video_reader.NumberOfFrames;
            obj.video_reader = VideoReader(obj.options.input_file);
            
            obj.m = obj.video_reader.Height;
            obj.n = obj.video_reader.Width;
        end
        
        function reset(obj)
            obj.current_frame = 0;
        end
        
        function buffer1 = read_batch(obj)
            if obj.current_frame >= obj.frame_count
                buffer1 = [];
                return;
            end
            
            n_elem_left = min(obj.buffer_size * obj.bin_size, ...
                obj.frame_count - obj.current_frame);
            buffer1 = uint8(zeros(obj.m, obj.n, 3, n_elem_left));
            
            for i = 1:n_elem_left
                obj.current_frame = obj.current_frame + 1;
                
                buffer1(:, :, :, i) =...
                    demosaic(obj.video_reader.readFrame, 'rggb');
            end
            
            buffer1 = obj.bin_buffer(buffer1);

        end
        
        function buffer1 = read_frames(obj, idx)
            
            assert(sum(idx >= obj.frame_count) == 0);
            
            n_elements = length(idx);
            
            buffer1 = uint8(zeros(obj.m, obj.n, 3, n_elements));
            
            for i = 1:n_elements
                obj.current_frame = obj.current_frame + 1;
                
                obj.video_reader.CurrentTime = ...
                    (idx(i) - 1) / obj.video_reader.FrameRate;
                
                buffer1(:, :, :, i) =...
                demosaic(obj.video_reader.readFrame(), 'rggb');
            end
            
            buffer1 = obj.bin_buffer(buffer1);
        end
        
        function result = has_batch(obj)
            result = obj.current_frame < obj.frame_count;
        end
        
        function left = batches_left(obj)
            left = ceil((obj.frame_count - obj.current_frame) / ...
                (obj.buffer_size * obj.bin_size));
        end
    end
end

