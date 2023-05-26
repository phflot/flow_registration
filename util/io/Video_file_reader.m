% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

classdef Video_file_reader < handle
    
    properties(GetAccess = public, SetAccess = protected)
        frame_count;
        n_channels;
        bitdepth;
        mat_data_type;
        input_file_name;
    end
    
    properties(Access = protected)
        downsampling_kernel;
        m;
        n;
    end
    
    properties(Access = public)
        buffer_size = 500;
        bin_size = 1;
    end
    
    properties(Abstract = true, GetAccess = public, Constant)
        datatype;
    end
    
    properties (Abstract = true)
        current_frame;
    end
    
    methods (Abstract = true)
        reset(obj);
        buffers = read_batch(obj);
        buffers = read_frames(obj, idx);
        result = has_batch(obj);
        left = batches_left(obj);
    end
    
    methods
        function set.buffer_size(obj, x)
            obj.buffer_size = x;
            obj.reset;
        end
        
        function set.bin_size(obj, x)
            obj.bin_size = x;
            obj.downsampling_kernel = 1/obj.bin_size * ...
                ones(1, 1, 1, obj.bin_size);
            obj.reset;
        end
        
        function width = get_width(obj)
            width = obj.n;
        end
        
        function height = get_height(obj)
            height = obj.m;
        end
    end
    
    methods (Static, Access = protected)
        function mat_data_type = get_mat_data_type(bitdepth, signed)
            if nargin < 2
                signed = false;
            end
            switch bitdepth
                case 8
                    if signed
                        mat_data_type = 'int8';
                    else
                        mat_data_type = 'uint8';
                    end
                case {10, 12, 16}
                    if signed
                        mat_data_type = 'int16';
                    else
                        mat_data_type = 'uint16';
                    end
                otherwise
                    mat_data_type = 'double';
            end
        end
    end

    methods (Access = protected)
        function buffer = bin_buffer(obj, buffer)
            if obj.bin_size > 1
                if size(buffer, 4) <= obj.bin_size
                    buffer = mean(buffer, 4);
                else
                    buffer = cast(convn(buffer, obj.downsampling_kernel, 'same'), obj.mat_data_type);
                    buffer = buffer(:, :, :, ceil(obj.bin_size / 2):obj.bin_size:end);
                end
            end
        end
    end
end

