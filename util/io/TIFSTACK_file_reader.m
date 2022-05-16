% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

% class for reading blocks from an multipage tif file with binning
classdef TIFSTACK_file_reader < Video_file_reader
    
    properties (Access = private)
        tif;
        deinterleave = 1;
    end
    
    properties(GetAccess = public, Constant)
        datatype = 'TIFSTACK';
    end
    
    properties
        current_frame = 0;
    end
    
    methods
        function obj = TIFSTACK_file_reader(input_file, buffer_size, ...
                bin_size, varargin)
            
            p = inputParser;
            addParameter(p, 'deinterleave', 1, @(x) isnumeric(x) & isscalar(x));
            parse(p, varargin{:});
            obj.deinterleave = p.Results.deinterleave;
            
            warning('off', 'imageio:tiffmexutils:libtiffWarning');
            
            if nargin > 1 && ~isempty(buffer_size)
                obj.buffer_size = buffer_size;
            end
            if nargin > 2 && ~isempty(bin_size)
                obj.bin_size = bin_size;
            end
            
            obj.tif = Tiff(input_file, 'r');
            obj.n_channels = obj.tif.getTag('SamplesPerPixel');
            obj.bitdepth = obj.tif.getTag('BitsPerSample');
            obj.mat_data_type = obj.get_mat_data_type(obj.bitdepth, ...
                obj.tif.getTag('SampleFormat') == Tiff.SampleFormat.Int);
            
            img_struct = imfinfo(input_file);
            obj.frame_count = length(img_struct);
            
            obj.m = obj.tif.getTag('ImageLength');
            % number of columns
            obj.n = obj.tif.getTag('ImageWidth');
            
            if obj.deinterleave > 1
                obj.n_channels = obj.deinterleave;
                obj.frame_count = obj.frame_count / obj.n_channels;
            end
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
            buffer = zeros(obj.m, obj.n, obj.n_channels, ...
                n_elem_left, obj.mat_data_type);
            
            for i = 1:n_elem_left
                obj.current_frame = obj.current_frame + 1;
                current_idx = obj.deinterleave * (obj.current_frame - 1) + 1;
                obj.tif.setDirectory(current_idx);
                
                if obj.deinterleave > 1
                    for j = 1:obj.deinterleave
                        obj.tif.setDirectory(current_idx + j - 1);

                        buffer(:, :, j, i) =...
                            cast(obj.tif.read, obj.mat_data_type);
                    end
                else
                    buffer(:, :, :, i) =...
                        cast(obj.tif.read, obj.mat_data_type);
                end
            end
            
            buffer = obj.bin_buffer(buffer);
        end
        
        function buffer = read_frames(obj, idx)
            
            assert(sum(idx > obj.frame_count) == 0);
            
            n_elements = length(idx);
            
            buffer = zeros(obj.m, obj.n, obj.n_channels, n_elements, ...
                obj.mat_data_type);
            
            for i = 1:n_elements
                current_idx = (obj.deinterleave * (idx(i) - 1)) + 1;
                obj.tif.setDirectory(current_idx);

                if obj.deinterleave > 1
                    for j = 1:obj.deinterleave
                        obj.tif.setDirectory(current_idx + j - 1);

                        buffer(:, :, j, i) =...
                            cast(obj.tif.read, obj.mat_data_type);
                    end
                else
                    buffer(:, :, :, i) =...
                        cast(obj.tif.read, obj.mat_data_type);
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
        
        function close(obj)
            close(obj.tif);
        end
    end
end

