% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

classdef IMG_file_reader
    %IMG_FILE_READER reads frames from a folder with images
    
    properties(GetAccess = public, Constant)
        datatype = 'IMG';
    end
    
    properties
        current_frame = 0;
    end
    
    methods
        function obj = IMG_file_reader(input_folder, buffer_size, bin_size)
            
            error('Implementation of an IMG file reader planned for the future!');
            
            assert(isfolder(input_folder));
            
            if (nargin < 3)
                bin_size = 1;
            end
            if (nargin < 2)
                buffer_size = 500;
            end
            
            obj.downsampling_kernel = 1/bin_size * ...
                ones(1, 1, 1, bin_size);
            
            obj.buffer_size = buffer_size;
            obj.bin_size = bin_size;
            
        end

        function reset(obj)
            obj.current_frame = 0;
        end
        
        function buffer = read_batch(obj)
            buffer = [];
        end
        
        function buffer = read_frames(obj, idx)
            
            buffer = [];
        end
        
        function result = has_batch(obj)
            result = false;
        end
        
        function left = batches_left(obj)
            left = 0;
        end        
    end
end

