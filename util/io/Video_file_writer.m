% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

classdef Video_file_writer < handle
    %VIDEO_FILE_WRITER Abstract class for the video writers
    
    properties (Access = protected)
        bitdepth = [];
        n_channels;
        m;
        n;
    end
    
    properties (Abstract = true, SetAccess = private, GetAccess = public)
        file_name;
    end
    
    methods (Abstract = true)
        success = write_frames;
        close;
    end
    
    methods (Access = protected)        
        function init(obj, frames) 
            assert(isnumeric(frames) && size(frames(1, 1, 1, 1, :), 5) == 1);
            
            tmp = whos('frames');
            obj.bitdepth = tmp.bytes / length(frames(:)) * 8;
            
            obj.m = size(frames, 1);
            obj.n = size(frames, 2);
            obj.n_channels = size(frames, 3);
        end
    end
    
    methods
        function delete(obj)
            obj.close();
        end
    end
end

