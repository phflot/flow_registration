% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

classdef MULTIFILE_file_writer < Video_file_writer
    %MULTIFILE_FILE_WRITER file writer that writes one file per channel
    
    properties (Access = private)
        file_writers = {};
        folder;
        file_type = 'TIFF';
    end
    
    properties (SetAccess = private, GetAccess = public)
        file_name;
    end
    
    properties (Access = private)
        writer_parameters = {};
    end
    
    methods
        function obj = MULTIFILE_file_writer(filename, file_type, varargin)
            [filepath, name, ext] = fileparts(filename);
            if isempty(ext)
                obj.folder = fullfile(filepath, name);
                obj.file_name = 'compensated';
            else
                obj.file_name = name;
                obj.folder = filepath;
            end
            
            obj.file_type = file_type;
            
            if nargin > 2
                obj.writer_parameters = varargin;
            end
        end
        
        function success = write_frames(obj, frames)
            success = false;
            if isempty(obj.bitdepth)
                obj.init(frames);
                
                for i = 1:obj.n_channels
                    obj.file_writers{i} = get_video_file_writer(...
                        fullfile(obj.folder, [obj.file_name '_ch' ...
                        num2str(i) '.' obj.file_type]), ...
                        obj.file_type, obj.writer_parameters{:});
                end
            end
            for i = 1:obj.n_channels
                obj.file_writers{i}.write_frames(frames(:, :, i, :));
            end
        end
        
        function close(obj)
            for i = 1:length(obj.file_writers)
                obj.file_writers{i}.close();
            end
        end
    end
end

