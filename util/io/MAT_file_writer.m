% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

classdef MAT_file_writer < Video_file_writer & DS_file_writer
    %MAT_FILE_WRITER Class for writing mat files
    
    properties (Access = private)
        mat_file;
        dataset_name;
    end
    
    properties (SetAccess = private, GetAccess = public)
        file_name;
    end
    
    methods
        function obj = MAT_file_writer(file_name, varargin)

            obj = obj@DS_file_writer([1, 2, 3], varargin{:});
            
            obj.mat_file = matfile(file_name, 'Writable', true);
        end
        
        function success = write_frames(obj, frames)
            success = false;
            if isempty(obj.bitdepth)
                obj.init(frames);
                
                for i = 1:obj.n_channels
                    ds = obj.get_ds_name(i, obj.n_channels);
                    obj.mat_file.(ds) = ...
                        permute(squeeze(frames(:, :, i, :)), ...
                        obj.dimension_ordering);
                end
            else
                for i = 1:obj.n_channels
                    ds = obj.get_ds_name(i, obj.n_channels);
                    
                    tmp = permute(squeeze(frames(:, :, i, :)), ...
                        obj.dimension_ordering);
                    
                    time_idx = find(obj.dimension_ordering == 3);
                    time_idx = time_idx(1);
                    current_frame = size(obj.mat_file.(ds), time_idx);
                    
                    idx = {1:obj.m, 1:obj.n, ...
                        current_frame+1:current_frame+size(frames, 4)};
                    
                    obj.mat_file.(ds)( ...
                        idx{obj.dimension_ordering(1)}, ...
                        idx{obj.dimension_ordering(2)}, ...
                        idx{obj.dimension_ordering(3)}) = tmp;
                end
            end
            success = true;
        end
        
        function close(~)
        end
    end
end

