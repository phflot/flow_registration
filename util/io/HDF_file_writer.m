% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

classdef HDF_file_writer < Video_file_writer & DS_file_writer
    %HDF5_FILE_WRITER Class for writing hdf5 files
    
    properties (SetAccess = private, GetAccess = public)
        file_name;
    end
    
    properties (Access = private)
        frame_counter = 1;
    end
    
    methods
        function obj = HDF_file_writer(file_name, varargin)
            
            obj = obj@DS_file_writer([2, 1, 3], varargin{:});
            
            obj.file_name = file_name;
        end
        
        function success = write_frames(obj, frames)
            success = false;
            [m, n, ~, t] = size(frames);
            if isempty(obj.bitdepth)
                obj.init(frames);
                
                datatype = obj.get_mat_data_type(obj.bitdepth, false);
                
                if (isfile(obj.file_name))
                    delete(obj.file_name);
                end
                
                for i = 1:obj.n_channels
                    ds = obj.get_ds_name(i, obj.n_channels);
                    
                    dataspace = [m, n, Inf];
                    chunksize = [m, n, 1];
                    
                    h5create(obj.file_name, ...
                        ['/' ds], ...
                        dataspace(obj.dimension_ordering), ...
                        'ChunkSize', chunksize(obj.dimension_ordering), ...
                        'Datatype', datatype);
                end
            end
            for i = 1:obj.n_channels
                ds = obj.get_ds_name(i, obj.n_channels);
                
                start_idx = [1, 1, obj.frame_counter];
                n_data = [m, n, t];

                h5write(obj.file_name, ['/' ds], ...
                    permute(squeeze(frames(:, :, i, :)), obj.dimension_ordering), ...
                    start_idx(obj.dimension_ordering), ...
                    n_data(obj.dimension_ordering));
            end
            obj.frame_counter = obj.frame_counter + t;
            success = true;
        end

        function close(~)
            
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
end

