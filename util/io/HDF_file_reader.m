% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

% class for reading blocks from a hdf5 file with binning
classdef HDF_file_reader < Video_file_reader & DS_file_reader
    %HDF_FILE_READER reads frames from an hdf file
    
    properties
        input_file;
        mat_file;
    end
    
    properties(GetAccess = public, Constant)
        datatype = 'HDF5'
    end
    
    properties(Access = ?DS_file_reader, Constant)
        %% add a new format here:
        known_datasetnames = {'mov', 'ch*', 'ch*_reg'};
        %% add the corresponding data representation here (currently only 3D)
        known_dims = {[2, 1, 3], [2, 1, 3], [2, 1, 3]};
    end
    
    properties
        current_frame = 0;
    end
    
    methods
        function obj = HDF_file_reader(input_file, ...
                buffer_size, bin_size, varargin)
            %HDF_FILE_READER
            
            obj = obj@DS_file_reader([2, 1, 3], varargin{:});
            
            obj.input_file = input_file;
            
            if nargin > 1 && ~isempty(buffer_size)
                obj.buffer_size = buffer_size;
            end
            if nargin > 2 && ~isempty(bin_size)
                obj.bin_size = bin_size;
            end
            
            hdf5_ds = h5info(input_file);
            hdf5_ds = hdf5_ds.Datasets;
            found_ds_names = cell(1, length(hdf5_ds));
            for i = 1:length(hdf5_ds)
                found_ds_names{i} = hdf5_ds(i).Name;
            end
            found_ds_names = sort(found_ds_names);
            
            % finding all 3D datasets:
            idx_3D = false(1, length(found_ds_names));
            dim_sum = zeros(1, length(found_ds_names));
            for i = 1:length(found_ds_names)
                datasize = h5info(input_file, ['/', ...
                    found_ds_names{i}]);
                datasize = datasize.Dataspace.Size;
                dim_sum(i) = sum(datasize);
                idx_3D(i) = length(datasize) == 3;
            end
            if ~any(idx_3D)
                error('no 3D dataset found in file %s', input_file);
            end

            % preparing the supplied or default dataset names:
            obj.process_datasets(found_ds_names, idx_3D, dim_sum);
            
            % verifying the detected datasets:
            if length(obj.dataset_names) > length(hdf5_ds)
                error('Too many dataset names specified!');
            end

            obj.n_channels = length(obj.dataset_names);
            for i = 1:obj.n_channels
                if obj.dataset_names{i}(1) == '/'
                    obj.dataset_names{i} = obj.dataset_names{i}(2:end);
                end

                tmp_name = obj.dataset_names{i};
                if ~any(strcmp(found_ds_names, tmp_name))
                    error('Dataset %s not found in file %s', ...
                        tmp_name, input_file);
                end
            end            
            
            for i = 1:obj.n_channels
                tmp_name = obj.dataset_names{i};
                datasize = h5info(input_file, ['/' tmp_name]);
                datasize = datasize.Dataspace.Size;
                if length(datasize) ~= 3
                    error('Number of dimensions for dataset %s is %i, 3D input expected!', ...
                        tmp_name, length(datasize));
                end
                if i == 1
                    obj.m = datasize(obj.dimension_ordering(1));
                    obj.n = datasize(obj.dimension_ordering(2));
                    obj.frame_count = datasize(obj.dimension_ordering(3));
                    obj.mat_data_type = class(...
                        h5read(input_file, ['/' tmp_name], [1, 1, 1], [1, 1, 1]));
                elseif obj.m ~= datasize(obj.dimension_ordering(1)) || ...
                        obj.n ~= datasize(obj.dimension_ordering(2)) || ...
                        obj.frame_count ~= datasize(obj.dimension_ordering(3)) || ...
                        ~strcmp(class(...
                        h5read(input_file, ['/' tmp_name], [1, 1, 1], [1, 1, 1])), obj.mat_data_type)
                    error('Dimensions or type for %s do not match!', ...
                        tmp_name);
                end
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
            
            data_start = [1, 1, obj.current_frame + 1];
            data_offset = [obj.m, obj.n, n_elem_left];

            for j = 1:obj.n_channels
                ds = ['/' obj.get_ds_name(j, obj.n_channels)];
                tmp = permute(h5read(obj.input_file, ds, ...
                    data_start(obj.dimension_ordering), ...
                    data_offset(obj.dimension_ordering)), ...
                    obj.dimension_ordering);
                buffer(:, :, j, :) = tmp;
            end

            obj.current_frame = obj.current_frame + n_elem_left;
            
            buffer = obj.bin_buffer(buffer);
        end
        
        function buffer = read_frames(obj, idx)
            assert(sum(idx > obj.frame_count) == 0);
            
            n_elements = length(idx);
            
            buffer = zeros(obj.m, obj.n, obj.n_channels, n_elements, ...
                obj.mat_data_type);
            
            data_offset = [obj.m, obj.n, 1];
            for i = 1:n_elements
                data_start = [1, 1, idx(i)];
                for j = 1:obj.n_channels
                    ds = ['/' obj.get_ds_name(j, obj.n_channels)];
                    tmp = permute(h5read(obj.input_file, ds, ...
                        data_start(obj.dimension_ordering), ...
                        data_offset(obj.dimension_ordering)), ...
                        obj.dimension_ordering);
                    buffer(:, :, j, i) = tmp;
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

