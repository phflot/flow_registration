% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

% class for reading blocks from a mat file with binning
classdef MAT_file_reader < Video_file_reader & DS_file_reader
    
    properties (Access = private)
        input_file;
        mat_file;
    end
    
    properties(GetAccess = public, Constant)
        datatype = 'MAT';
    end
    
    properties(Access = ?DS_file_reader, Constant)
        %% add a new format here:
        known_datasetnames = {'ch*_reg', 'buffer*'};
        %% add the corresponding data representation here (currently only 3D)
        known_dims = {[1, 2, 3], [1, 2, 3]};
    end
    
    properties
        current_frame = 0;
    end
    
    methods
        function obj = MAT_file_reader(input_file, ...
                buffer_size, bin_size, varargin)
            
            obj = obj@DS_file_reader([1, 2, 3], varargin{:});
            
            obj.input_file = input_file;
            
            if nargin > 1 && ~isempty(buffer_size)
                obj.buffer_size = buffer_size;
            end
            if nargin > 2 && ~isempty(bin_size)
                obj.bin_size = bin_size;
            end
            
            % previously in the open function:
            obj.mat_file = matfile(obj.input_file);
            found_ds_names = fieldnames(obj.mat_file);
            found_ds_names(strcmp(found_ds_names, 'Properties')) = [];
            found_ds_names = sort(found_ds_names);
%             % testing all variables, looking for the most variables
%             % that share the largest dimensions
%             % heuristic assumes that the spatial dimensions are
%             if nargin < 3
% 
%             end
                        
%             videoSize = size(obj.mat_file.(obj.variable_names{1}));
%             
%             if (length(variableSize) == 3)
%                 [obj.m, obj.n, obj.frame_count] = ...
%                     size(obj.mat_file.(obj.variable_names{1}));
%                 obj.n_channels = 1;
%             else
%                 [obj.m, obj.n, obj.n_channels, obj.frame_count] = ...
%                     size(obj.mat_file.(obj.variable_names{1}));
%             end
%             
%             for i = 2:length(obj.variable_names)
%                 assert(sum(videoSize ~= size(obj.mat_file.(obj.variable_names{1}))) == 0);
%             end
                        % identifying the datasets:
                        

            % finding all 3D datasets:
            idx_3D = false(1, length(found_ds_names));
            dim_sum = zeros(1, length(found_ds_names));
            for i = 1:length(found_ds_names)
                datasize = size(obj.mat_file.(found_ds_names{i}));
                dim_sum(i) = sum(datasize);
                idx_3D(i) = length(datasize) == 3;
            end
            if ~any(idx_3D)
                error('no 3D dataset found in file %s', input_file);
            end

            % preparing the supplied or default dataset names:
            obj.process_datasets(found_ds_names, idx_3D, dim_sum);
            
            % verifying the detected datasets:
            if length(obj.dataset_names) > length(found_ds_names)
                error('Too many dataset names specified!');
            end

            if (isstring(obj.dataset_names) || ...
                    ischar(obj.dataset_names))
                obj.dataset_names = char(obj.dataset_names);
                obj.dataset_names = {obj.dataset_names};
            end

            if length(obj.dataset_names) > length(found_ds_names)
                error('Too many dataset names specified!');
            end

            obj.n_channels = length(obj.dataset_names);
            for i = 1:obj.n_channels

                tmp_name = obj.dataset_names{i};
                if ~any(strcmp(found_ds_names, tmp_name))
                    error('Dataset %s not found in mat file %s', ...
                        tmp_name, input_file);
                end
            end
            
            for i = 1:obj.n_channels
                tmp_name = obj.dataset_names{i};
                datasize = size(obj.mat_file.(tmp_name));
                if length(datasize) ~= 3
                    error('Number of dimensions for dataset %s is %i, 3D input expected!', ...
                        tmp_name, length(datasize));
                end
                if i == 1
                    obj.m = datasize(obj.dimension_ordering(1));
                    obj.n = datasize(obj.dimension_ordering(2));
                    obj.frame_count = datasize(obj.dimension_ordering(3));
                    obj.mat_data_type = class(obj.mat_file.(tmp_name));
                elseif obj.m ~= datasize(obj.dimension_ordering(1)) || ...
                        obj.n ~= datasize(obj.dimension_ordering(2)) || ...
                        obj.frame_count ~= datasize(obj.dimension_ordering(3)) || ...
                        ~strcmp(class(obj.mat_file.(tmp_name)), obj.mat_data_type)
                    error('Dimensions or type for %s do not match!', ...
                        tmp_name);
                end
            end
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
            data_end = [obj.m, obj.n, obj.current_frame + n_elem_left];

            for j = 1:obj.n_channels
                tmp = permute(obj.mat_file.(obj.dataset_names{j})...
                    (data_start(obj.dimension_ordering(1)):data_end(obj.dimension_ordering(1)), ...
                    data_start(obj.dimension_ordering(2)):data_end(obj.dimension_ordering(2)), ...
                    data_start(obj.dimension_ordering(3)):data_end(obj.dimension_ordering(3))), ...
                    obj.dimension_ordering);
                buffer(:, :, j, :) = tmp;
            end

            obj.current_frame = obj.current_frame + n_elem_left;
            
            if obj.bin_size > 1
                buffer = cast(convn(buffer, obj.downsampling_kernel, 'same'), obj.mat_data_type);
                buffer = buffer(:, :, :, ceil(obj.bin_size / 2):obj.bin_size:end);
            end
        end
                
        function buffer = read_frames(obj, idx)
            assert(sum(idx > obj.frame_count) == 0);
            
            n_elements = length(idx);
            
            buffer = zeros(obj.m, obj.n, obj.n_channels, n_elements, ...
                obj.mat_data_type);
            
            idx_full{1} = 1:obj.m;
            idx_full{2} = 1:obj.n;
            idx_full{3} = idx;

            for j = 1:obj.n_channels
                tmp = permute(obj.mat_file.(obj.dataset_names{j})...
                    (idx_full{obj.dimension_ordering(1)}, ...
                    idx_full{obj.dimension_ordering(2)}, ...
                    idx_full{obj.dimension_ordering(3)}), ...
                    obj.dimension_ordering);
                buffer(:, :, j, :) = tmp;
            end
            
            if obj.bin_size > 1
                buffer = cast(convn(buffer, obj.downsampling_kernel, 'same'), obj.mat_data_type);
                if (obj.bin_size >= size(buffer, 4))
                    buffer = buffer(:, :, :, round(size(buffer, 4) / 2));
                else
                    buffer = buffer(:, :, :, ceil(obj.bin_size / 2):obj.bin_size:end);
                end
            end
        end
        
        function result = has_batch(obj)
            result = obj.current_frame < obj.frame_count;
        end
        
        function left = batches_left(obj)
            left = ceil((obj.frame_count - obj.current_frame) / ...
                (obj.buffer_size * obj.bin_size));
        end
        
        function reset(obj)
            obj.current_frame = 0;
        end
    end
    
    methods (Access = private)
        
        function [channel_names, dim_ordering] = get_dim_ordering(obj)

            vars = fieldnames(obj.mat_file);
            found_matrix = false;
            for i = 1:length(vars)
                if ~isnumeric(obj.mat_file.(vars{i}))
                    continue;
                end
                n_dims = ndims(obj.mat_file.(vars{i}));
                dims = size(obj.mat_file.(vars{i}));

                found_matrix = true;

                stds = zeros(1, n_dims);
                for j = 1:n_dims
                    stds(j) = mean(squeeze(std(tmp, 0, j)), [1, 2]);
                end
                % find the two most similar standard deviations:
                [sorted, idx_s] = sort(stds);
                sorted_grad = diff(sorted);
                [~, idx] = min(sorted_grad);

                if (idx_s(idx) == 1 || idx_s(idx) == 2)
                    spatial_dims = [1, 2];
                    if (n_dims == 4)
                        channel_dim = 3;
                        temporal_dim = 4;
                    elseif (n_dims == 3)
                        temporal_dim = 3;
                    end
                else                    
                    % test
                    spatial_dims = [idx_s(idx + 1), idx_s(idx)];
                    temporal_dim = setdiff(1:n_dims, spatial_dims);
                    if (n_dims == 4)
                        channel_dim = setdiff(1:n_dims, ...
                            [temporal_dim, spatial_dims]);
                    end
                end

            end
            if ~found_matrix
                error('Could not find data in the specified mat file!');
            end
        end
        
        function [success, dim_ordering] = get_dim_ordering_single(...
                obj, variable)
            
            success = true;
            channel_dim = [];
            n_dims = ndims(obj.mat_file.(variable));
            dims = size(obj.mat_file.(variable));

            switch n_dims
                case 2
                    success = false;
                    return;
                case 3
                    tmp = obj.mat_file.(variable)...
                        (1:min(dims(1), 512), ...
                         1:min(dims(2), 512), ...
                         1:min(dims(3), 512));
                case 4
                    tmp = obj.mat_file.(variable)...
                        (1:min(dims(1), 512), ...
                         1:min(dims(2), 512), ...
                         1:min(dims(3), 512), ...
                         1:min(dims(4), 512));
                otherwise
                    success = false;
                    return;
            end
            
            
            stds = zeros(1, n_dims);
            for j = 1:n_dims
                stds(j) = mean(squeeze(std(tmp, 0, j)), [1, 2]);
            end
            % find the two most similar standard deviations:
            [sorted, idx_s] = sort(stds);
            sorted_grad = diff(sorted);
            [~, idx] = min(sorted_grad);

            if (idx_s(idx) == 1 || idx_s(idx) == 2)
                spatial_dims = [1, 2];
                if (n_dims == 4)
                    channel_dim = 3;
                    temporal_dim = 4;
                elseif (n_dims == 3)
                    temporal_dim = 3;
                end
            else                    
                % test
                spatial_dims = [idx_s(idx + 1), idx_s(idx)];
                temporal_dim = setdiff(1:n_dims, spatial_dims);
                [~, temporal_dim_idx] = max(...
                    size(obj.mat_file.(variable), temporal_dim));
                temporal_dim = temporal_dim(temporal_dim_idx(1));
                if (n_dims == 4)
                    channel_dim = setdiff(1:n_dims, ...
                        [temporal_dim, spatial_dims]);
                end
            end
            
            dim_ordering = [spatial_dim, channel_dim, temporal_dim];
        end
    end
end

