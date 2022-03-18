% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

classdef TSERIESH5_file_writer < Video_file_writer
    %TSERIESH5_FILE_WRITER Class for writing tseries h5 files that are
    %compatible with Begonia
    
    properties (SetAccess = private, GetAccess = public)
        file_name;
    end
    
    properties (Access = private)
        frame_counter = 1;
        parameters;
    end
    
    methods
        function obj = TSERIESH5_file_writer(file_name, varargin)
            [path, filename, ~] = fileparts(file_name);
            obj.parameters = varargin;
            obj.file_name = fullfile(path, [filename '.h5']);
        end
        
        function success = write_frames(obj, frames)
            success = false;
            [m, n, n_channels, t] = size(frames);
            if isempty(obj.bitdepth)
                obj.init(frames);
                
                datatype = obj.get_mat_data_type(obj.bitdepth, false);
                
                if (isfile(obj.file_name))
                    delete(obj.file_name);
                end
                
                dataspace = [m, n, n_channels, Inf];
                chunksize = [m, n, n_channels, 1];

                h5create(obj.file_name, ...
                    '/recording', ...
                    dataspace, ...
                    'ChunkSize', chunksize, ...
                    'Datatype', datatype);
                %set_h5_tseries_attributes(obj.file_name, ...
                %    obj.parameters);
            end
                
            start_idx = [1, 1, 1, obj.frame_counter];
            n_data = [m, n, n_channels, t];

            h5write(obj.file_name, '/recording', ...
                frames(end:-1:1, :, :, :), start_idx, n_data);
            
            obj.frame_counter = obj.frame_counter + t;
            
            tmp = obj.parameters;
            tmp{end + 1} = 'frame_count';
            tmp{end + 1} = obj.frame_counter - 1;
            tmp{end + 1} = 'slices';
            tmp{end + 1} = obj.frame_counter - 1;
            
            % final update of the metadata
            set_h5_tseries_attributes(obj.file_name, ...
                tmp{:});
            
            success = true;
        end

        function close(obj)
            
            tmp = obj.parameters;
            tmp{end + 1} = 'frame_count';
            tmp{end + 1} = obj.frame_counter - 1;
            tmp{end + 1} = 'slices';
            tmp{end + 1} = obj.frame_counter - 1;
            
            % final update of the metadata
            set_h5_tseries_attributes(obj.file_name, ...
                tmp{:});
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

