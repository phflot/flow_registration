% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

classdef MULTICHANNEL_file_reader < Video_file_reader
    %MULTICHANNEL_FILE_READER Generic multichannel reader, that reads and
    %bins blocks from multiple different video_file_readers and returns
    %arrays with the highest bitrate of the different file readers. 
    
    properties (Access = private)
        filereaders;
        options;
        channel_idx = [1, 2];
    end
    
    properties(GetAccess = public, Constant)
        datatype = 'MULTICHANNEL';
    end
    
    properties
        current_frame = 0;
    end
    
    methods
        function obj = MULTICHANNEL_file_reader(input_files, ...
                buffer_size, bin_size, varargin)
            %MULTICHANNEL_FILE_READER Construct an instance of this class
            %   Detailed explanation goes here

            if (nargin < 3)
                bin_size = 1;
            end
            if (nargin < 2)
                buffer_size = 500;
            end
            
            obj.buffer_size = buffer_size;
            obj.bin_size = bin_size;
            
            obj.filereaders = {};
            obj.n_channels = 0;
            different_bits = false;
            for i = 1:length(input_files)
                filereader = get_video_file_reader(input_files{i}, ...
                    buffer_size, bin_size, varargin{:});
                obj.filereaders{end+1} = filereader;
                obj.n_channels = obj.n_channels + filereader.n_channels;
                if i == 1
                    obj.bitdepth = filereader.bitdepth;
                    obj.mat_data_type = filereader.mat_data_type;
                    obj.m = filereader.m;
                    obj.n = filereader.n;
                    obj.frame_count = filereader.frame_count;
                else
                    if obj.m ~= filereader.m || ...
                            obj.n ~= filereader.n
                        error('Resolution of the inputfiles do not match!');
                    end
                    if obj.frame_count ~= filereader.frame_count
                        error('Frame count of the inputfiles do not match!');
                    end
                    if ~strcmp(obj.mat_data_type, filereader.mat_data_type)
                        obj.mat_data_type = double;
                    end
                    if obj.bitdepth ~= filereader.bitdepth
                        different_bits = true;
                        obj.bitdepth = max(obj.bitdepth, filereader.bitdepth);
                    end
                end
            end

            obj.input_file_name = "";
            for i = 1:length(input_files)
                obj.input_file_name = strcat(obj.input_file_name, ...
                    obj.filereaders.input_file_name);
            end


            if different_bits
                % if bitdepth is different, taking the highest bitrate:
                warning("Different bitrate in the different channels, using %i bits and double output", obj.bitdepth);
            end
        end
        
        function buffer = read_batch(obj)
            idx = 1;
            for i = 1:length(obj.filereaders)
                tmp = obj.filereaders{i}.read_batch;
                if i == 1
                    t = size(tmp, 4);
                    buffer = zeros(obj.m, obj.n, obj.n_channels, t, obj.mat_data_type);
                end
                buffer(:, :, idx:(idx ...
                    + obj.filereaders{i}.n_channels - 1), :) = tmp;
                idx = idx + obj.filereaders{i}.n_channels;
            end
        end
        
        function buffer = read_frames(obj, idx)
            ch_idx = 1;
            for i = 1:length(obj.filereaders)
                tmp = obj.filereaders{i}.read_frames(idx);
                if i == 1
                    t = size(tmp, 4);
                    buffer = zeros(obj.m, obj.n, obj.n_channels, t, obj.mat_data_type);
                end
                buffer(:, :, ch_idx:(ch_idx ...
                    + obj.filereaders{i}.n_channels - 1), :) = tmp;
                ch_idx = ch_idx + obj.filereaders{i}.n_channels;
            end
        end
        
        function result = has_batch(obj)
            result = obj.filereaders{1}.has_batch;
        end
        
        function left = batches_left(obj)
            left = obj.filereaders{1}.batches_left;
        end
        
        function reset(obj)
            for i = 1:length(obj.filereaders)
                obj.filereaders{i}.reset;
            end
        end
        
        function close(obj)
            for i = 1:length(obj.filereaders)
                obj.filereaders{i}.close;
            end
        end
    end
end

