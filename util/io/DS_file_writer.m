% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

classdef DS_file_writer < handle
    % utility class to wrap handling of the dataset name specification
    % for the file writers with multiple channels. 
    % formats are 'ch*_reg' for ch1_reg, ch2_reg, ...
    % cell {ch1_reg, ch2_reg} for the same output or 
    % 'ch1_reg' for single channel fixed dataset names
    
    properties (Access = protected)
        dataset_names;
        dimension_ordering;
    end
    
    methods (Access = protected)
        function obj = DS_file_writer(default_order, varargin)
            p = inputParser;
            addParameter(p, 'dataset_names', [], @(x) isstring(x) || ...
                ischar(x) || iscell(x));
            addParameter(p, 'dimension_ordering', default_order, @(x) ...
                isvector(x) && isnumeric(x) && length(x) == 3 );
            parse(p, varargin{:})
            obj.dataset_names = p.Results.dataset_names;
            obj.dimension_ordering = p.Results.dimension_ordering;
            
            if ~isempty(obj.dataset_names)
                if iscell(obj.dataset_names)
                    for i = 1:length(obj.dataset_names)
                        assert(isstring(obj.dataset_names{i}) || ischar(obj.dataset_names{i}));
                        obj.dataset_names{i} = char(obj.dataset_names{i});
                        if obj.dataset_names{i}(1) == '/'
                            obj.dataset_names{i} = obj.dataset_names{i}(2:end);
                        end
                    end
                else
                    obj.dataset_names = char(obj.dataset_names);
                    if obj.dataset_names(1) == '/'
                        obj.dataset_names = obj.dataset_names(2:end);
                    end
                end
            end
        end
        
        function ds = get_ds_name(obj, channel_id, n_channels)
            
            if ~isempty(obj.dataset_names)
                if iscell(obj.dataset_names)
                    assert(length(obj.dataset_names) == n_channels);
                    ds = obj.dataset_names{channel_id};
                else
                    tmp = split(obj.dataset_names, '*');
                    if n_channels == 1 && ~contains(obj.dataset_names, '*')
                        ds = obj.dataset_names;                        
                    else
                        if length(tmp) == 1
                            ds = [tmp{1} num2str(channel_id)];
                        else
                            ds = [tmp{1} num2str(channel_id) tmp{2}];
                        end
                    end
                end
            else
                ds = ['ch' num2str(channel_id)];
            end
        end
    end
end

