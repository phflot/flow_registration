% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

classdef OF_batchprocessor < handle
    %OF_BATCHPROCESSOR Class that takes filenames or OF_option objects
    % and does batch motion compensation on all elements. It offers the
    % option to use different references, or the same reference for each
    % recording. 
    
    properties
        reference_mode = 'same';
        reference = [];
        registration_jobs = {};
        options;
    end
    
    methods
        function obj = OF_batchprocessor(options, varargin)
            obj.options = options;
            
            p = inputParser;

            addParameter(p, 'reference_mode', obj.reference_mode, ...
                @(x) (isstring(x) || ischar(x)) && ...
                (strcmp(x, 'same') || strcmp(x, 'single')));
            parse(p, varargin{:});
            
            obj.reference_mode = p.Results.reference_mode;
        end
        
        function success = append_job(obj, registration_job)
            success = false;
            obj.registration_jobs{end + 1} = registration_job;
        end
        
        function success = start_all(obj)
            success = true;
            for i = 1:length(obj.registration_jobs)
                try
                    option = copy(obj.options);
                    option.output_path = fullfile(...
                        obj.options.output_path, ['batch' num2str(i)]);
                    option.input_file = obj.registration_jobs{i};
                    
                    fprintf('starting job %i\n', i);
                    
                    if (strcmp(obj.reference_mode, 'same'))
                        if i == 1
                            obj.reference = compensate_recording(option);
                        else
                            compensate_recording(option, obj.reference);
                        end
                    else
                        compensate_recording(option);
                    end
                catch error
                    success = false;
                    warning('Job %i failed with message \n%s', i, ...
                        error.message);
                end
            end
        end
    end
end

