% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function convert2tseriesh5(input_file, output_file, varargin)
%CONVERT_TO_TSERIESH5 takes an input path or file reader and an output file
% and converts it to a Begonia compatible h5 file

    if(~isa(input_file, 'Video_file_reader'))
        file_reader = get_video_file_reader(input_file, 100, 1);
    else
        file_reader = input_file;
    end
    if (isa(file_reader, 'MDF_file_reader'))
        varargin{end+1} = 'mdf_reference';
        varargin{end+1} = file_reader;
    end
    
    file_writer = get_video_file_writer(output_file, 'BEGONIA', ...
        varargin{:});
    counter = 1;
    while file_reader.has_batch
        file_writer.write_frames(file_reader.read_batch());
        fprintf('Written batch %i from %i\n', counter, ...
            counter + file_reader.batches_left);
        counter = counter + 1;
    end
    delete(file_writer);
end

