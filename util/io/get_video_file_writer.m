% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function video_writer = get_video_file_writer(filename, ...
    output_format, varargin)
%GET_VIDEO_FILE_WRITER constructs video file writers

    if (nargin < 3)
        ds_name = [];
    end

    switch output_format
        case 'TIFF'
            video_writer = TIFSTACK_file_writer(filename, varargin{:});
        case 'MAT'
            video_writer = MAT_file_writer(filename, varargin{:});
        case 'HDF5'
            video_writer = HDF_file_writer(filename, varargin{:});
        case 'MULTIFILE_TIFF'
            video_writer = MULTIFILE_file_writer(filename, 'TIFF', varargin{:});
        case 'MULTIFILE_MAT'
            video_writer = MULTIFILE_file_writer(filename, 'MAT', varargin{:});
        case 'MULTIFILE_HDF5'
            video_writer = MULTIFILE_file_writer(filename, 'HDF5', varargin{:});
        case 'CAIMAN_HDF5'
            % multifile hdf5 reader with /mov as dataset to be readable
            % by bigread2
            video_writer = MULTIFILE_file_writer(filename, 'HDF5', ...
                'dataset_names', '/mov');
        case 'BEGONIA'
            video_writer = TSERIESH5_file_writer(filename, varargin{:});
        otherwise
            error('No valid output file in the options set!');
    end
end