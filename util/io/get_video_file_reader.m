% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function video_file_reader = get_video_file_reader(input_file, ...
    buffer_size, bin_size, varargin)
%GET_VIDEO_WRITER Facory function that returns a file reader

    datatype = check_input(input_file);
    
    if (nargin < 3 || isempty(bin_size))
        bin_size = 1;
    end
    if (nargin < 2 || isempty(buffer_size))
        buffer_size = 500;
    end
    
    datatypes = Datatypes;
    switch datatype
        case datatypes.TIFSTACK
            video_file_reader = ...
                TIFSTACK_file_reader(input_file, ...
                buffer_size, bin_size, varargin{:});
        case datatypes.MDF
            video_file_reader = MDF_file_reader(input_file, ...
                buffer_size, bin_size, varargin{:});
        case datatypes.VIDEO
            video_file_reader = AVI_file_reader(input_file, ...
                buffer_size, bin_size, varargin{:});
        case datatypes.MAT
            video_file_reader = MAT_file_reader(input_file, ...
                buffer_size, bin_size, varargin{:});
        case datatypes.MULTICHANNEL
            video_file_reader = MULTICHANNEL_file_reader(...
                input_file, buffer_size, bin_size, varargin{:});
        case datatypes.IMG
            video_file_reader = IMG_file_reader(input_file, ...
                buffer_size, bin_size, varargin{:});
        case datatypes.HDF5
            % todo: factory method here:
            video_file_reader = get_HDF_file_reader(input_file, ...
                buffer_size, bin_size, varargin{:});
        otherwise
            error('File format currently not supported');
    end
end

function hdf_reader = get_HDF_file_reader(input_file, ...
                buffer_size, bin_size, varargin)
    try 
        hdf_reader = HDF_file_reader(input_file, ...
            buffer_size, bin_size, varargin{:});
        return;
    catch
    end
    hdf_reader = HDF_file_reader_4D(input_file, ...
            buffer_size, bin_size, varargin{:});
end

function datatype = check_input(x)
    datatypes = Datatypes;
    datatype = datatypes.NONE;
    
    if (isempty(x))
        datatype = datatypes.NONE;
        return;
    end

    if iscell(x) && (isstring(x{1}) || ischar(x{1}))
%         is_multichannel = true;
        for i = 1:length(x)
            check_input(x{i});
%             is_multichannel = is_multichannel & check_input(x{i});
        end
        datatype = datatypes.MULTICHANNEL;
        return;
    end

    if (exist(x, 'file') == 2)
        [~, ~, ext] = fileparts(x);

        [result, type] = check_extension(ext);
        if ~result
            result = H5F.is_hdf5(x);
            if result
                datatype = datatypes.HDF5;
            else
                try
                    VideoReader(x);
                    datatype = datatypes.VIDEO;
                    result = true;
                catch
                end
            end
        else
            datatype = type;
        end

        if (~result)
            error('File format currently not supported!');
        end
    end
    if isfolder(x)
        % checking for IMG files:
        dir_content = dir(x);
        % checking if at least one image present:
        IMG_ready = false;
        counter = 1;
        while (counter < length(dir_content))
            try
                imfinfo(fullfile(...
                    dir_content(counter).folder, ...
                    dir_content(counter).name));
                IMG_ready = true;
                break;
            catch
            end
            counter = counter + 1;
        end
        if (~IMG_ready)
            warning('Folder does not contain images!');
        end
        datatype = datatypes.IMG;
    end
    if ~isfolder(x) && ~(exist(x, 'file') == 2)
        error('Input File or Folder %s does not exist!', x);
    end
end