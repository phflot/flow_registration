% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

classdef TIFSTACK_file_writer < Video_file_writer
    %TIFSTACK_FILE_WRITER Class for writing tiffstacks
    
    properties (Access = private)
        tif_tagstruct;
        tif_file;
    end
    
    properties (SetAccess = private, GetAccess = public)
        file_name;
    end
    
    methods
        function obj = TIFSTACK_file_writer(file_name, varargin)
            
            warning('off', 'MATLAB:imagesci:Tiff:missingExtraSamples');
            
            obj.tif_tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
            obj.tif_tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            obj.tif_tagstruct.Compression = Tiff.Compression.None;
            
            obj.file_name = file_name;
            obj.tif_file = Tiff(file_name, 'w8');
        end
        
        function success = write_frames(obj, frames)
            success = false;
            if isempty(obj.bitdepth)
                obj.init(frames);
                
                obj.tif_tagstruct.ImageLength = obj.m;
                obj.tif_tagstruct.ImageWidth = obj.n;
                obj.tif_tagstruct.SamplesPerPixel = obj.n_channels;
                obj.tif_tagstruct.BitsPerSample = obj.bitdepth;
                datatype = class(frames);
                
                if  obj.bitdepth == 64
                    obj.tif_tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
                elseif any([8, 16, 32] == obj.bitdepth)
                    if datatype(1) == 'u'
                        obj.tif_tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
                    else
                        obj.tif_tagstruct.SampleFormat = Tiff.SampleFormat.Int;
                    end
                else
                    error('bitdepth / typename not supported by tiff!');
                end
            end
            
            for i = 1:size(frames, 4)
                setTag(obj.tif_file, obj.tif_tagstruct);
                obj.tif_file.write(frames(:, :, :, i));
                obj.tif_file.writeDirectory();
            end
            success = true;
        end
        
        function close(obj)
            obj.tif_file.close();
        end
    end
end

