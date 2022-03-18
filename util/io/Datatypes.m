% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

classdef Datatypes
    %REGISTRATION_DATATYPES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, Constant)
        NONE = 'NONE';
        MAT = 'MAT';
        HDF5 = 'HDF5';
        MDF = 'MDF';
        TIFSTACK = 'TIFSTACK';
        IMG = 'IMG';
        VIDEO = 'VIDEO'
        MULTICHANNEL = 'MULTICHANNEL';
    end
end

