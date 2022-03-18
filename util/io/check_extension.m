% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function [result, ext] = check_extension(ext)
    supported_extensions = {'.MDF', '.tif', '.tiff', ...
        '.mat'};
    ext_map = {'MDF', 'TIFSTACK', 'TIFSTACK', 'MAT'};

    result = false;
    if isempty(ext)
        return;
    end
    for i = 1:length(supported_extensions)
        result = ...
            result | strcmpi(ext, supported_extensions{i});
        if result
            idx = i;
            ext = ext_map{idx};
            return;
        end
    end
end

