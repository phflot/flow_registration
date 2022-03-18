% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function img_filt = imgaussfilt3_multichannel(img, options, offset)

    if nargin < 3
        offset = [0, 0, 0];
    end
    
    img_filt = zeros(size(img), class(img));
    for i = 1:size(img, 3)
        img_filt(:, :, i, :) = imgaussfilt3(squeeze(img(:, :, i, :)), ...
            options.get_sigma_at(i) + offset);
    end
end

