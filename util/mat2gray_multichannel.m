% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function c_out = mat2gray_multichannel(c, ref)

    c_out = zeros(size(c), 'double');
    
    if nargin < 2
        for i = 1:size(c, 3)
            c_out(:, :, i, :) = mat2gray(c(:, :, i, :));
        end
    else
        for i = 1:size(c, 3)
            ch_ref = ref(:, :, i);
            min_ref = min(double(ch_ref(:)));
            max_ref = max(double(ch_ref(:)));
            
            tmp = double(c(:, :, i, :));
            tmp = (tmp - min_ref) / (max_ref - min_ref);
            
            c_out(:, :, i, :) = tmp;
        end
    end
end

