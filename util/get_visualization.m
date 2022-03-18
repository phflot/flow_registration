% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function img = get_visualization(ch1, ch2, ...
    scaling_left, scaling_right, reference_left, reference_right, inverted)

    if (nargin < 7)
        inverted = false;
    end

    if (nargin < 6)
        ch1 = ch1 - min(ch1(:));
        ch1 = ch1 ./ max(ch1(:));

        ch2 = ch2 - min(ch2(:));
        ch2 = ch2 ./ max(ch2(:));
    else
        ch1 = ch1 - min(reference_left(:));
        ch1 = ch1 ./ (max(reference_left(:)) - min(reference_left(:)));

        ch2 = ch2 - min(reference_right(:));
        ch2 = ch2 ./ (max(reference_right(:)) - min(reference_left(:)));
    end

    if (nargin <= 4)
        scaling_left = [1, 0];
        scaling_right = [1, 0];
    end


    ch1 = scaling_left(1) * ch1 - scaling_left(2);
    ch2 = scaling_right(1) * ch2 - scaling_right(2);
    
    img(:, :, 1, :) = ch1;
    img(:, :, 2, :) = (ch1 + ch2) * 0.5;
    img(:, :, 3, :) = ch2;
    
%     img = 3 * img - 0.2;
    img(img < 0) = 0;
    img(img > 1) = 1;
    
    if inverted
        tmp = img(:, :, 3, :);
        img(:, :, 3, :) = img(:, :, 1, :);
        img(:, :, 1, :) = tmp;
    end
end

