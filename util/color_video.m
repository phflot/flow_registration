% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function color_video = color_video(vid, cmap_name)
    if nargin < 2
        cmap = colormap('parula');
    else
        cmap = colormap(cmap_name);
    end
    
    [~, idx] = imquantize(vid, linspace(min(vid(:)), max(vid(:)), ...
        size(cmap, 1)));
    
    [m, n, t] = size(vid);
    color_video = zeros(m, n, 3, t, class(vid));
    
    for i = 1:size(vid, 3)
        color_video(:, :, :, i) = ind2rgb(idx(:, :, i), cmap);
    end
end

