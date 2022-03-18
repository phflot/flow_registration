% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function y = get_time_course(vid, rois)
    [m, n, t] = size(vid);

    y = zeros(t, length(rois));
    
    for i = 1:length(rois)
        BW = roipoly(zeros(m, n), ...
            rois{i}.mnCoordinates(:, 1) + 0.5, rois{i}.mnCoordinates(:, 2) + 0.5);
        y(:, i) = double(sum(sum(vid .* BW, 1), 2)) ./ sum(BW(:));
    end
end