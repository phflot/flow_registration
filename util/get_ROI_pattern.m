% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function BW = get_ROI_pattern(rois, m, n)

    BW = false(m, n);
    
    for i = 1:length(rois)
        BW = BW | roipoly(zeros(m, n), ...
            rois{i}.mnCoordinates(:, 1) + 0.5, rois{i}.mnCoordinates(:, 2) + 0.5);
    end
end