% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function d = warpingDepth(eta, levels, m, n)
    min_dim = min(m, n);
    warpingdepth = 0;
    d = warpingdepth;

    for i = 1:levels
        warpingdepth = warpingdepth + 1;
        min_dim = min_dim * eta;
        if (round(min_dim) < 10	)
            break;
        end
        d = warpingdepth;
    end
end