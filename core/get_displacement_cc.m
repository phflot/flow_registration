% Author   : Philipp Flotho
% Copyright 2023 by Philipp Flotho, All rights reserved.

function w = get_displacement_cc( fixed, moving, varargin )
% computes rigid displacements via cross correlation
    [m, n, n_channels] = size(fixed);

    xoffsets = zeros(1, n_channels);
    yoffsets = zeros(1, n_channels);
    
    for k = 1:n_channels
        c = normxcorr2(fixed(:, :, k), moving(:, :, k));
    
        [max_c, imax] = max(abs(c(:)));
        [ypeak, xpeak] = ind2sub(size(c),imax(1));
    
        yoffsets(k) = ypeak - m;
        xoffsets(k) = xpeak - n;
    end
    
    w = ones(m, n, 2);
    u = mean(xoffsets);
    v = mean(yoffsets);
    
    w = w .* reshape([u, v], [1, 1, 2]);
end

