% Author   : Philipp Flotho
% Copyright 2022 by Philipp Flotho, All rights reserved.

function energy = get_energy(f, ref, sigma)
    if nargin < 3
        sigma = [2, 2, 0.0001];
    elseif isscalar(sigma)
        sigma = [sigma, sigma, 0.0001];
    end
    [fx, fy] = gradient(imgaussfilt3(double(f), sigma));
    [refx, refy] = gradient(imgaussfilt3(double(ref), sigma));

    energy = abs(fx - refx) + abs(fy - refy);
    
    tmp = sort(energy(:));

    energy = sum(tmp(end-ceil(length(tmp) / 200):end));
end
