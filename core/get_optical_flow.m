% Author   : Philipp Flotho
% Copyright 2023 by Philipp Flotho, All rights reserved.

function w = get_optical_flow(c, varargin)

    m = size(c, 1);
    n = size(c, 2);
    t = size(c, 4);

    w = zeros(m, n, 2, t, 'double');
    w_tmp = zeros(m, n, 2, 'double');

    for i = 1:t-1
        w_tmp = get_displacement(c(:, :, :, i), ...
            c(:, :, :, i+1), 'w_init', w_tmp,  varargin{:});

        w(:, :, :, i) = w_tmp;
    end
end

