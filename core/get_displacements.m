% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function w = get_displacements(c, c_ref, varargin)

    m = size(c, 1);
    n = size(c, 2);
    t = size(c, 4);

    if t == 1
        w = get_displacement(c_ref, ...
            c, varargin{:});
    else
        w = zeros(m, n, 2, t, 'double');

        parfor i = 1:t
            w_tmp = get_displacement(c_ref, ...
                c(:, :, :, i), varargin{:});
    
            w(:, :, :, i) = w_tmp;
        end
    end
end

