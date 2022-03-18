% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function [ registered ] = imregister_wrapper( f2, u, v, f1 )

    if (nargin < 4)
        f1 = f2;
    end
     
    w = zeros([size(u) 2], 'double');
    w(:, :, 1) = u;
    w(:, :, 2) = v;
    
    registered = imregister_wrapper_w(f2, w, f1);
end