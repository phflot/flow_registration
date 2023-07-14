% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function [ c_comp, idx ] = compensate_sequence_uv( c, c_ref, w, interpolation_method )
%COMPENSATE_SEQUENCE_UV compensates a sequence w.r.t. c_ref

    c_comp = zeros(size(c), class(c));
    [height, width, ~, n_frames] = size(c);
    idx = false(height, width, 1, n_frames);
    %idx = zeros(size(c, 1), size(c, 2), 1, 'bo')

    if nargin < 4
        interpolation_method = 'cubic';
    end
    
    for i = 1:n_frames
        [c_comp(:, :, :, i), idx(:, :, :, i)] = imregister_wrapper_w(c(:, :, :, i), ...
            w(:, :, :, i), c_ref, interpolation_method);
    end
    idx = uint8(idx) * 255;
end
