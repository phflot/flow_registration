% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function [ c_comp ] = compensate_sequence_uv( c, c_ref, w, interpolation_method )
%COMPENSATE_SEQUENCE_UV compensates a sequence w.r.t. c_ref

    c_comp = zeros(size(c), class(c));

    if nargin < 4
        interpolation_method = 'cubic';
    end
    
    for i = 1:size(c, 4)
        c_comp(:, :, :, i) = imregister_wrapper_w(c(:, :, :, i), ...
            w(:, :, :, i), c_ref, interpolation_method);
    end
end

