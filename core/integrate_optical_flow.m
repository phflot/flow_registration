% Author   : Philipp Flotho
% Copyright 2023 by Philipp Flotho, All rights reserved.

function [ c_comp] = integrate_optical_flow(w, interpolation_method)
%COMPENSATE_SEQUENCE_UV compensates a sequence w.r.t. c_ref
    
    c_comp = zeros(size(w), "double");
    c_ref = w(:, :, :, 1);

    if nargin < 4
        interpolation_method = 'cubic';
    end
    
    for i = 2:size(w, 4)
        w_tmp = w(:, :, :, i-1);
        % w_tmp = imregister_wrapper_w(w(:, :, :, i-1), ...
        %     w(:, :, :, i-1), c_ref, interpolation_method);
        [c_comp(:, :, :, i), ~] = imregister_wrapper_w(w_tmp, ...
            c_comp(:, :, :, i-1), c_ref, interpolation_method);
        c_comp(:, :, :, i) = c_comp(:, :, :, i) + c_comp(:, :, :, i-1);
    end
end
