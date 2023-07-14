% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function [c_reg, w, idx] = compensate_sequence( c, c_ref, c_raw, c_ref_raw, varargin )
%COMPENSATE_SEQUENCE compensates a sequence w.r.t. c_ref

    if ndims(c) == 3
        c = reshape(c, size(c, 1), size(c, 2), 1, size(c, 3));
        c_raw = reshape(c_raw, size(c_raw, 1), size(c_raw, 2), 1, size(c_raw, 3));
    end

    w = get_displacements(c, c_ref, varargin{:});
    
    [c_reg, idx] = compensate_sequence_uv(c_raw, c_ref_raw, w);
end