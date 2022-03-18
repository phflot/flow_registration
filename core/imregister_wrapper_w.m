% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function [ registered ] = imregister_wrapper_w( f2, w, f1, interpolation_method )

    if (nargin < 3)
        f1 = f2;
    end

    if (nargin < 4)
        interpolation_method = 'cubic';
    end

    if verLessThan('matlab', '9.10')
        size_f2 = size(f2);
        size_f1 = size(f1);
        size_w = size(w);        
        
        assert(sum(size_f2(1:2) == size_w(1:2)) == 2 && ...
               sum(size_f1(1:2) == size_w(1:2)) == 2 && ...
               sum(size(f1) == size(f1)) == ndims(f1), ...
               'imregister sizes do not match, size \n f1 = (%i, %i) f2 = (%i, %i) u = (%i, %i)', ...
               size(f1, 1), size(f1, 2), size(f2, 1), size(f2, 2), ...
               size(w, 1), size(w, 2));
    else
        assert(sum(size(f2, [1, 2]) == size(w, [1, 2])) == 2 && ...
               sum(size(f1, [1, 2]) == size(w, [1, 2])) == 2 && ...
               sum(size(f2) == size(f1)) == ndims(f1), ...
               'imregister sizes do not match, size \n f1 = (%i, %i) f2 = (%i, %i) u = (%i, %i)', ...
               size(f1, 1), size(f1, 2), size(f2, 1), size(f2, 2), ...
               size(w, 1), size(w, 2));
    end
    
    registered = zeros(size(f2), 'double');
    for i = 1:size(f2, 3)
        registered(:, :, i) = imwarp(double(f2(:, :, i)), w, interpolation_method, ...
            'FillValues', nan);
    end
    
    idx = isnan(registered);
    registered(idx) = f1(idx);
%     registered = cast(registered, class(f2));
end

