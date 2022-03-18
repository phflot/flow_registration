% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function rgb = multispectral_mapping(img)

    [m, n, n_bands] = size(img);
    rgb = zeros(m, n, 3);

    switch n_bands
        case 1
            rgb = mat2gray(repmat(img, 1, 1, 3));
        case 2
            rgb(:, :, 1) = mat2gray(img(:, :, 2));
            rgb(:, :, 2) = mat2gray(img(:, :, 1));
        case 3
            rgb = mat2gray(img);
        otherwise
            coeff = pca(reshape(img, m * n, n_bands)', ...
                'NumComponents', 3);
            for i = 1:3
                rgb(:, :, i) = mat2gray(reshape(coeff(:, i), m, n));
            end
    end
end

