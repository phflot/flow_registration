% Copyright 2021 by Philipp Flotho, All rights reserved.
% Author   : Philipp Flotho

function [ fx, fy ] = my_imgradient( f )
%MY_IMGRADIENT Straight forward central differences implementation
% assuming a zero padded input array
    fx = zeros(size(f));
    fy = zeros(size(f));
    
    fy(2:end-1, :) = 0.5 * (f(3:end, :) - f(1:end - 2, :));
    fx(:, 2:end-1) = 0.5 * (f(:, 3:end) - f(:, 1:end - 2));
end

