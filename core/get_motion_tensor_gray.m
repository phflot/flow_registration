% Author   : Philipp Flotho
% Copyright 2023 by Philipp Flotho, All rights reserved.

function [ J11, J22, J33, J12, J13, J23] = get_motion_tensor_gray...
    (f1, f2, hx, hy)

    f1 = padarray(f1, [1 1], 'symmetric');
    f2 = padarray(f2, [1 1], 'symmetric');

    [fx1, ~ ] = gradient(f1, hx, hy);
    [fx2, ~ ] = gradient(f2, hx, hy);
    
    fx = 0.5 * (fx1 + fx2);
    ft = f2 - f1;
    
    fx = padarray(fx(2:end-1, 2:end-1), [1 1], 'symmetric');
    ft = padarray(ft(2:end-1, 2:end-1), [1 1], 'symmetric');

    J11 = fx .* fx;
    J22 = fy .* fy;
    J33 = ft .* ft;
    J12 = fx .* fy;
    J13 = fx .* ft;
    J23 = fy .* ft;
    
    J11 = set_boundary0(J11);
    J22 = set_boundary0(J22);
    J33 = set_boundary0(J33);
    J12 = set_boundary0(J12);
    J13 = set_boundary0(J13);
    J23 = set_boundary0(J23);
end

function f = set_boundary0(f)
    f(:, 1) = 0;
    f(:, end) = 0;
    f(1, :) = 0;
    f(end, :) = 0;
end

function f = set_boundary(f)
    f(:, 1) = f(:, 3);
    f(:, end) = f(:, end - 2);
    f(1, :) = f(3, :);
    f(end, :) = f(end - 2, :);
end

function [fxx, fyy] = gradient2(f, hx, hy)
    fxx = zeros(size(f));
    fyy = zeros(size(f));

    fxx(2:end-1, 2:end-1) = ...
        (f(2:end-1, 1:end-2) - 2 * f(2:end-1, 2:end-1) ...
        + f(2:end-1, 3:end)) ./ hx^2;
    
    fyy(2:end-1, 2:end-1) = ...
        (f(1:end-2, 2:end-1) - 2 * f(2:end-1, 2:end-1) ...
        + f(3:end, 2:end-1)) ./ hy^2;
    
%     fyy = padarray(fyy(3:end-2, 3:end-2), [2, 2], 'replicate');
%     fxx = padarray(fxx(3:end-2, 3:end-2), [2, 2], 'replicate');
end
