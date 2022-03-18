% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function [ J11, J22, J33, J12, J13, J23] = get_motion_tensor_gc...
    (f1, f2, hx, hy)

    f1 = padarray(f1, [1 1], 'symmetric');
    f2 = padarray(f2, [1 1], 'symmetric');
%     f1 = set_boundary(f1);
%     f2 = set_boundary(f2);

    [fx1, ~ ] = gradient(f1, hx, hy);
    [fx2, ~ ] = gradient(f2, hx, hy);
    
    fx = 0.5 * (fx1 + fx2);
    ft = f2 - f1;
    
    fx = padarray(fx(2:end-1, 2:end-1), [1 1], 'symmetric');
    ft = padarray(ft(2:end-1, 2:end-1), [1 1], 'symmetric');
    
%     fx = set_boundary(fx);
%     ft = set_boundary(ft);
    [ ~ , fxy] = gradient(fx, hx, hy);
    [fxt, fyt] = gradient(ft, hx, hy);
    
%     f1 = set_boundary(f1);
%     f2 = set_boundary(f2);
    [fxx1, fyy1] = gradient2(f1, hx, hy);
    [fxx2, fyy2] = gradient2(f2, hx, hy);
    fxx = 0.5 * (fxx1 + fxx2);
    fyy = 0.5 * (fyy1 + fyy2);
    
    % dataterm normalization from Zimmer et al.:
    reg_x = 1 ./ (sqrt(fxx .* fxx + fxy .* fxy).^2 + 0.000001);
    reg_y = 1 ./ (sqrt(fxy .* fxy + fyy .* fyy).^2 + 0.000001);
    
    % gray value constancy with normalized dataterm:
    J11 = reg_x .* fxx .* fxx + reg_y .* fxy .* fxy;
    J22 = reg_x .* fxy .* fxy + reg_y .* fyy .* fyy;
    J33 = reg_x .* fxt .* fxt + reg_y .* fyt .* fyt;
    J12 = reg_x .* fxx .* fxy + reg_y .* fxy .* fyy;
    J13 = reg_x .* fxx .* fxt + reg_y .* fxy .* fyt;
    J23 = reg_x .* fxy .* fxt + reg_y .* fyy .* fyt;
    
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
