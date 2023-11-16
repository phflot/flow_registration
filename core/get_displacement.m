% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function w = get_displacement( fixed, moving, varargin )
% computes the displacements

    alpha = [2, 2];
    update_lag = 10;
    iterations = 20;
    
    min_level = 0;
        
    levels = 50;
    eta = 0.75;
    
    a_smooth = 0.5;
    
    [m, n, n_channels] = size(fixed);
    
    u_init = zeros(m, n);
    v_init = zeros(m, n);
    
    weight = ones(1, n_channels, 'double') / n_channels;

    const_assumption = "gc";
    
    a_data = 0.45 * ones(1, n_channels);
    
    for k = 1:length(varargin)
        if ~isa(varargin{k}, 'char')
            continue;
        end
        switch varargin{k}
            case 'weight'
                weight = varargin{k + 1};
            case 'alpha'
                alpha = varargin{k + 1};
                if length(alpha) == 1
                    alpha = alpha .* ones(1, 2);
                end
            case 'eta'
                eta = varargin{k + 1};
            case 'levels'
                levels = varargin{k + 1};       
            case 'update_lag'
                update_lag = varargin{k + 1};
            case 'iterations'
                iterations = varargin{k + 1};
            case 'uv'
                u_init = varargin{k + 1};
                v_init = varargin{k + 2};
            case 'a_data'
                a_data = varargin{k + 1};
                if (length(a_data) == 1)
                    a_data = a_data * ones(1, n_channels);
                end
            case 'a_smooth'
                a_smooth = varargin{k + 1};
            case 'min_level'
                min_level = varargin{k + 1};
            case 'constancy_assumption'
                const_assumption = varargin{k + 1};
            otherwise
                % fprintf(['could not parse input argument ' varargin{k} '\n']);
        end
    end
        
    switch const_assumption
        case 'gc'
            get_motion_tensor = @get_motion_tensor_gc;
        case 'gray'
            get_motion_tensor = @get_motion_tensor;
        otherwise
            fprintf("Could not parse constancy assumption, setting to default gc!\n");
            get_motion_tensor = @get_motion_tensor_gc;
    end


    f1_low = double(fixed);
    f2_low = double(moving);
    
    method = 'bicubic';
    
    max_level_y = warpingDepth(eta, levels, m, m);
    max_level_x = warpingDepth(eta, levels, n, n);
    
    max_level = min(max_level_x, max_level_y) * 4;
    
    max_level_y = min(max_level_y, max_level);
    max_level_x = min(max_level_x, max_level);
    
    local_weight = ndims(weight) == ndims(fixed) && sum(size(weight) == size(fixed)) == ndims(fixed);
    weight_level = weight;
    
    if max(max_level_x, max_level_y) <= min_level
        min_level = max(max_level_x, max_level_y) - 1;
    end
    if min_level < 0
        min_level = 0;
    end
        
    for i = max(max_level_x, max_level_y):-1:min_level       
        level_size = round([m * eta^(min(i, max_level_y)), ...
            n * eta^(min(i, max_level_x))]);
        
        f1_level = imresize(f1_low, ...
            level_size, ...
            method, 'Colormap', 'original', 'Antialiasing', true);
        f2_level = imresize(f2_low, ...
            level_size, ...
            method, 'Colormap', 'original', 'Antialiasing', true);
        
        if local_weight
            weight_level = padarray(imresize(weight, ...
                level_size, ...
                method, 'Colormap', 'original', 'Antialiasing', true), ...
                [1 1], 0.0);
        end
        
        hx = m / size(f1_level, 1);
        hy = n / size(f1_level, 2);
        
        if i == max(max_level_x, max_level_y)
            u = add_boundary(imresize(u_init, level_size, method, 'Colormap', 'original')); 
            v = add_boundary(imresize(v_init, level_size, method, 'Colormap', 'original')); 
            tmp = double(f2_level);
        else
            u = add_boundary(imresize(u(2:end - 1, 2:end - 1), level_size, method, 'Colormap', 'original'));
            v = add_boundary(imresize(v(2:end - 1, 2:end - 1), level_size, method, 'Colormap', 'original'));
            try
                tmp = imregister_wrapper(double(f2_level), ...
                    double(u(2:end-1, 2:end-1)) / hx, ...
                    double(v(2:end-1, 2:end-1)) / hy, double(f1_level));
            catch err
                disp(err.message);
                error("Error using imregister for compensating flow increments, try increasing alpha!");
            end
        end
        
        J11 = zeros([level_size + 2, n_channels]);
        J22 = zeros([level_size + 2, n_channels]);
        J33 = zeros([level_size + 2, n_channels]);
        J12 = zeros([level_size + 2, n_channels]);
        J13 = zeros([level_size + 2, n_channels]);
        J23 = zeros([level_size + 2, n_channels]);
        for j = 1:n_channels
            [ J11(:, :, j), J22(:, :, j), J33(:, :, j), ...
                J12(:, :, j), J13(:, :, j), J23(:, :, j)] = ...
                get_motion_tensor(...
                f1_level(:, :, j), tmp(:, :, j), hx, hy);
        end
        
        if i == min_level
            alpha_scaling = 1;
        else
            alpha_scaling = eta.^(-0.5 * i);
        end
        
        [du, dv] = level_solver(J11, J22, J33, J12, J13, J23, weight_level, ... 
            u, v, alpha * alpha_scaling, iterations, update_lag, 0, a_data, a_smooth, hx, hy);
        
        if min(level_size > 5)
            du(2:end-1, 2:end-1) = medfilt2(du(2:end-1, 2:end-1), [5 5], 'symmetric');
            dv(2:end-1, 2:end-1) = medfilt2(dv(2:end-1, 2:end-1), [5 5], 'symmetric');
        end
        
        u = u + du;
        v = v + dv;
    end
    
    w = zeros([size(u) - 2, 2], 'double');
    w(:, :, 1) = u(2:end-1, 2:end-1);
    w(:, :, 2) = v(2:end-1, 2:end-1);
    
    if min_level > 0
        w = imresize(w, [m, n]);
    end
end

function [du, dv] = OF_solver_GPU(J11, J22, J33, J12, J13, J23, weight_level, ... 
            u, v, alpha, iterations, update_lag, verbose, a_data, a_smooth, hx, hy)
        
    u_gpu = gpuArray(u);
    v_gpu = gpuArray(v);
        
    du = gpuArray(zeros(size(u)));
    dv = gpuArray(zeros(size(u)));
    psi_data = gpuArray(ones(size(J11) - [2, 2, 0]));
    
    [m, n] = size(du);
    c1 = gpuArray(2:m-1);
    c2 = gpuArray(2:n-1);
    
    tmp1 = alpha(1) / hx.^2;
    tmp2 = alpha(2) / hy.^2;
    tmp3 = 2 * tmp1 + 2 * tmp2;

    kernel = [0, tmp1, 0;
              tmp2, 0, tmp2;
              0, tmp1, 0];
          
    j11 = gpuArray(J11);
    j22 = gpuArray(J22);
    j33 = gpuArray(J33);
    j12 = gpuArray(J12);
    j23 = gpuArray(J23);
    j13 = gpuArray(J13);
    
    for iteration_counter = 0:iterations-1
        
        if mod(iteration_counter, update_lag) == 0
            % update non-linearities
            psi_data = j11 .* du.^2 + j22 .* dv.^2 + j23 .* dv + ...
                2 * j12 .* du .* dv + 2 * j13 .* du + j23 .* dv + j33;
            psi_data(psi_data < 0) = 0;
            psi_data = a_data(1) * (psi_data + 0.00001).^(a_data(1) - 1);
        end
        
        set_boundary(du);
        set_boundary(dv);
        

%         l = 1:n-2;
%         r = 3:n;
%         t = 1:m-2;
%         b = 3:m;        
%         
%         smooth_u = 0;
%         smooth_u = smooth_u + tmp1 * (u(c1, l) + du(c1, l) - u(c1, c2));
%         smooth_u = smooth_u + tmp1 * (u(c1, r) + du(c1, r) - u(c1, c2));
%         smooth_u = smooth_u + tmp2 * (u(t, c2) + du(t, c2) - u(c1, c2));
%         smooth_u = smooth_u + tmp2 * (u(b, c2) + du(b, c2) - u(c1, c2));

        denom_u = 0; %tmp3;
        denom_v = denom_u;

%         smooth_u = conv2(u + du, kernel, 'valid') - tmp3 * u(c1, c2);
%         smooth_v = conv2(v + dv, kernel, 'valid') - tmp3 * v(c1, c2);
        
%         smooth_v = 0;
%         smooth_v = smooth_v + tmp1 * (v(c1, l) + dv(c1, l) - v(c1, c2));
%         smooth_v = smooth_v + tmp1 * (v(c1, r) + dv(c1, r) - v(c1, c2));
%         smooth_v = smooth_v + tmp2 * (v(t, c2) + dv(t, c2) - v(c1, c2));
%         smooth_v = smooth_v + tmp2 * (v(b, c2) + dv(b, c2) - v(c1, c2));

        num_u = sum(- psi_data .* (J13 + J12 .* dv), 3);
        denom_u = denom_u + sum(psi_data .* J11, 3);
        denom_v = denom_v + sum(psi_data .* J22, 3);
        
        du(c1, c2) = (1 - 1.95) * du(c1, c2) + ...
            1.95 * ((num_u(c1, c2, :)) ./ denom_u(c1, c2));
        
        num_v = sum(- psi_data .* (J23 + J12 .* du), 3);
        
        dv(c1, c2) = (1 - 1.95) * dv(c1, c2) + ...
            1.95 * ((num_v(c1, c2, :)) ./ denom_v(c1, c2));
    end 
    dv = (1-tmp1) * (dv + v_gpu) + tmp1 * imgaussfilt(v_gpu + dv, 10) - v_gpu;
    du = (1-tmp1) * (du + u_gpu) + tmp1 * imgaussfilt(u_gpu + du, 10) - u_gpu;
    
    du = gather(du);
    dv = gather(dv);
end

function f = add_boundary(f)
    f = padarray(f, [1 1]);
    
    f = set_boundary(f);
end

function f = set_boundary(f)
    f(:, 1) = f(:, 2);
    f(:, end) = f(:, end - 1);
    f(1, :) = f(2, :);
    f(end, :) = f(end - 1, :);
end
