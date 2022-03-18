% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function [avg_psnr, avg_mse, avg_std] = get_metrics(ch1, ch2, ref_idx)
%GET_METRICS Calculates the metrics as reported in the paper

    [~, ~, t] = size(ch1);

    sigma = [3 3 0.0001];
    b = 25;
    
    ch1_low = imgaussfilt3(double(ch1), sigma);
    ch2_low = imgaussfilt3(double(ch2), sigma);
    
    notIdx = setdiff(1:t, ref_idx);
    
    ch1_ref = mean(ch1_low(:, :, ref_idx), 3);
    ch2_ref = mean(ch2_low(:, :, ref_idx), 3);
    
    ch1_low = ch1_low(:, :, notIdx);
    ch2_low = ch2_low(:, :, notIdx);
    
    avg_psnr = zeros(1, length(notIdx));
    avg_mse = zeros(1, length(notIdx));
    
    for i = 1:length(notIdx)
        avg_psnr(i) = 0.5 * (...
            psnr(ch1_ref(b:end-b, b:end-b), ch1_low(b:end-b, b:end-b, i)) + ...
            psnr(ch2_ref(b:end-b, b:end-b), ch2_low(b:end-b, b:end-b, i)));
        
        avg_mse(i) = 0.5 * (...
            immse(ch1_ref(b:end-b, b:end-b), ch1_low(b:end-b, b:end-b, i)) + ...
            immse(ch2_ref(b:end-b, b:end-b), ch2_low(b:end-b, b:end-b, i)));
    end
    
    avg_std = 0.5 * (...
        mean(mean(std(ch1_low, 0, 3), 1), 2) + ...
        mean(mean(std(ch2_low, 0, 3), 1), 2));
end

