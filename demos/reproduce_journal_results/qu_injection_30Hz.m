% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

% script for the complete compensation of the injection 30Hz sequence. Please
% load the dataset into data/
run('../../set_path.m');

input_file = 'data/injection_30hz.HDF';
output_folder = 'injection_30Hz_fast_result';

if ~isfolder(output_folder)
    mkdir(output_folder);
end

if ~isfile(input_file)
    error('Please download the benchmark data into the data/ folder first!');
end

run('../../set_path.m');

vid_reader = get_video_file_reader(input_file, 500, 1);

ref_idx = 1:500;

options = OF_options(...
    'input_file', vid_reader, ... % input path
    'output_path', output_folder, ... % results folder
    'output_format', 'HDF5', ...
    'quality_setting', 'quality', ...
    'sigma', [1, 1, 0.5; ... % gauss kernel size channel 1
              1, 1, 0.5], ... % gauss kernel size channel 2
    'weight', [1.15, 0.85], ...
    'buffer_size', 2500, ... % size of blocks for the parallel evaluation (larger takes more memory)
    'reference_frames', ref_idx ...
    );

% % NoRMCorre Parameters:
% options_nr = NoRMCorreSetParms(...
%     'd1', 512, 'd2', 512, ...
%     'grid_size', [64, 64], ...
%     'max_dev', [25, 25, 1], 'max_shift', 25, 'use_parallel', true, ...
%     'shifts_method', 'FFT', 'upd_template', false);

compensate_recording(options);

fprintf('starting visualization...\n');

vid_comp = get_video_file_reader(fullfile(output_folder, ...
    'compensated.HDF5'), 2500);

vid_reader.reset;
raw = vid_reader.read_batch;
compensated = vid_comp.read_batch;
avg = mean(raw, 4);

clear output
output = get_visualization(double(raw(:, :, 1, :)), double(raw(:, :, 2, :)), ...
    [1, 0], [1, 0], avg(:, :, 1), avg(:, :, 2));

output(:, end+1:end+size(avg, 2), :, :) = get_visualization(double(compensated(:, :, 1, :)), double(compensated(:, :, 2, :)), ...
    [1, 0], [1, 0], avg(:, :, 1), avg(:, :, 2));

implay(output, 60);

disp('computing the metrics...');
[flow_psnr, flow_mse, flow_std] = get_metrics(squeeze(compensated(:, :, 1, :)), ...
    squeeze(compensated(:, :, 2, :)), ref_idx);
[raw_psnr, raw_mse, raw_std] = get_metrics(squeeze(raw(:, :, 1, :)), ...
    squeeze(raw(:, :, 2, :)), ref_idx);

fprintf('Raw PSNR = %f, flow PSNR = %f\n', mean(raw_psnr), mean(flow_psnr));
fprintf('MSE performance factor = %f (figure 2 (A))\n', mean(raw_mse) / mean(flow_mse));
fprintf('STD performance factor = %f (figure 2 (B))\n', raw_std / flow_std);
