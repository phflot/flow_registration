% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.
% Example that downloads the jupiter demo data and demonstrates the minimal motion compensation config:

clear;
run('../set_path.m');

% preparing the data:
output_folder = 'jupiter_demo';
if ~isfolder(output_folder)
    mkdir(output_folder)
end
input_file = fullfile(output_folder, 'jupiter.tiff');
if (~exist(input_file, 'file'))
    websave(input_file, ...
        'https://cloud.hiz-saarland.de/s/JpHyczRSMDbLwzP/download');
end

options = OF_options(...
    'input_file', input_file, ...
    'output_path', fullfile(output_folder, 'hdf5_comp_minimal'), ... 
    'output_format', 'HDF5', ...
    'alpha', 4, ... 
    'quality_setting', 'balanced', ...
    'output_typename', '', ...
    'n_references', 2, ...
    'min_frames_per_reference', 20, ...
    'reference_frames', 1:250 ...
    );
compensate_recording(options);

vid = get_video_file_reader('jupiter_demo/hdf5_comp_minimal/compensated.HDF5');
vid = vid.read_frames(1:vid.frame_count);
energy = zeros(1, size(vid, 4));
ref = squeeze(mean(vid(:, :, :, 100:200), 4));
parfor i = 1:size(vid, 4)
    energy(i) = get_energy(vid(:, :, :, i), ref);
end

implay(mat2gray(vid));
load('jupiter_demo/hdf5_comp_minimal/statistics.mat')
load('jupiter_demo/hdf5_comp_minimal/reference_frames.mat')
ref1 = energy;
ref2 = energy;
ref2(idx == 1) = nan;
ref1(idx == 2) = nan;
figure;
subplot(2, 1, 1);
hold on
plot(ref1);
plot(ref2);
title("Energy");
subplot(2, 1, 2);
imshowpair(reference_frames{1}, reference_frames{2}, 'montage');
title("Reference frames:")

