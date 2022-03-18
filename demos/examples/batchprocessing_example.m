% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

clear;
run('../../set_path.m');

file1 = 'file1.mdf';
file2 = 'file2.mdf';

options = OF_options(...
    'output_path', 'batch_processing_results/', ... % results folder
    'output_format', 'TIFF', ... % output file format: HDF5, MAT, 
                             ... % TIFF or MULTIFILE_HDF5, ... to generate multiple files
                             ... % or CAIMAN_HDF5 for CAIMAN support
    'alpha', 1.5, ... % smoothness parameter
    'sigma', [1, 1, 0.1; ...  % gauss kernel size channel 1
              1, 1, 0.1], ... % gauss kernel size channel 2
    'quality_setting', 'balanced', ... % set the quality out of 'fast', 'medium' or 'quality'
    'bin_size', 5, ... % binning over 5 frames from the 30 hz data
    'buffer_size', 24, ... % size of blocks for the parallel evaluation (larger takes more memory)
    'reference_frames', 100:200 ...
    );

% saving the options to txt file (for archiving):
options.save_options('results/options.json');

batchprocessor = OF_batchprocessor(...
    options, 'reference_mode', 'same'); % reference_mode can be 'same' for 
                                        % the same reference for all jobs 
                                        % or single, where the reference is
                                        % fetched for each job individually

batchprocessor.append_job(file1);
batchprocessor.append_job(file2);

% % syntax for all mdf files in folder for example:
% files = dir(fullfile('test', '*.mdf'));
% for i = 1:length(files)
%     batchprocessor.append_job(fullfile(...
%         files(i).folder, files(i).name));
% end

batchprocessor.start_all;