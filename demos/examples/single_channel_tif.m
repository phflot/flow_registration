% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

clear;
run('../../set_path.m');

options = OF_options(...
    'input_file', 'image.tif', ... % input path
    'output_path', 'results/', ... % results folder
    'output_format', 'TIFF', ... % output file format: HDF5, MAT or TIFF
    'alpha', 1.5, ... % smoothness parameter
    'sigma', [2, 2, 0.1; ...  % gauss kernel size channel 1
              2, 2, 0.1], ... % gauss kernel size channel 2
    'quality_setting', 'balanced', ... % set the quality out of 'fast', 'medium' or 'quality'
    'bin_size', 1, ... % binning over 5 frames from the 30 hz data
    'buffer_size', 22, ... % size of blocks for the parallel evaluation (larger takes more memory)
    'reference_frames', 100:200 ...
    );

% saving the options to txt file (for archiving):
options.save_options('results/options.json');

compensate_recording(options);