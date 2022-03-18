% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

clear;
run('../../set_path.m');

ch1 = 'ch1.TIF';
ch2 = 'ch2.TIF';

options = OF_options(...
    'input_file', {ch1, ch2}, ... % input path
    'output_path', 'results_caiman/', ... % results folder
    'output_format', 'CAIMAN_HDF5', ... % output file format: HDF5, MAT, 
                                    ... % TIFF or MULTIFILE_HDF5, ... to generate multiple files
                                    ... % or CAIMAN_HDF5 for CAIMAN support
    'alpha', 1.5, ... % smoothness parameter
    'sigma', [1, 1, 0.1; ...  % gauss kernel size channel 1
              1, 1, 0.1], ... % gauss kernel size channel 2
    'quality_setting', 'balanced', ... % set the quality out of 'fast', 'medium' or 'quality'
    'bin_size', 1, ... % binning over 5 frames from the 30 hz data
    'buffer_size', 24, ... % size of blocks for the parallel evaluation (larger takes more memory)
    'reference_frames', 100:200 ...
    );

% saving the options to txt file (for archiving):
options.save_options('results/options.json');

compensate_recording(options);
