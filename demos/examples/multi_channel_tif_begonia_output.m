% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

clear;
run('../../set_path.m');

ch1 = 'ch1.TIF';
ch2 = 'ch2.TIF';

results_folder = 'results_begonia/';

options = OF_options(...
    'input_file', {ch1, ch2}, ... % input path
    'output_path', results_folder, ... % results folder
    'output_format', 'BEGONIA', ... % output file format: HDF5, MAT, 
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

% add this, to set begonia metaparameters (here 500hz and experiment name):
options.output_file_writer = get_video_file_writer(fullfile(...
    results_folder, 'compensated.h5'), 'BEGONIA', 'dt', 1/500, 'name', 'experiment name');

% saving the options to txt file (for archiving):
options.save_options('results/options.json');

compensate_recording(options);
