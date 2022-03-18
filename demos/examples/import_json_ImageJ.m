% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

% Example to demonstrate how to load a json settings file generated with ImageJ

clear;
run('../../set_path.m');

options = OF_options;
options.load_options('options.json');

options.output_format = 'HDF5';
options.input_file = 'inputfile.tiff';
options.output_path = 'results';

compensate_recording(options);