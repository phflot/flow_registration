Flow-Registration

Compile and install the toolbox by running the script set_path.m

The folder demos/ contains the standalone jupiter demos, examples on how to use the code
and scripts to reproduce the journal results. 

The Code has been tested with MATLAB R2018a and R2021a.

To get started, define an OF_options object 

options = OF_options(...
    'input_file', 'inputfile.tiff/mdf/hdf5/MAT', ...
    'output_path', 'results, ... 
    'output_format', 'HDF5', ...
    'quality_setting', 'balanced', ... % default behaviour is 'quality'
    'reference_frames', 1:100 ...
    );

and run

compensate_recording(options);

to compensate the input file and write the results into the results folder. 