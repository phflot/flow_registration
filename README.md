# Flow-Registration: Optical flow based motion compensation / video stabilization / registration for 2-photon imaging data

Toolbox for the compensation and stabilization of multichannel microscopy videos. The code is written in Matlab, Java (IJ Plugin) and C++. The preprint of our paper can be found on [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.07.25.453381v1) and the project website with video results [here](https://www.snnu.uni-saarland.de/flow-registration/).

![Fig1](img/bg.jpg)


## Installation

To install the toolbox clone the repo or download the most recent release and run the ```set_path.m``` script. With ```savepath``` the toolbox will be permanently available in MATLAB.

Please [contact us](mailto:Philipp.Flotho@uni-saarland.de) for more details.

## Requirements

To run this toolbox, MATLAB 2018b onwards with configured C++ compiler (check [here](https://www.mathworks.com/support/requirements/supported-compilers.html) for supported compilers) is required.

## Getting started

This repository contains the demo scripts ```demos/jupiter.m``` and ```demos/jupiter_minimal_example.m``` which run out of the box and compensate the jitter in an amateur recording of a meteor impact on jupiter. The folder ```demos/examples``` contains examples that illustrate use cases of the toolbox and ```demos/reproduce_journal_results``` contains scripts that replicate the evaluations from our paper.

The plugin supports most of the commonly used file types such as HDF5, tiff stacks and matlab mat files. To run the motion compensation, the options need to be defined into a ```OF_options``` object such as

```
options = OF_options(...
    'input_file', 'input.hdf', ... % input path
    'output_path', results_folder, ... % results folder
    'output_format', 'MAT', ...
    'alpha', [0.5, 0.5], ... % smoothness parameter
    'sigma', [0.5, 0.5, 0.1; ... % gauss kernel size channel 1
              0.5, 0.5, 0.1], ... % gauss kernel size channel 2
    'weight', [1, 1], ...
    'levels', 15, ... % solver levels
    'eta', 0.86, ... % pyramid stepsize
    'iterations', 25, ... % outer iterations (the larger the better the result, but slower)
    'bin_size', 5, ... % binning over 5 frames from the 30 hz data
    'buffer_size', 500, ... % size of blocks for the parallel evaluation (larger takes more memory)
    'output_typename', [], ...
    'reference_frames', 1:5 ...
    );
```

The object is then passed to ```compensate_recording(options)``` to run the motion compensation on ```input.hdf``` into ```results_folder```. To run ```compensate_recording``` with default parameters, only the input file and output path need to be specified.

## Dataset

The dataset which we used for our evaluations is available as [2-Photon Movies with Motion Artifacts](https://www.datadryad.org).

## Citation

Details on the method and video results can be found [here](https://www.snnu.uni-saarland.de/flow-registration/).

If you use parts of this code or the plugin for your work, please cite
  
> P. Flotho, S. Nomura, B. Kuhn and D. J. Strauss, “Software for Non-Parametric Image Registration of 2-Photon Imaging Data,” J Biophotonics, 2022. [doi:https://doi.org/10.1002/jbio.202100330](https://doi.org/10.1002/jbio.202100330)

BibTeX entry
```
@article{flotea2022a,
    author = {Flotho, P. and Nomura, S. and Kuhn, B. and Strauss, D. J.},
    title = {Software for Non-Parametric Image Registration of 2-Photon Imaging Data (in press)},
    year = {2022},
  journal = {J Biophotonics},
  doi = {https://doi.org/10.1002/jbio.202100330}
}
```