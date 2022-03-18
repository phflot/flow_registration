Please download the "2-Photon Movies with Motion Artifacts" from Dryad and put the 
hdf5 files into a folder data/

The scripts reproduce the supplemental video results from our work
P. Flotho, S. Nomura, B. Kuhn and D. J. Strauss, “Software for non-parametric Image Registration of 2-Photon Imaging Data”

Please feel free to try the "fast" versions (indicated with fast_ prefix) of the scripts with approximated solutions. 
On 6 Cores and SSD, they should take around 1min to complete and the results are comparable to the complete solutions. 
The complete solutions (qu_ prefix) will take around 3-4min. For a compromise between quality and speed,
set the quality setting to 'medium' or min_level to 4. 

For better memory efficiency during the compensation of the 30Hz files, the buffer size could be reduced. 