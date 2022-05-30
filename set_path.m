% Author   : Philipp Flotho, modification: E. Hennestad
% Copyright 2022 by Philipp Flotho, All rights reserved.

currentdir = fileparts( mfilename('fullpath') );

addpath(fullfile(currentdir, 'core'));
addpath(fullfile(currentdir, 'util'));
addpath(fullfile(currentdir, 'util', 'io'));

run(fullfile(currentdir, 'core', 'make.m'));