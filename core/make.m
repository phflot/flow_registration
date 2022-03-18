% Author   : Philipp Flotho
% Copyright 2021 by Philipp Flotho, All rights reserved.

function make(varargin)  
    clear functions;
    clear mex;

    if exist(fullfile(pwd, ['level_solver.' mexext]), 'file')
        delete(fullfile(pwd, ['level_solver.' mexext]));
    end
    if exist(fullfile(pwd, ['level_solver.' mexext '.pdb']), 'file')
        delete(fullfile(pwd, ['level_solver.' mexext '.pdb']));
    end
    
    mex(varargin{:}, 'level_solver.cpp');
end