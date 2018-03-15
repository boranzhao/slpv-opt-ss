%% install 

filename = mfilename('fullpath');
[pathstr,name,ext] = fileparts(filename);
addpath(pathstr,[pathstr '\data'], [pathstr '\PSOmatlab'], [pathstr '\PSOmatlab\private'], [pathstr '\PSOmatlab\testfcns']);
