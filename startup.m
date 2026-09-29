%% FlutterWing startup: add the repo root and all its subfolders to the MATLAB path
% Run from the repo root (run startup.m) or from anywhere via run('/path/to/flutterwing/startup.m').
projRoot = fileparts(mfilename("fullpath"));
addpath(genpath(projRoot));
disp(['FlutterWing: added ', projRoot, ' and its subfolders to the MATLAB path'])
