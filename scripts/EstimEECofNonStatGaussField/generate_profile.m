%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%	This script produces a .mat containing the paths for use in the
%%%%	other scripts to make it easily reproducible.
%%%%    Moreover, it sets up all other required steps like compiling the
%%%%    mex files
%%%%
%%%%    Output: profile.mat
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Fabian Telschow
%--------------------------------------------------------------------------

% clear workspace
clear all
close all

% change accordingly to your folder structure
path_main = '/home/drtea/matlabToolboxes/HPE/';
path_pics = strcat( path_main, 'pics/EstimEECofNonStatGaussField/' );
path_data = strcat( path_main, 'data/EstimEECofNonStatGaussField/' );
path_results = strcat( path_main, 'results/EstimEECofNonStatGaussField/' );

% path to matlab toolboxes
path_spm12 = '/home/drtea/matlabToolboxes/spm12/';

% change to main directory
cd( path_main )
mkdir 'pics'
mkdir 'data'
mkdir 'results'
mkdir( path_pics )
mkdir( path_data )
mkdir( path_results )


save( strcat( path_main, '/scripts/EstimEECofNonStatGaussField/paths.mat' ),...
                         'path_main', 'path_pics', 'path_spm12', ...
                         'path_data', 'path_results' )
% precompute random fields
run( 'scripts/EstimEECofNonStatGaussField/Precompute_RandomFields.m' );
                                         
% compile the c++ code
cd( strcat( path_main,'/code/EEC/csource/') )
mex EulerCharCrit_c.cpp
cd( path_main )