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
%--------------------------------------------------------------------------
%   Author: Fabian Telschow
%--------------------------------------------------------------------------

% clear workspace
clear all
close all

%------ change accordingly to your local folder structure
% path to git repository HPE
path_main = '/home/fabian/Seafile/Code/matlabToolboxes/HPE/';
% path to spm12 toolbox
path_spm12 = '/home/fabian/Seafile/Code/matlabToolboxes/spm12/';

%------ derived paths and save them
path_pics = strcat( path_main, 'pics/EstimEECofNonStatGaussField/' );
path_data = strcat( path_main, 'data/EstimEECofNonStatGaussField/' );
path_results = strcat( path_main, 'results/EstimEECofNonStatGaussField/' );
path_RFT = strcat( path_main, 'code/');

% save the path identifiers
save( strcat( path_main, '/scripts/EstimEECofNonStatGaussField/paths.mat' ),...
                         'path_main', 'path_pics', 'path_spm12', ...
                         'path_data', 'path_results', 'path_RFT' )

%------ produce necessary folder substructure
cd( path_main )
mkdir 'pics'
mkdir 'data'
mkdir 'results'
mkdir( path_pics )
mkdir( path_data )
mkdir( path_results )

%------ precompute random fields for speed up
run( 'scripts/EstimEECofNonStatGaussField/Precompute_RandomFields.m' );
                                         
%------ compile the c++ code
cd( strcat( path_main,'/code/EEC/csource/') )
mex EulerCharCrit_c.cpp
cd( path_main )