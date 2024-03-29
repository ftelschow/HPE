%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%	This script produces a .mat containing the paths for use in the
%%%%	other scripts for the 3D LKC estimation project to make it easily
%%%%    reproducible on any machine.
%%%%
%%%%    Output: profile.mat
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%   Author: Fabian Telschow
%--------------------------------------------------------------------------
%------ prepare workspace
% clear workspace
clear all
close all

%------ change accordingly to your local folder structure
path_main       = '/home/drtea/matlabToolboxes/HPE/';
path_RFTtoolbox = '/home/drtea/matlabToolboxes/RFTtoolbox/';
path_spm12      = '/home/drtea/matlabToolboxes/spm12/';

%------ derived paths and save them
path_pics    = fullfile( path_main, 'pics', 'LKCestimation3D' );
path_data    = fullfile( path_main, 'data', 'LKCestimation3D' );
path_tmp     = fullfile( path_data, 'tmp' );
path_results = fullfile( path_main, 'results', 'LKCestimation3D' );

% save the path identifiers
save( fullfile( path_main, 'scripts',   'LKCestimation3D', 'paths.mat' ),...
                           'path_main', 'path_pics', 'path_tmp',...
                           'path_data', 'path_results',...
                           'path_RFTtoolbox', 'path_spm12' )

%------ produce necessary folder substructure
cd( path_main )
mkdir 'pics'
mkdir 'data'
mkdir 'results'
cd pics
mkdir 'LKCestimation3D'
cd ..
cd data
mkdir 'LKCestimation3D'
cd LKCestimation3D
mkdir tmp
cd ..
cd ..
cd results
mkdir 'LKCestimation3D'
cd ..
                     
%------ precompute random fields for speed up
cd( fullfile( path_main, 'scripts', 'LKCestimation3D' ) )
% runtime for 4e3, 50^3-dim fields for FWHM [3 6 12 15] approx. 45min
run( 'Precompute_RandomFields3D.m' );
