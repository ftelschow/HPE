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
% clear workspace
clear all
close all

%------ change accordingly to your folder structure
path_main       = '/home/drtea/matlabToolboxes/HPE/';
path_RFTtoolbox = '/home/drtea/matlabToolboxes/RFTtoolbox/';
path_spm12      = '/home/drtea/matlabToolboxes/spm12/';

%------ related paths and save them
path_pics    = fullfile( path_main, 'pics', 'LKCestimation3D' );
path_data    = fullfile( path_main, 'data', 'LKCestimation3D' );
path_results = fullfile( path_main, 'results', 'LKCestimation3D' );

% save the path identifiers
save( strcat( path_main, '/scripts/LKCestimation3D/paths.mat' ),...
                         'path_main', 'path_pics', ...
                         'path_data', 'path_results',...
                         'path_RFTtoolbox', 'path_spm12' )

%------ change to main directory and produce necessary folder substructure
cd( path_main )
mkdir 'pics'
mkdir 'data'
mkdir 'results'
cd pics
mkdir 'LKCestimation3D'
cd ..
cd data
mkdir 'LKCestimation3D'
cd ..
cd results
mkdir 'LKCestimation3D'
cd ..
                     
%------ precompute random fields for speed up
cd( strcat( path_main, 'scripts/LKCestimation3D/' ) )
run( 'Precompute_RandomFields3D.m' );
