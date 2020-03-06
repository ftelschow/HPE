%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%    Precompute different types of random fields over 3D domain for 
%%%%    studying the performance of LKC estimators
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: This script precomputes the random fields used in the matlab
% script Simulation_LKCestimators3D.m
%__________________________________________________________________________
% REFERENCES:
%
%__________________________________________________________________________
% AUTHOR: Fabian Telschow (ftelschow@ucsd.edu)
%__________________________________________________________________________
%------ prepare workspace
% clean history
clear
close all
% load generated profile containing the paths
load('paths.mat')

% create matlab paths to search for functions
addpath(genpath(path_main));
addpath(genpath(path_RFTtoolbox));

%------ parameters for computation of random fields
% number of random fields divided into N computations of Nbatch for memory
% purposes. A total of N*Nbatch random fields is saved
N = 400;
Nbatch = 10;

% random field parameters
D    = 3;
T    = 50;
Dim  = repmat(T, [1 D]);
FWHM = [ 3 6 12 15 ];

%------ general constants
band2FWHM = sqrt(8*log(2));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%------------------------- Generate random fields -------------------------
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%--------------------------------------------------------------------------
%------ Isotropic Gaussian fields
tic
for f = FWHM
    % Get true LKC
    alpha = 1 / ( 4 * (f/band2FWHM)^2 );
    switch D
        case 1
            LKC = ( T - 1 ) * sqrt( 2 * alpha );
        case 2
            LKC = [ 2 * ( T - 1 ); ( T - 1 )^2 ]...
                  .* [ sqrt( 2 * alpha ); 2 * alpha ] ;
        case 3
            LKC = [ 3 * ( T - 1 ); 3 * ( T - 1 )^2; ( T - 1 )^3 ]...
                  .* [ sqrt( 2 * alpha ); 2 * alpha; ( 2 * alpha )^( 3/2 ) ] ;
    end
    % generate random fields
    rfs = [];
    for k = 1:Nbatch
            rfs_tmp = genRF( Dim, 1, f, N, 1 ); % option 0 not working, ask Sam!
            rfs = cat(D+1, rfs, permute(reshape(rfs_tmp, [N, Dim]), [2:(D+1) 1]));
        clear rfs_tmp 
    end

    save( fullfile( path_data, strcat("Isotropic_FWHM", num2str(f), "_T",...
                    num2str(T), "_D", num2str(D), "_Nfields",...
                    num2str(Nbatch*N) )), 'rfs', 'LKC', '-v7.3' )

end
toc

% 
%--------------------------------------------------------------------------
%------ Scale space Gaussian fields
rng(23)
params = 4:0.2:15;
tic
    % generate random fields
    rfs = [];
    for k = 1:Nbatch
            rfs_tmp = generateField( N, 1, [T T length(params)], ...
                              "scale-space", params );
            rfs = cat(D+1, rfs, rfs_tmp);
        clear rfs_tmp 
    end
toc
save( fullfile( path_data, strcat("ScaleSpace__T",...
                num2str(T), "_D", num2str(D), "_Nfields",...
                num2str(N*Nbatch) )), 'rfs', 'LKC', 'params', '-v7.3' )