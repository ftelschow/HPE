%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%    Precompute different types of random fields Hermite Projector
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: This script precomputes the random fields used in the matlab
% script Simulation_LKCestimators.m. This procedure will speed up the
% simulations quite considerably.
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
addpath( genpath( path_RFT ) );

%------ define parameters for generation of random fields
% number of generated random fields
Ndata = 1e4;
% Field parameters
TYPE  = [ "isotropic" "scale-space" "nongauss" ];
D     = 2;           % Dimension of domain (including scale)
T     = 50;          % size of domain
nu    = 5;           % bandwidth parameter for isotropic and nongauss
gamma = 4:.2:15;     % bandwiths for scale-space

save( strcat(path_data, 'RandomFields_params') )

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%------------------ Generate the random fields ----------------------------
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this bit requires approximatively 70seconds on a standard laptop with
% 16gB ram and 4 cores.
%__________________________________________________________________________
% set seed to make sure the output is always the same
rng(42)
%------ isotropic
tic
type = "isotropic";
[eps, L ] = generateField( Ndata, 1, T*ones([1 D]), type, nu );
outputname = strcat( path_data, 'RandomFields_', 'Ndata', num2str(Ndata),'_',...
                     type, '_D', num2str(D), 'T', num2str(T),...
                    '_params', num2str(nu), '.mat' );
save(outputname, 'eps', 'L','-v7.3')

%-- scale-space
type = "scale-space";
outputname = strcat( path_data, 'RandomFields_', 'Ndata', num2str(Ndata),'_',...
                     type, '_D', num2str(D), 'T', num2str(T),...
                    '_params', num2str(gamma(1)),...
                    '_', num2str(gamma(end)), '.mat' );
[eps, L] = generateField( Ndata, 1, T*ones([1 D]), type, gamma );
save(outputname, 'eps', 'L','-v7.3')

%-- non-gauss
type = "nongauss";
outputname = strcat( path_data, 'RandomFields_', 'Ndata', num2str(Ndata),'_',...
                     type, '_D', num2str(D), 'T', num2str(T),...
                    '_params', num2str(nu), '.mat' );
[eps, L] = generateField( Ndata, 1, T*ones([1 D]), type, nu );
save(outputname, 'eps', 'L','-v7.3')
toc