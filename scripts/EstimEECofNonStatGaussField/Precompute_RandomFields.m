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
TYPE  = ["isotropic" "scale-space" "nongauss" "ng_scale"];
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

%-- non-gauss
type = "ng_scale";
outputname = strcat( path_data, 'RandomFields_', 'Ndata', num2str(Ndata),'_',...
                     type, '_D', num2str(D), 'T', num2str(T),...
                    '_params', num2str(gamma(1)),...
                    '_', num2str(gamma(end)), '.mat' );
[eps, L] = generateField( Ndata, 1, T*ones([1 D]), type, gamma );
save(outputname, 'eps', 'L','-v7.3')
toc

%%
%-- non-stat-Gauß expsquare
%
addpath(genpath("/home/fabian/Seafile/Code/matlabToolboxes/RFTtoolbox"))

D      = 2;
T      = 21;
Dim    = [T,T];
pars   = [[2.5, 0.05]; [2.5, -0.05]];

type = "nonstatgauss_exp";
outputname = strcat( path_data, 'RandomFields_', 'Ndata', num2str(Ndata),'_',...
                     type, '_D', num2str(D), 'T', num2str(T), '.mat' );

resadd = 1;
dx = 1/(resadd + 1);
xx = -10:10;

pad     = 20;
mask    = true( Dim );
mask    = logical( pad_vals( mask, pad ) );
mask_hr = mask_highres(mask, resadd);
x       = (-(10+pad)-dx):dx:((10+pad)+dx);
[X,Y]   = meshgrid(x, x);


data = wfield( mask_hr, Ndata);

%[min(pars(1,1) * exp(pars(1,2) * xx)), max(pars(1,1) * exp(pars(1,2) * xx))] 
%[min(pars(2,1) * exp(pars(2,2) * xx)), max(pars(2,1) * exp(pars(2,2) * xx))] 

val1 = exp( -(x - Y).^2 ./ (sqrt(2) * pars(1,1) * exp(pars(1,2) * x)).^2) * dx;
val2 = exp( -(x - Y).^2 ./ (sqrt(2) * pars(2,1) * exp(pars(2,2) * x)).^2) * dx;

sd_val = sqrt(pi) ./ sqrt( (pars(2,1) * pars(1,1)) * (exp(pars(1,2) * X)  .* exp(pars(2,2) * Y) )) * pi;

tic
for k = 1:data.fibersize
    data.field(:,:,k) = val1' * data.field(:,:,k) * val2 ./ sd_val;
end
ende = length(data.xvals{1});

data = data((pad + pad*resadd+1):(ende-((pad + pad*resadd))), (pad + pad*resadd+1):(ende-((pad + pad*resadd))),:);
data = data./std(data);
eps  = data.field;
L = [12.53, 39.25];

save(outputname, 'eps', 'L','-v7.3')
toc


%%
%-- non-stat-non-Gauß expsquare
%
addpath(genpath("/home/fabian/Seafile/Code/matlabToolboxes/RFTtoolbox"))

D      = 2;
T      = 21;
Dim    = [T,T];
pars   = [[2.5, 0.05]; [2.5, -0.05]];

type = "nonstatnongauss_exp";
outputname = strcat( path_data, 'RandomFields_', 'Ndata', num2str(Ndata),'_',...
                     type, '_D', num2str(D), 'T', num2str(T), '.mat' );


Nsubj  = 100;
resadd = 1;
dx = 1/(resadd + 1);
xx = -10:10;

pad = 20;
mask    = true( Dim );
mask    = logical( pad_vals( mask, pad ) );
mask_hr = mask_highres(mask, resadd);
x       = (-(10+pad)-dx):dx:((10+pad)+dx);
[X,Y]   = meshgrid(x, x);


data = wfield( mask_hr, Ndata, 'T', 3);

%[min(pars(1,1) * exp(pars(1,2) * xx)), max(pars(1,1) * exp(pars(1,2) * xx))] 
%[min(pars(2,1) * exp(pars(2,2) * xx)), max(pars(2,1) * exp(pars(2,2) * xx))] 

val1 = exp( -(x - Y).^2 ./ (sqrt(2) * pars(1,1) * exp(pars(1,2) * x)).^2) * dx;
val2 = exp( -(x - Y).^2 ./ (sqrt(2) * pars(2,1) * exp(pars(2,2) * x)).^2) * dx;

sd_val = sqrt(pi) ./ sqrt( (pars(2,1) * pars(1,1)) * (exp(pars(1,2) * X)  .* exp(pars(2,2) * Y) )) * pi;

tic
for k = 1:data.fibersize
    data.field(:,:,k) = val1' * data.field(:,:,k) * val2 ./ sd_val;
end
ende = length(data.xvals{1});

data = data((pad + pad*resadd+1):(ende-((pad + pad*resadd))), (pad + pad*resadd+1):(ende-((pad + pad*resadd))),:);
data = data./std(data);
eps  = data.field;
L = [12.53, 39.25];

save(outputname, 'eps', 'L','-v7.3')
toc