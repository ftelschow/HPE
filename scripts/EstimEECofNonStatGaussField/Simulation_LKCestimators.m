%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%  Simulation of different LKC estimators in 2D
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: This script provides the LKC simulation for Telschow (2020).
%              It contains a comparison of the LKC estimation in a theoretical
%              (mean known to be 0 and variance known to be 1) and the
%              usual experimental setting, where these are unknown.
%              Moreover, it shows that LKC estimation is "unbiased" as
%              smoothness increases.
%__________________________________________________________________________
% REFERENCES:
%   - Telschow et al, (2020)
%     "Estimation of Expected Euler Characteristic Curves of
%      Nonstationary Smooth Gaussian Random Fields."
%      arXiv preprint arXiv:1908.02493 (2020).
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

%------ define parameters for this script
%-- global parameters for the simulation
method   = [0 0 1 0]; % [ "HPE" "bHPE" "warp" "worsley" ]; % LKC estimation method
Nvec     = [10 30 50 75 100 150 200];         % sample sizes
Nsim     = 1e3;                               % number of simulations
casNums  = [1 2]; % which cases for standardization ['theory' 'demeanvar1' 'var1' 'demean']

%-- parameters for methods
% bHPE
Mboot    = 5e3; % number of bootstrap replicates

%-- load parameters from simulated fields
load(strcat(path_data,'/RandomFields_params.mat'))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%---------- Simulation of LKC estimation using different estimators -------
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We compare LKC estimation from Gaussian, non-gaussian and non-isotropic
% data with our estimator to currently used estimators.
% Note that all simulations are under 10 minutes except for the bHPE. It 
% requires to estimate EC curves for each bootstrap sample, which makes it
% slower despite C++ implementation.
%--------------------------------------------------------------------------
% set seed to make sure the output is always the same
rng(23)

for type = TYPE
    switch type
        case "isotropic"
            dataname   = strcat( path_data, '/RandomFields_', 'Ndata', num2str(Ndata),'_',...
                     type, '_D', num2str(D), 'T', num2str(T),...
                    '_params', num2str(nu), '.mat' );
            outputname = strcat( type, '_D', num2str(D), 'T', num2str(T),...
                    '_params', num2str(nu));
        case "scale-space"
            dataname   = strcat( path_data, '/RandomFields_', 'Ndata', num2str(Ndata),'_',...
                     type, '_D', num2str(D), 'T', num2str(T),...
                    '_params', num2str(gamma(1)),...
                    '_', num2str(gamma(end)), '.mat' );
            outputname = strcat( type, '_D', num2str(D), 'T', num2str(T),...
                    '_params', num2str(gamma(1)), '_', num2str(gamma(end)));
        case "nongauss"
            dataname   = strcat( path_data, '/RandomFields_', 'Ndata', num2str(Ndata),'_',...
                     type, '_D', num2str(D), 'T', num2str(T),...
                    '_params', num2str(nu), '.mat' );
            outputname = strcat( type, '_D', num2str(D), 'T', num2str(T),...
                    '_params', num2str(nu));
    end
    
    % load precomputed data
    load(dataname)

%% %%%% HPE    
if method(1) % This simulations run in under 10 minutes
    % initialize containers for the LKCs
    LKCherm     = struct( 'hatn', zeros(2, Nsim, length(Nvec), length(casNums)), ...
                          'hatmean', zeros(2, length(Nvec), length(casNums)),...
                          'hatstd', zeros(2, length(Nvec), length(casNums)) );
    %%%% Estimate LKC using hermite estimator
    for i = 1:length(Nvec)
        tic
        N = Nvec(i);
        for n = 1:Nsim
            for cont_cas = casNums
                % standardize the data
                switch cont_cas
                    case 1 % theory case, no standardization
                        eps1 = eps(:,:,randsample(1:Ndata, N));
                    case 2 % experiment case, demean and standardize samp-variance to 1
                        eps1 = standardize(eps(:,:,randsample(1:Ndata, N)), 1, 1);
                    case 3 % variance known case, demean data
                        eps1 = standardize(eps(:,:,randsample(1:Ndata, N)), 1, 0);
                    case 4 % mean known case, standardize samp-variance to 1
                        eps1 = standardize(eps(:,:,randsample(1:Ndata, N)), 0, 1);
                end
                % Estimate the LKC
                tmp = LKCestim_HermProjExact( eps1, D, -666, 1, ...
                                              1, "Gaussian", 4, "C" );
                LKCherm.hatn(:,n,i,cont_cas)  = tmp.hatn;
            end
        end
        LKCherm.hatmean(:,i,:) = mean( LKCherm.hatn(:,:,i,:), 2 );
        LKCherm.hatstd(:,i,:)  = std( LKCherm.hatn(:,:,i,:), 0, 2 );
        %%%% save results
        save( strcat( path_results,'/simLKCherm_maxN', num2str(max(Nvec)),...
                      outputname), 'dataname', 'LKCherm', 'L', 'Nvec' )

        toc
    end

    clear LKCherm eps1 cont_cas N i
end

%% %%%% bHPE
if method(2) % This simulations runs in approximate 7 hours on a standard laptop
    % initialize containers for the LKCs
     LKChermB     = struct( 'hatn', zeros(2, Nsim, length(Nvec), length(casNums)), ...
                          'hatmean', zeros(2, length(Nvec), length(casNums)),...
                          'hatstd', zeros(2, length(Nvec), length(casNums)) );
    %%%% Estimate LKC using hermite estimator
    for i = 1:length(Nvec)
        tic
        N = Nvec(i);
        for n = 1:Nsim
            for cont_cas = casNums
                % standardize the data
                switch cont_cas
                    case 1 % theory case, no standardization
                        eps1 = eps(:,:,randsample(1:Ndata, N));
                    case 2 % experiment case, demean and standardize samp-variance to 1
                        eps1 = standardize(eps(:,:,randsample(1:Ndata, N)), 1, 1);
                    cdcase 3 % variance known case, demean data
                        eps1 = standardize(eps(:,:,randsample(1:Ndata, N)), 1, 0);
                    case 4 % mean known case, standardize samp-variance to 1
                        eps1 = standardize(eps(:,:,randsample(1:Ndata, N)), 0, 1);
                end

                % Estimate the LKC
                tmp = LKCestim_HermProjExact( eps1, D, -666, Mboot, ...
                                              1, "Gaussian", 4, "C" );
                LKChermB.hatn(:,n,i,cont_cas)  = tmp.hatn;
            end
        end
        LKChermB.hatmean(:,i,:) = mean( LKChermB.hatn(:,:,i,:), 2 );
        LKChermB.hatstd(:,i,:)  = std( LKChermB.hatn(:,:,i,:), 0, 2 );
        %%%% save results
        save( strcat( path_results,'/simLKChermB_maxN', num2str(max(Nvec)),...
                      outputname), 'dataname', 'LKChermB', 'L', 'Nvec' )

        toc
    end

    clear LKChermB eps1 cont_cas N i
end

%% %%%% WarpE
if method(3) % This simulations run in a couple of minutes
    %initialize containers for the LKCs
    LKCwarp    = struct( 'hatn', zeros(2, Nsim, length(Nvec), length(casNums)), ...
                         'hatmean', zeros(2, length(Nvec), length(casNums)),...
                         'hatstd', zeros(2, length(Nvec), length(casNums)) );
    %%%% Estimate LKC using hermite estimator
    for i = 1:length(Nvec)
        tic
        N = Nvec(i);
        for cont_cas = casNums
            % standardize the data
            switch cont_cas
                case 1
                    center    = 0;
                    normalize = 0;
                case 4
                    center    = 1;
                    normalize = 0;
                case 3
                    center    = 0;
                    normalize = 1;
                case 2
                    center    = 1;
                    normalize = 1;
            end
            tmp = zeros(D,Nsim);
            %%%% Estimate LKC using Taylor Worsley
            for m = 1:Nsim
                tmp(:,m) = LKCestim_warp( eps(:,:,randsample(1:Ndata, N)), ...
                                          center, normalize );
            end
            clear m

            % Estimate the LKC
            LKCwarp.hatn(:,:,i,cont_cas)  = tmp;
            LKCwarp.hatmean(:,i,cont_cas) = mean( LKCwarp.hatn(:,:,i,cont_cas), 2 );
            LKCwarp.hatstd(:,i,cont_cas)  = std( LKCwarp.hatn(:,:,i,cont_cas), 0, 2 );
        end
        %%%% save results
        save( strcat( path_results,'/simLKCwarp_maxN', num2str(max(Nvec)),...
                      outputname), 'dataname', 'LKCwarp', 'L', 'Nvec' )

        toc
    end
    clear LKCwarp tmp con_cas N i

end

%% %%%% IsotE
if method(4) % This simulations run in a couple of minutes
    if ~strcmp(type,"scale-space")
        mask = ones([T T]);
    else
        mask = ones([T length(gamma)]);
    end
    % initialize containers for the LKCs
    LKCiso    = struct( 'hatn', zeros(2, Nsim, length(Nvec), length(casNums)), ...
                        'hatmean', zeros(2, length(Nvec), length(casNums)),...
                        'hatstd', zeros(2, length(Nvec), length(casNums)) );
    %%%% Estimate LKC using hermite estimator
    for i = 1:length(Nvec)
        tic
        N = Nvec(i);
        for cont_cas = casNums
            % standardize the data
            switch cont_cas
                case 1
                    typ = [0 0 1];
                    df = N;
                case 2
                    typ = [1 1 1];
                    df = N-1;
                case 3
                    typ = [1 0 1];
                    df = N-1;
                case 4
                    typ = [0 1 1];
                    df = N;
            end
            
            %%%% Estimate LKC using sam's (direct isotropic appraoch)           
            tmp = zeros(D, Nsim);
            for m = 1:Nsim
                % Estimate resels using SPM
                tmp(:,m) = LKCestim_iso( eps(:,:,randsample(1:Ndata, N)),...
                                                   mask, typ, df );
            end
            clear m fwhm_est

            % Estimate the LKC
            LKCiso.hatn(:,:,i,cont_cas)  = tmp;
            LKCiso.hatmean(:,i,cont_cas) = mean( LKCiso.hatn(:,:,i,cont_cas), 2 );
            LKCiso.hatstd(:,i,cont_cas)  = std( LKCiso.hatn(:,:,i,cont_cas), 0, 2 );
        end
        toc
    end

        %%%% save results
        save( strcat( path_results,'/simLKCiso_maxN', num2str(max(Nvec)),...
                      outputname), 'dataname', 'LKCiso', 'L', 'Nvec' )

    clear eps1 con_cas N i
end
toc

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%------------- Simulation of dependence of BIAS on smoothness -------------
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% As discussed in the article the data generation induces a bias compared
% to the theoretical values. We show that for larger smoothing this "bias"
% vanishes. Note especially that the true LKCs for the discrete convolution
% process are different from the theoretical process and therefore the term
% bias is slightly misleading.
% We discuss this effect only in the theoretical framework and only with
% the HPE. But it generalizes to all other estimators as well.
%--------------------------------------------------------------------------
%------ prepare workspace
% clean history
clear
close all
% load generated profile containing the paths
load('paths.mat')

% create matlab paths to search for functions
addpath(genpath(path_main));

cd( path_main )

%------ define parameters for this script
%-- global parameters for the simulation
Nvec     = [10 75];    % sample sizes
Nsim     = 5e2;        % number of simulations

%%%% Parameters for the fields
% Field parameters
L0    = 1;           % EC of domain
D     = 2;           % Dimension of domain (including scale)
T     = 50;          % size of domain
dim   = [T T];
nuVec = [2 3 5 6 7]; % parameter for isotropic

% Outputname
outputname = strcat( 'D', num2str(D),'T', num2str(T),...
                     'Nsim', num2str(Nsim), '_params', num2str(nuVec(1)),...
                     '_', num2str(nuVec(end)) );
                 
% compute the true LKCs
L1 = sqrt(2)./nuVec.*sum((dim-1))/2;
L2 = 1./(2*nuVec.^2)*prod(dim-1);

L  = [L1; L2];

%% ------- Simulate the dependence on smoothness using HPE and WarpE ------
% set seed to make sure the output is always the same
rng(666)
%loop twice to save working memory
for nbatch = 1:2
    % initialize containers for the LKCs
    LKCherm = struct( 'hatn', zeros(2, Nsim, length(Nvec), length(nuVec)), ...
                      'hatmean', zeros(2, length(Nvec), length(nuVec)),...
                      'hatstd', zeros(2, length(Nvec), length(nuVec)) );  
    LKCwarp = LKCherm;
    %%%% Estimate LKC using hermite estimator
    for nu_i = 1:length(nuVec)
        tic
        % generate data
        f = generateField( max(Nvec), Nsim, dim, "isotropic", nuVec(nu_i) );

        for i = 1:length(Nvec)
            N    = Nvec(i);

            % Estimate the LKC
            tmp  = LKCestim_HermProjExact( squeeze(f(:,:,1:N,:)), 2, -666, 1,...
                                           [-1 1], 1, "Gaussian", 4, "C");
            LKCherm.hatn(:,:,i,nu_i)  = tmp.hatn;
            LKCherm.hatmean(:,i,nu_i) = mean( LKCherm.hatn(:,:,i,nu_i), 2 );
            LKCherm.hatstd(:,i,nu_i)  = std( LKCherm.hatn(:,:,i,nu_i), 0, 2 );

            tmp = zeros( D, Nsim );
            %%%% Estimate LKC using Taylor Worsley
            for m = 1:Nsim
                tmp(:,m) = LKCestim_warp( squeeze(f(:,:,1:N,m)), ...
                                          1, 1, 0, 1 );
            end
            clear m

            % Estimate the LKC
            LKCwarp.hatn(:,:,i,nu_i)  = tmp;
            LKCwarp.hatmean(:,i,nu_i) = mean( LKCwarp.hatn(:,:,i,nu_i), 2 );
            LKCwarp.hatstd(:,i,nu_i)  = std( LKCwarp.hatn(:,:,i,nu_i), 0, 2 );
        end
        toc
    end

    if nbatch == 1
        LKCherm1 = LKCherm;
        LKCwarp1 = LKCwarp;
    else
        tmp = [ LKCherm1.hatn LKCherm.hatn ];
        tmp1 = [ LKCwarp1.hatn LKCwarp.hatn ];
        
        LKCherm = struct( 'hatn', zeros(2, Nsim, length(Nvec), length(nuVec)), ...
                      'hatmean', zeros(2, length(Nvec), length(nuVec)),...
                      'hatstd', zeros(2, length(Nvec), length(nuVec)) );  
        LKCwarp = LKCherm;
        LKCherm.hatn    = tmp;
        LKCherm.hatmean = squeeze(mean( tmp, 2 ));
        LKCherm.hatstd  = squeeze(std( tmp, 0, 2 ));
        LKCwarp.hatn    = tmp;    
        LKCwarp.hatmean = squeeze(mean( tmp1, 2 ));
        LKCwarp.hatstd  = squeeze(std( tmp1, 0, 2 ));
        
        %%%% save results
        save( strcat( path_results, 'simSmoothnessDependenceLKC_maxN', ...
                  num2str(max(Nvec)), outputname),...
          'LKCherm','LKCwarp', 'L', 'Nvec', 'nuVec', 'T', 'Nsim' )
    end
end