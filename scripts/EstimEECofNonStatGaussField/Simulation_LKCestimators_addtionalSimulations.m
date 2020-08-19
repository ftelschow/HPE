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
Nvec     = [10 50 100 200];         % sample sizes
Nsim     = 1e3;                     % number of simulations
casNums  = [1 2]; % which cases for standardization ['theory' 'demeanvar1' 'var1' 'demean']

%-- parameters for methods
% bHPE
Mboot    = 5e3; % number of bootstrap replicates

%-- load parameters from simulated fields
load(strcat(path_data,'/RandomFields_params.mat'))

%% load the isotropic random fields, results for other types of fields can 
% be obtained by choosing TYPE(2) or TYPE(3) in the next line
for type = TYPE(1)
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
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%---------- Simulation of bHPE (using multinomial multipliers) -------
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We compare LKC estimation using the bHPE and different multipliers. We
% compare the Gaussian, Rademacher and two different multinomial multipliers
%--------------------------------------------------------------------------
% set seed to make sure the output is always the same
rng(42)

% This simulations runs in approximate 3 hours on a standard laptop
% initialize containers for the LKCs
 LKChermB  = struct( 'hatn', zeros(2, Nsim, length(Nvec), length(casNums)), ...
                    'hatmean', zeros(2, length(Nvec), length(casNums)),...
                    'hatstd', zeros(2, length(Nvec), length(casNums)) );
 LKChermBn = struct( 'hatn', zeros(2, Nsim, length(Nvec), length(casNums)), ...
                    'hatmean', zeros(2, length(Nvec), length(casNums)),...
                    'hatstd', zeros(2, length(Nvec), length(casNums)) );
 LKChermBm = struct( 'hatn', zeros(2, Nsim, length(Nvec), length(casNums)), ...
                     'hatmean', zeros(2, length(Nvec), length(casNums)),...
                     'hatstd', zeros(2, length(Nvec), length(casNums)) );
 LKChermBr = struct( 'hatn', zeros(2, Nsim, length(Nvec), length(casNums)), ...
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
                case 3 % varia1nce known case, demean data
                    eps1 = standardize(eps(:,:,randsample(1:Ndata, N)), 1, 0);
                case 4 % mean known case, standardize samp-variance to 1
                    eps1 = standardize(eps(:,:,randsample(1:Ndata, N)), 0, 1);
            end

            % Estimate the LKCs using Gaussian multipliers
            tmp = LKCestim_HermProjExact( eps1, D, -666, Mboot, ...
                                          1, "Gaussian", 4, "C" );
            LKChermB.hatn(:,n,i,cont_cas)  = tmp.hatn;

            % Estimate the LKCs using naive Multinomial multipliers
            tmp = LKCestim_HermProjExact( eps1, D, -666, Mboot, ...
                                          1, "Naive", 4, "C" );
            LKChermBn.hatn(:,n,i,cont_cas)  = tmp.hatn;


            % Estimate the LKCs using multinomial multipliers
            tmp = LKCestim_HermProjExact( eps1, D, -666, Mboot, ...
                                          1, "Multinomial", 4, "C" );
            LKChermBm.hatn(:,n,i,cont_cas)  = tmp.hatn;

            % Estimate the LKCs using rademacher multipliers
            tmp = LKCestim_HermProjExact( eps1, D, -666, Mboot, ...
                                          1, "Rademacher", 4, "C" );
            LKChermBr.hatn(:,n,i,cont_cas)  = tmp.hatn;
        end
    end        

    LKChermB.hatmean(:,i,:,:) = mean( LKChermB.hatn(:,:,i,:), 2 );
    LKChermB.hatstd(:,i,:,:)  = std( LKChermB.hatn(:,:,i,:), 0, 2 );

    LKChermBn.hatmean(:,i,:,:) = mean( LKChermBn.hatn(:,:,i,:), 2 );
    LKChermBn.hatstd(:,i,:,:)  = std( LKChermBn.hatn(:,:,i,:), 0, 2 );

    LKChermBm.hatmean(:,i,:,:) = mean( LKChermBm.hatn(:,:,i,:), 2 );
    LKChermBm.hatstd(:,i,:,:)  = std( LKChermBm.hatn(:,:,i,:), 0, 2 );

    LKChermBr.hatmean(:,i,:,:) = mean( LKChermBr.hatn(:,:,i,:), 2 );
    LKChermBr.hatstd(:,i,:,:)  = std( LKChermBr.hatn(:,:,i,:), 0, 2 );
    %%%% save results
    save( strcat( path_results,'simLKChermBm_maxN', num2str(max(Nvec)),...
                  outputname), 'dataname', 'LKChermB', 'LKChermBm',...
                  'LKChermBn', 'LKChermBr', 'L', 'Nvec' )

    toc
end

clear LKChermBm LKChermBr eps1 cont_cas N i


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%---------- Simulation of HPE to show dependence on connectivity ----------
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We compare LKC estimation using the HPE using 4-connectivity and 8-
% connectivity for our standard example of a isotropic Gaussian field.
%--------------------------------------------------------------------------
% set seed to make sure the output is always the same
rng(23)

% This simulations runs in approximate 7 hours on a standard laptop
% initialize containers for the LKCs
LKCherm4  = struct( 'hatn', zeros(2, Nsim, length(Nvec), length(casNums)), ...
                    'hatmean', zeros(2, length(Nvec), length(casNums)),...
                    'hatstd', zeros(2, length(Nvec), length(casNums)) );
LKCherm8 = struct( 'hatn', zeros(2, Nsim, length(Nvec), length(casNums)), ...
                   'hatmean', zeros(2, length(Nvec), length(casNums)),...
                   'hatstd', zeros(2, length(Nvec), length(casNums)) );

%%%% Estimate LKC using hermite estimator
for i = 1:length(Nvec)
    tic
    N = Nvec(i);
    for n = 1:Nsim
        for cont_cas = casNums(1)
            % standardize the data
            switch cont_cas
                case 1 % theory case, no standardization
                    eps1 = eps(:,:,randsample(1:Ndata, N));
                case 2 % experiment case, demean and standardize samp-variance to 1
                    eps1 = standardize(eps(:,:,randsample(1:Ndata, N)), 1, 1);
                case 3 % varia1nce known case, demean data
                    eps1 = standardize(eps(:,:,randsample(1:Ndata, N)), 1, 0);
                case 4 % mean known case, standardize samp-variance to 1
                    eps1 = standardize(eps(:,:,randsample(1:Ndata, N)), 0, 1);
            end

            % Estimate the LKCs using Gaussian multipliers
            tmp = LKCestim_HermProjExact( eps1, D, -666, 1, ...
                                          1, "NaN", 4, "C" );
            LKCherm4.hatn(:,n,i,cont_cas)  = tmp.hatn;

            % Estimate the LKCs using naive Multinomial multipliers
            tmp = LKCestim_HermProjExact( eps1, D, -666, 1, ...
                                          1, "NaN", 8, "C" );
            LKCherm8.hatn(:,n,i,cont_cas)  = tmp.hatn;
        end
    end        

    LKCherm4.hatmean(:,i,:,:) = mean( LKCherm4.hatn(:,:,i,:), 2 );
    LKCherm4.hatstd(:,i,:,:)  = std( LKCherm4.hatn(:,:,i,:), 0, 2 );

    LKCherm8.hatmean(:,i,:,:) = mean( LKCherm8.hatn(:,:,i,:), 2 );
    LKCherm8.hatstd(:,i,:,:)  = std( LKCherm8.hatn(:,:,i,:), 0, 2 );

    %%%% save results
    save( strcat( path_results,'sim_varyingCC_LKCherm_maxN', num2str(max(Nvec)),...
                  outputname), 'dataname', 'LKCherm4', 'LKCherm8', 'L',...
                  'Nvec' )

    toc
end

clear LKCherm4 LKCherm8 eps1 cont_cas N i

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%--------- Simulation of HPE to show dependence on discetization ----------
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This simulation compares the sensitivity of the .
%--------------------------------------------------------------------------
% set seed to make sure the output is always the same
rng(666)

resaddVec = [ 0 1 3 5 7 ];
drVec     = [ 6 4 3 2 1 ];

T      = 49;
Nsim   = 1e3;
N      = 100;
    
dimp  = [145 145];
sigma = 15;

% Precompute random fields
Ndata = 1e3;
tic
[ eps, L ] = generateField( Ndata, 1, dimp, "isotropic", sigma );
toc

% This simulations runs in approximate 7 hours on a standard laptop
% initialize containers for the LKCs
LKCherm4  = struct( 'hatn', zeros( 2, Nsim, length(resaddVec) ), ...
                    'hatmean', zeros( 2, length(resaddVec) ),...
                    'hatstd', zeros( 2, length(resaddVec) ) );
LKCherm8 = struct( 'hatn', zeros( 2, Nsim, length(resaddVec) ), ...
                   'hatmean', zeros( 2, length(resaddVec) ),...
                   'hatstd', zeros( 2, length(resaddVec) ) );

%%%% Estimate LKC using hermite estimator
tic
for n = 1:Nsim
    % Generate random field at high resolution
    eps1 = eps( :, :, randsample( 1:Ndata, N ) );
    
    for i = 1:length(resaddVec)
        dr = drVec(i);
        eps11 = eps1( 1:dr:end, 1:dr:end, : );
        
        % Estimate the LKCs using Gaussian multipliers
        tmp = LKCestim_HermProjExact( eps11, D, -666, 1, ...
                                      1, "NaN", 4, "C" );
        LKCherm4.hatn(:,n,i)  = tmp.hatn;

        % Estimate the LKCs using naive Multinomial multipliers
        tmp = LKCestim_HermProjExact( eps11, D, -666, 1, ...
                                      1, "NaN", 8, "C" );
        LKCherm8.hatn(:,n,i)  = tmp.hatn;
    end        
end
toc

for i = 1:length(resaddVec)
    LKCherm4.hatmean(:,i) = mean( LKCherm4.hatn(:,:,i), 2 );
    LKCherm4.hatstd(:,i)  = std( LKCherm4.hatn(:,:,i), 0, 2 );

    LKCherm8.hatmean(:,i) = mean( LKCherm8.hatn(:,:,i), 2 );
    LKCherm8.hatstd(:,i)  = std( LKCherm8.hatn(:,:,i), 0, 2 );
end

%%%% Plot the resolution figures
scale=0.9;
for dr = drVec( [ 1 2 4 5 ] )
    figure, clf
        WidthFig = scale*500;
        HeightFig = scale*450;
            set(groot, 'defaultAxesTickLabelInterpreter','latex');
            set(groot, 'defaultLegendInterpreter','latex');
            set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
            set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig])
            I = flip( eps( 1:dr:end, 1:dr:end, 1 ), 1 );
            imagesc( I );
            axis equal tight
            colorbar;
            h = title( strcat("Resolution ", num2str(round(1/dr, 3)) ) );
            set(h, 'Interpreter', 'latex');
            set(gca, 'FontSize', 20)
            set(gca,'Ydir','reverse')
            set(gcf,'papersize',[12 12])
%             xticks([10 20 30 40])
%             xticklabels( {'10' '20' '30' '40'} )
%             yticks([0 10 20 30 40])
%             yticklabels( {'50' '40' '30' '20' '10'} )
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
     print(strcat(path_pics,'varyingResolution_LKCherm_res', num2str(dr) ), '-dpng')
end

%%%% save results
save( strcat( path_results, 'sim_varyingResolution_LKCherm.mat'),...
              'LKCherm4', 'LKCherm8', 'L', 'Nvec', 'drVec' )

clear LKCherm4 LKCherm8 eps1 cont_cas N i