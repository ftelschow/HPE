%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%  This script plots the fields in
%%%%  Telschow, Fabian, et al.
%%%%  "Estimation of Expected Euler Characteristic Curves of
%%%%   Nonstationary Smooth Gaussian Random Fields."
%%%%  arXiv preprint arXiv:1908.02493 (2020).
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Fabian Telschow, Armin Schwartzman
%--------------------------------------------------------------------------
%------ prepare workspace
% clean history
clear
close all

% create matlab paths to search for functions
addpath(genpath('/home/drtea/matlabToolboxes/HPE'));
addpath(genpath("/home/fabian/Seafile/Code/matlabToolboxes/RFTtoolbox"))

% load generated profile containing the paths
load('paths.mat')

%-- load parameters from simulated fields
load(strcat(path_data,'/RandomFields_params.mat'))

%--------------------------------------------------------------------------
% define parameters for this script
TYPE      =  "nonstatgauss_exp"; % "isotropic"; % "scale-space"; %   %  ["nonstatnongauss_exp",];
Msim      = 1000;         % number of simulations
Nsubj     = [10 50 100] ; % number of subjects/sample size
FWERlevel = 0.05;         % 1-confidence level
u = -6:0.005:6;            % evaluation grid

%--------------------------------------------------------------------------
% Load precomputed fields
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
        case "ng_scale"
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
        case "nonstatnongauss_exp"
            dataname   = strcat( path_data, 'RandomFields_', 'Ndata', num2str(Ndata),'_',...
                     type, '_D', num2str(D), 'T21', '.mat' );
            outputname = strcat( type, '_D', num2str(D), 'T21');
            T = 43;
        case "nonstatgauss_exp"
            dataname   = strcat( path_data, 'RandomFields_', 'Ndata', num2str(Ndata),'_',...
                     type, '_D', num2str(D), 'T21', '.mat' );
            outputname = strcat( type, '_D', num2str(D), 'T21');
            T = 43;
    end
    
    % load precomputed data
    load(dataname)
    if size(L,1) == 1
        L = L';
    end
end
                 

%--------------------------------------------------------------------------
% Global figure settings
sfont = 20;
addf  = 5;
scale = 3.5/12;

% color scheme for colorblind
Vibrant    = [[0 119 187];... % blue
              [51 187 238];...% cyan
              [0 153 136];... % teal
              [238 119 51];...% orange
              [204 51 17];... % red
              [238 51 119];...% magenta
              [187 187 187]...% grey
              ]/255;
% scale for grey lines
greyScale = 0.55;


%--------------------------------------------------------------------------
%% prepare the random fields for the simulation
% Transform array to field Class
field =  Field(eps, 2);

% Get the EC curves of the data
ECs   = ECcurve( field, u);

% estimate the LKCs using the HPE
LKCs = LKC_HP_est(field, 0, 0);

% true EEC
[ ~, EECtrue] = EEC( u, L', 1, 'Z');

% get "true" standard errors
av_EEC_se_hat_true = std(ECs).field / sqrt(Ndata);
[ ~, ~, HPE_EEC_se_hat_true, ~ ]  =  EEC( u, LKCs.hatL1', 1, 'Z');

%--------------------------------------------------------------------------
%% Simulate pointwise and simultaneous confidence bands
rng(42)

coverage_ptw_av  = zeros([length(u), Msim, length(Nsubj)]);
coverage_ptw_HPE = zeros([length(u), Msim, length(Nsubj)]);
coverage_scb_HPE = zeros([length(u), Msim, length(Nsubj)]);
coverage_ptw_av_true  = zeros([length(u), Msim, length(Nsubj)]);
coverage_ptw_HPE_true = zeros([length(u), Msim, length(Nsubj)]);
coverage_scb_HPE_true = zeros([length(u), Msim, length(Nsubj)]);

for n = 1:length(Nsubj)
    Nsample = Nsubj(n);
    % Get the correct quantiles
    z_a        = norminv(1-FWERlevel/2);
    z_a_finite = tinv(1-FWERlevel/2, Nsample-1);
    q_a        = sqrt(chi2inv( 1 - FWERlevel, 2));
    q_a_finite = sqrt(2*(Nsample-1)*finv(1 - FWERlevel, 2, Nsample-2) / (Nsample-2));
    for m = 1:Msim
        % Get a sample
        sample_vec = randsample(1:field.fibersize, Nsample);

        % Get the HPE EEC estimator
        LKCm     = LKCs.hatL1(:,sample_vec);
        hatSIGMA = cov(LKCm');
        [ ~, HPE_EEC, HPE_EEC_se_hat, C_hat ]  =  EEC( u, LKCm', 1, 'Z');

        % Get the average EEC estimator
        ECsm          = ECs(:,sample_vec);
        av_EEC        = mean(ECsm).field;
        av_EEC_se_hat = std(ECsm).field / sqrt(Nsample);

        % Get the error fields for the two estimators
        av_Error       = abs(av_EEC - EECtrue.field) ./  av_EEC_se_hat;
        HPE_Error      = abs(HPE_EEC - EECtrue) ./ HPE_EEC_se_hat;
        av_Error_true  = abs(av_EEC - EECtrue.field) ./  (av_EEC_se_hat_true * sqrt(Ndata) / sqrt(Nsample));
        HPE_Error_true = abs(HPE_EEC - EECtrue) ./ (HPE_EEC_se_hat_true * sqrt(Ndata) / sqrt(Nsample));

        % Get the coverage
        coverage_ptw_av(:,m,n)  = av_Error <= z_a_finite;
        coverage_ptw_HPE(:,m,n) = HPE_Error.field <= z_a_finite;
        coverage_scb_HPE(:,m,n) = HPE_Error.field <= q_a_finite;
        coverage_ptw_av_true(:,m,n)  = av_Error_true <= z_a;
        coverage_ptw_HPE_true(:,m,n) = HPE_Error_true.field <= z_a;
        coverage_scb_HPE_true(:,m,n) = HPE_Error_true.field <= q_a;
    end
end

%%%% save results
save( strcat( path_results, 'sim_ConfidenceBands_', outputname),...
          'coverage_ptw_av','coverage_ptw_HPE', 'coverage_scb_HPE',...
          'coverage_ptw_av_true','coverage_ptw_HPE_true', 'coverage_scb_HPE_true',...
          'Msim', 'Nsubj', 'q_a', 'z_a', 'type', 'FWERlevel', 'u')


