%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%    Data analysis fMRI - Different LKC estimators
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:  This script provides the application of data analysis of
%               the Moran data presented in Telschow et al (2020).
% REQUIREMENTS: Matlab spm12 Toolbox
%               https://www.fil.ion.ucl.ac.uk/spm/software/spm12/ 
%__________________________________________________________________________
% REFERENCES:
%   - Telschow et al, (2020)
%     "Estimation of Expected Euler Characteristic Curves of
%      Nonstationary Smooth Gaussian Random Fields."
%      arXiv preprint arXiv:1908.02493 (2020).
%__________________________________________________________________________
% AUTHOR: Fabian Telschow (ftelschow@ucsd.edu), Armin Schwartzman
%__________________________________________________________________________
%------ prepare workspace
% clean history
clear
close all
% load generated profile containing the paths
load('paths.mat')

% create matlab paths to search for functions
addpath(genpath(path_main));
addpath(genpath(path_spm12));

%------  Processing choices
% smoothing of data/statistic yes/no
smoothT    = 0;
smoothData = 1;
% bandwidth for smoothing
FWHM       =  1.6*sqrt(2*log(2))*2;

% Inference parameters
FWER = 0.05;
CER = 1;

% Parameters for EC curves
du = 0.01;
u = -8:0.01:8;

% Load Moran data and put into appropriate buckets
load( strcat(path_data,'/sub049_regr_data_fwhm0.mat'))
Y    = double(regr_data.fMRI);
mY   = regr_data.mn_mask; % mean of data within mask
mask = regr_data.mn_mask > 0;
X    = regr_data.X;
c    = regr_data.c;

% invert the order in the z-direction
maxzcoord = size(mask,3);
Y    = Y(:,:,maxzcoord:-1:1,:);
mask = mask(:,:,maxzcoord:-1:1);
mY   = mY(:,:,maxzcoord:-1:1);

% Get constants from the data set
n  = size(X, 1);
sY = size(mask); 
D  = length(sY);

% compute EC of mask
L0 = EulerChar( mask, 0.5, D );

% Color scheme for plots
Vibrant    = [[0 119 187];... % blue
              [51 187 238];...% cyan
              [0 153 136];... % teal
              [238 119 51];...% orange
              [204 51 17];... % red
              [238 51 119];...% magenta
              [187 187 187]...% grey
              ]/255;

% Container for results
HPE = struct('hatLKC', ones([1, D]), 'u', u, 'hatEEC', struct, ...
             'u_FWER', -666, 'u_FWER_se', -666, 'u_CER', -666, 'u_CER_se', -666,...
             'Cluster', [-666 -666]);
bHPE = struct('hatLKC', ones([1, D]), 'u', u, 'hatEEC', struct, ...
             'u_FWER', -666, 'u_FWER_se', -666, 'u_CER', -666, 'u_CER_se', -666,...
             'Cluster', [-666 -666]);
SPM  = struct('hatLKC', ones([1, D]), 'u', u, 'hatEEC', struct, ...
             'u_FWER', -666, 'u_FWER_se', -666, 'u_CER', -666, 'u_CER_se', -666,...
             'Cluster', [-666 -666]);

% set seed to make sure the output is always the same
rng(42)
%% %%%% Preprocessing before data and compute T statistic
% Slice to visualize
nslice = 22;
% plot data
figure
imagesc(Y(:,:,nslice,100)); axis square
colorbar;

% smooth fitts using spm_smooth with Gaussian kernel with specified FWHM
% and plot smoothed data
tic
if smoothData
    Ys = smoothfMRIvolume(Y, FWHM);
else
    Ys = Y;    
end
toc

figure
imagesc(Ys(:,:,nslice,100))
colorbar;

% remove zeros from the boundaries to make pics look more pretty later
cut1 = 11:sY(1)-10;
cut2 = 9:sY(1)-8;
cut3 = 3:36;
nslice = nslice-2;

Ys   = Ys(cut1, cut2, cut3, :);
mask = mask(cut1, cut2, cut3);
mY   = mY(cut1, cut2, cut3);
sY   = size(mask);

% fit GLM to the data, note the the residuals are standardized
[betahat, fitts, residuals, sigma2hat, df, T] = fitGLM2fMRIvolume(Ys, X, c);
clear Ys

% plot a slice of the fit
figure
imagesc(fitts(:,:,nslice,100)); axis square
colorbar;


%% %%%%% mask the residuals and the T-field
tic
if smoothT
    T = smoothfMRIvolume(T, FWHM);
    % get the asymptotical correct residuals for the smoothed Wald
    % statistic
    residuals = smoothfMRIvolume(residuals, FWHM);
    % normalize the T-field to have asymptotically variance 1 and
    % normalize the residuals as well
    hatSigmaSmooth = sqrt( sum(residuals.^2, 4) / (n-1) );
    T = T ./ hatSigmaSmooth;
    for i =1:n
       residuals(:,:,:,i) = residuals(:,:,:,i) ./ hatSigmaSmooth;
    end
else
    T = T;    
end
toc

for i =1:n
    tmp                        = residuals(:,:,:,i);
    tmp(~mask)                 = -Inf;
    residuals(:,:,:,i)         = tmp;
end

T(~mask)  = 0;

figure
imagesc(T(:,:,nslice)); axis square
colorbar;
set(gca,'FontSize',20)
set(findall(gcf,'type','text'),'FontSize',20)

%% Plot the standard deviation map
tmp = sigma2hat;
tmp(~mask) = 0;
figure(5), hold on, set(gcf, 'Position', [ 100 100 1500 1500]) 
subplot(2,3,1), imagesc(tmp(:,:,30)), colorbar;
subplot(2,3,2), imagesc(tmp(:,:,20)), colorbar;
subplot(2,3,3), imagesc(tmp(:,:,25)), colorbar;

% find variance outlier
tmp = sigma2hat(logical(mask) );
Ioutlier = find( sigma2hat>2.5*median(tmp(:)) );
mask2 = mask;
mask2(Ioutlier) = 0;

tmp = sigma2hat;
tmp(~mask2) = 0;
%figure(6), set(gcf, 'Position', [ 100 100 1800 500]) 
subplot(2,3,4), imagesc(tmp(:,:,30)), colorbar;
subplot(2,3,5), imagesc(tmp(:,:,20)), colorbar;
subplot(2,3,6), imagesc(tmp(:,:,25)), colorbar;
hold off
clear mask2 tmp Ioutlier 

% close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%--------------- Estimate the LKCs with different procedures --------------
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---- Estimate LKCs using HPE
% L1 = 13.5, L2 = 261.3, L3= 650.7
hatLherm = LKCestim_HermProjExact( residuals, D, mask, 1, u );
hatLherm.hatn
%---- Estimate LKCs using HPE
% L1 = 13.2, L2 = 266.7, L3= 670.1
% is based of generating a sample)
hatLhermB = LKCestim_HermProjExact( residuals, D, mask, 1e4, u );
hatLhermB.hatn

%---- Estimate LKCs using SPM
% transform data and mask into a NIFTI file in order to process it in SPM
nii_img = make_nii(double(mask) );
save_nii(nii_img, strcat(path_data,'/MoranData_mask.nii') );
clear nii_img

tmp = strcat(path_data, '/tmp');
mkdir(tmp)

for j = 1:n
    nii_img = make_nii( residuals(:,:,:,j) );
    if j<10
        save_nii(nii_img, strcat(path_data, '/tmp/00', num2str(j), 'im.nii' ));
    elseif j<100
        save_nii(nii_img, strcat(path_data, '/tmp/0', num2str(j), 'im.nii' ));
    else
        save_nii(nii_img, strcat(path_data, '/tmp/', num2str(j), 'im.nii' ));
    end
end
clear j

dinfo = dir(strcat(path_data, '/tmp/', '*im.nii'));
names = {dinfo.name};

% Estimate resels using SPM
[FWHM2,~,R] = spm_est_smoothness( strcat(path_data, '/tmp/', char(names)), ...
                                  strcat(path_data,'/MoranData_mask.nii'),...
                                  [n df] );

% remove .nii files and unnecessary variables
clear dinfo names
system( strcat('rm', " ", path_data, '/tmp/', '*.nii' ) );

% Transform resels into LKCs
% L1 = 35.8, L2 = 315.3, L3 = 669.0
SPM.hatLKC = (R(2:D+1) .* [(4*log(2))^(1/2) (4*log(2)) (4*log(2))^(3/2)])';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%--------------- Estimate the EEC and thresholds and inference ------------
%--------------- using the LKC estimates ----------------------------------
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%% Thresholds from LKC estimates
uvals = -6:0.01:6;

%------ HPE
HPE.hatLKC = hatLherm.hatn;
se.N      = size( hatLherm.hat1, 2 );
se.LKCcov = hatLherm.Sigma_hat;

% FWER
[u_FWER, u_FWER_se] = get_EECthreshold( FWER, uvals, HPE.hatLKC, L0,...
                                              "gaussian", se);
HPE.u_FWER    = u_FWER;
HPE.u_FWER_se = u_FWER_se;

% Get number of clusters in excursion set of u_FWER
thresh         = HPE.u_FWER;
T_thresh       = T.*(T >= thresh)./(T >= thresh);
HPE.Cluster(1) = EulerChar(T_thresh, thresh, 3);

% CER
[u_CER, u_CER_se] = get_EECthreshold( CER, uvals, HPE.hatLKC, L0,...
                                              "gaussian", se);
HPE.u_CER    = u_CER;
HPE.u_CER_se = u_CER_se;

% Get number of clusters in excursion set of u_CER
thresh          = HPE.u_CER;
T_thresh        = T.*(T >= thresh)./(T >= thresh);
HPE.Cluster(2) = EulerChar(T_thresh, thresh, 3);


%------ bHPE
bHPE.hatLKC = hatLhermB.hatn; 
se.N        = size( hatLhermB.hat1, 2 );
se.LKCcov   = hatLhermB.Sigma_hat;

% FWER
[u_FWER, u_FWER_se] = get_EECthreshold( FWER, uvals, bHPE.hatLKC, L0,...
                                              "gaussian", se);
bHPE.u_FWER    = u_FWER;
bHPE.u_FWER_se = u_FWER_se;

% Get number of clusters in excursion set of u_FWER
thresh          = bHPE.u_FWER;
T_thresh        = T.*(T >= thresh)./(T >= thresh);
bHPE.Cluster(1) = EulerChar(T_thresh, thresh, 3);

% CER
[u_CER, u_CER_se] = get_EECthreshold( CER, uvals, bHPE.hatLKC, L0,...
                                              "gaussian", se);
bHPE.u_CER    = u_CER;
bHPE.u_CER_se = u_CER_se;

% Get number of clusters in excursion set of u_CER
thresh          = bHPE.u_CER;
T_thresh        = T.*(T >= thresh)./(T >= thresh);
bHPE.Cluster(2) = EulerChar(T_thresh, thresh, 3);


%------ spm12
se = struct();

% FWER
[u_FWER, u_FWER_se] = get_EECthreshold( FWER, uvals, SPM.hatLKC, L0,...
                                              "gaussian", se);
SPM.u_FWER    = u_FWER;
SPM.u_FWER_se = u_FWER_se;

% Get number of clusters in excursion set of u_FWER
thresh          = SPM.u_FWER;
T_thresh        = T.*(T >= thresh)./(T >= thresh);
SPM.Cluster(1)  = EulerChar(T_thresh, thresh, 3);

% CER
[u_CER, u_CER_se] = get_EECthreshold( CER, uvals, SPM.hatLKC, L0,...
                                              "gaussian", se);
SPM.u_CER    = u_CER;
SPM.u_CER_se = u_CER_se;

% Get number of clusters in excursion set of u_CER
thresh         = SPM.u_CER;
T_thresh       = T.*(T >= thresh)./(T >= thresh);
SPM.Cluster(2) = EulerChar(T_thresh, thresh, 3);

save(strcat( path_data, "results_fMRI_smoothdata", num2str(smoothData),...
            "smoothT", num2str(smoothT)), 'u', 'HPE', 'bHPE', 'SPM')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%--------------- Produce figures shown in the manuscript ------------------
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------ EEC Figures
load(strcat( path_data, "results_fMRI_smoothdata", num2str(smoothData),...
            "smoothT", num2str(smoothT)))

tic
EC  = EulerChar(residuals, u, 3);
toc

[~, EEC, ~] = HermiteEEC(EC, u, D, [0; 0; 2], L0);

% Nonparametric method
EC_bar = squeeze(mean(EC, 2));
EC_bar_se_hat = sqrt(diag(cov(EC(:,:)')) / n);
EC_bar_conf_hat = cat(3, EC_bar - 1.96*EC_bar_se_hat, EC_bar + 1.96*EC_bar_se_hat);

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

sfont = 20;
addf  = 5;
scale = 8/12;

xvec       = [-6, -4, -2, 0, 2, 4, 6];
xtickcell  = {'-6', '-4', '-2', '0', '2', '4', '6'};
yvec       = [-100 -40 0 40 75];
ytickcell  = {'-100', '-40', '0', '40', '75'};

% Plot EC curves
figure(5), clf, hold on, set(gca, 'FontSize', sfont)
h = title('EC curves', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
plot(u, EC, 'color', Vibrant(end,:)), plot(u, EC_bar, 'color', Vibrant(1,:), 'LineWidth', 2)
plot(u, EC_bar_conf_hat(:,1,1), '--', 'color', Vibrant(1,:), 'LineWidth', 2),
plot(u, EC_bar_conf_hat(:,1,2), '--', 'color', Vibrant(1,:), 'LineWidth', 2)

h = xlabel('u', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
h = ylabel('EC', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');

% Customize x and y axis
xlim([xvec(1) xvec(end)])
xticks(xvec)
xticklabels(xtickcell)
ylim([yvec(1) yvec(end)])
yticks(yvec)
yticklabels(ytickcell)

    % Print figure to file
    set(gcf,'papersize',[12 12*scale])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(strcat(path_pics,'fMRI_EC_curves.png'), '-dpng')
hold off

figure(6), clf, hold on, set(gca, 'FontSize', sfont)
h = title('Smooth EC curves', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
plot(u, EEC.hat1, 'color', Vibrant(end,:))
plot(u, EEC.hatn, 'color', Vibrant(1,:), 'LineWidth', 2)
plot(u, EEC.conf_hat(:,1,1), '--', 'color', Vibrant(1,:), 'LineWidth', 2)
plot(u, EEC.conf_hat(:,1,2), '--', 'color', Vibrant(1,:), 'LineWidth', 2)
h = xlabel('u', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
h = ylabel('EC', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');

xlim([xvec(1) xvec(end)])
xticks(xvec)
xticklabels(xtickcell)
ylim([yvec(1) yvec(end)])
yticks(yvec)
yticklabels(ytickcell)

    % Print figure to file
    set(gcf,'papersize',[12 12*scale])
    fig = gcf;
    fig.PaperPositionMode = 'auto'
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(strcat(path_pics,'fMRI_EC_curves_smooth.png'), '-dpng')
hold off

% figure(7), clf
% title('Covariance function of smooth EC curves')
% imagesc(u, u, EECherm.cov), axis square tight xy, colorbar


%% Activation map figures
% thresh = 3.9156;    % Height threshold from peakDetRF paper
thresh = bHPE.u_FWER;
%thresh = u_CER;
T_thresh = T.*(T >= thresh)./(T >= thresh);
EulerChar(T_thresh, thresh, 3)

figure(8)
RGB = anatomy(T_thresh(:,end:-1:1,30:-1:6), mY(:,end:-1:1,30:-1:6), [], 'hot');
montage(permute(RGB, [1 2 4 3]))
axis xy
camroll(90)
% thresh = 3.9156;    % Height threshold from peakDetRF paper

thresh = bHPE.u_CER;
T_thresh = T.*(T >= thresh)./(T >= thresh);
EulerChar(T_thresh, thresh, 3)
figure(9)
RGB = anatomy(T_thresh(:,end:-1:1,30:-1:6), mY(:,end:-1:1,30:-1:6), [], 'hot');
montage(permute(RGB, [1 2 4 3]))
axis xy
camroll(90)