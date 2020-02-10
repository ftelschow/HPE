%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%    Data analysis cosmology, studying the CMBR
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:  This script provides the data analysis of Cosmic Microwave 
%               Background Radiation presented in Telschow et al (2020).
% REQUIREMENTS: 
%__________________________________________________________________________
% REFERENCES:
%   - Telschow et al, (2020)
%     "Estimation of Expected Euler Characteristic Curves of
%      Nonstationary Smooth Gaussian Random Fields."
%      arXiv preprint arXiv:1908.02493 (2020).
%__________________________________________________________________________
% AUTHOR: Fabian Telschow (ftelschow@ucsd.edu), Armin Schwarzman
%__________________________________________________________________________
%------ prepare workspace
% clean history
clear
close all
% load generated profile containing the paths
load('paths.mat')

% create matlab paths to search for functions
addpath(genpath('/home/drtea/matlabToolboxes/HPE'));

% color scheme for colorblind
Vibrant    = [[0 119 187];... % blue
              [51 187 238];...% cyan
              [0 153 136];... % teal
              [238 119 51];...% orange
              [204 51 17];... % red
              [238 51 119];...% magenta
              [187 187 187]...% grey
              ]/255;

% Load data (created and provided by Pratyush Prana)
filename = strcat(path_data, 'ffp8_smth3_g256_fullDgm_masked.mat');
load(filename);

% Parameters
topo_type = {'betti0', 'betti1', 'EC'};
cmb_type = {'nilc'};
du = 0.01;
u = -5:du:5;   % excursion thresholds

N_topo = length(topo_type);
N_cmb = length(cmb_type);
L_u = length(u);
n = size(topo_sim, 3);

% Extract EC
EC_cmb = -single(topo_cmb(:,3));
EC_sim = -single(squeeze(topo_sim(:,3,:)));     % Minus sign probably typo


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%--------- Analysis of CMBR using the HPE and the nonparametric EEC -------
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We compare nonparametric EEC estimation and the HPE of EEC in order to
% draw conclusions about the physical model of the CMBR. It turns out that
% the observed CMBR is an unlikely event, but physically not significant,
% since 5*sigma is the standard and we have only 2*sigma evidence against
% it.
%__________________________________________________________________________
% estimate LKC, EEC from the real CMBR
L0 = EC_cmb(1);
[LKC_cmb, EEC_cmb, ~] = HermiteEEC(EC_cmb, u, 2, [0; 2], L0);
LKC_cmb.hatn
LKC_cmb.se_hat

%------ HPE method
[LKC_sim, EEC_sim, rho] = HermiteEEC(EC_sim, u, 2, [0; 2], L0);
LKC_sim.hatn
LKC_sim.se_hat
std(LKC_sim.hat1, 1, 2)
% pointwise prediction bands
EEC_sm_pred1_hat = cat(2, EEC_sim.hatn - sqrt(n)*EEC_sim.se_hat,...
                         EEC_sim.hatn + sqrt(n)*EEC_sim.se_hat );
EEC_sm_pred2_hat = cat(2, EEC_sim.hatn - 2*sqrt(n)*EEC_sim.se_hat,...
                         EEC_sim.hatn + 2*sqrt(n)*EEC_sim.se_hat );
EEC_sm_pred3_hat = cat(2, EEC_sim.hatn - 3*sqrt(n)*EEC_sim.se_hat,...
                         EEC_sim.hatn + 3*sqrt(n)*EEC_sim.se_hat );

%------ Nonparametric method
EC_bar = squeeze(mean(EC_sim, 2));
EC_bar_se_hat = sqrt(diag(cov(EC_sim(:,:)')) / n);
% pointwise prediction bands
EC_bar_pred2_hat = cat(2, EC_bar - 2*sqrt(n)*EC_bar_se_hat,...
                         EC_bar + 2*sqrt(n)*EC_bar_se_hat );
EC_bar_pred3_hat = cat(2, EC_bar - 3*sqrt(n)*EC_bar_se_hat,...
                         EC_bar + 3*sqrt(n)*EC_bar_se_hat );

% Inference threshold
FWER = 0.05;
CER = 1;
[ind, u0] = crossing(EEC_sim.hatn, u, FWER, 'linear');
u_FWER = max(u0)
tau = [0; diff(EEC_sim.hatn)]/du;
u_FWER_se = sqrt(EEC_sim.cov(max(ind), max(ind))/(n * tau(max(ind))^2))

[ind, u0] = crossing(EEC_sim.hatn, u, CER, 'linear');
u_CER = max(u0)
u_CER_se = sqrt(EEC_sim.cov(max(ind), max(ind))/(n * tau(max(ind))^2))


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%------------- Plot EEC figures containing the prediction bands -----------
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We compare nonparametric EEC estimation and the HPE of EEC in order to
% draw conclusions about the physical model of the CMBR. It turns out that
% the observed CMBR is an unlikely event, but physically not significant,
% since 5*sigma is the standard and we have only 2*sigma evidence against
% it.
%__________________________________________________________________________
WidthFig  = 800;
HeightFig = 550;

figure(1), clf, hold on, set(gca, 'FontSize', 25)
% Define size and location of the figure [xPos yPos WidthFig HeightFig]
set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig])
set(gcf, 'defaultAxesTickLabelInterpreter','latex');
set(gcf, 'defaultLegendInterpreter','latex');

h = title('EC curves');
set(h, 'Interpreter', 'latex');

plot(u, EC_sim, 'Color', 0.85*[1 1 1], 'LineWidth', 0.5)
plot(u, EC_bar, 'Color', Vibrant(1,:), 'LineWidth', 2.5)
plot(u, EC_bar_pred2_hat, 'LineStyle', '--', 'LineWidth', 2.5, 'Color', Vibrant(1,:))
plot(u, EC_bar_pred3_hat, 'LineStyle', '--', 'LineWidth', 2.5, 'Color', Vibrant(2,:))
plot(u, EC_cmb, 'Color', Vibrant(5,:), 'LineWidth', 1.5)

h = xlabel('u');
set(h, 'Interpreter', 'latex');
xticks([-5 -2.5 0 2.5 5])
xticklabels( {'-5', '-2.5', '0', '2.5', '5'} )
h = ylabel('EC');
yticks([-150 -100 -50 0 50 100 150])
yticklabels( {'-150' '-100' '-50' '0' '50' '100' '150'} )
set(h, 'Interpreter', 'latex');
axis([-5 5 -170 170])
hold off

figure(2), clf, hold on, set(gca, 'FontSize', 25)
% Define size and location of the figure [xPos yPos WidthFig HeightFig]
set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig])
set(gcf, 'defaultAxesTickLabelInterpreter','latex');
set(gcf, 'defaultLegendInterpreter','latex');

h = title('EC curves');
set(h, 'Interpreter', 'latex');

title('Smoothed EC curves')
plot(u, EEC_sim.hat1, 'color', 0.85*[1 1 1])
plot(u, EEC_sim.hatn, 'Color', Vibrant(1,:), 'LineWidth', 2.5)
plot(u, EEC_sm_pred2_hat, 'LineStyle', '--', 'LineWidth', 2.5, 'Color', Vibrant(1,:))
plot(u, EEC_sm_pred3_hat, 'LineStyle', '--', 'LineWidth', 2.5, 'Color', Vibrant(2,:))
plot(u, EEC_cmb.hat1, 'Color', Vibrant(5,:), 'LineWidth', 1.5)

h = xlabel('u');
set(h, 'Interpreter', 'latex');
xticks([-5 -2.5 0 2.5 5])
xticklabels( {'-5', '-2.5', '0', '2.5', '5'} )
h = ylabel('EC');
yticks([-150 -100 -50 0 50 100 150])
yticklabels( {'-150' '-100' '-50' '0' '50' '100' '150'} )
set(h, 'Interpreter', 'latex');
axis([-5 5 -170 170])
hold off