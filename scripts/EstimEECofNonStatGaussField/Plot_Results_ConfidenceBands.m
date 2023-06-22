%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%  This script is producing the figures for the confidence band
%%%%  simulations for the HPE paper
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: It simulates the performance of confidence sets for the EEC
%              derived from the plug-in HPE and the simple EC sample
%              average and plots the results.
%__________________________________________________________________________
% REFERENCES:
%   - Telschow et al, (2020)
%     "Estimation of Expected Euler Characteristic Curves of
%      Nonstationary Smooth Gaussian Random Fields."
%      arXiv preprint arXiv:1908.02493 (2020).
%
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
addpath(genpath(path_main));

type =  "nonstatgauss_exp"; %  'isotropic'; % "scale-space"; %

% Global figure settings
sfont = 20;
addf  = 5;
%scale = 3.5/12;

WidthFig = 400;
HeightFig = 350;

HighContr  = [[221, 170,  51];...   % yellow
              [187,  85, 102];...   % red
              [  0,  68, 136]]/255; % blue
Vibrant    = [[0 119 187];... % blue
              [51 187 238];...% cyan
              [0 153 136];... % teal
              [238 119 51];...% orange
              [204 51 17];... % red
              [238 51 119];...% magenta
              [187 187 187]...% grey
              ]/255;

colMat = Vibrant([1 3 4 5],:);

%--------------------------------------------------------------------------
% Load simulation results
D = 2;
T  = 50;
nu = 5;
gamma   = [4 15];

switch type
    case "isotropic"
        outputname = strcat( type, '_D', num2str(D), 'T', num2str(T),...
                '_params', num2str(nu));
    case "scale-space"
        outputname = strcat( type, '_D', num2str(D), 'T', num2str(T),...
                '_params', num2str(gamma(1)), '_', num2str(gamma(end)));
    case "ng_scale"
        outputname = strcat( type, '_D', num2str(D), 'T', num2str(T),...
                '_params', num2str(gamma(1)), '_', num2str(gamma(end)));
    case "nongauss"
        outputname = strcat( type, '_D', num2str(D), 'T', num2str(T),...
                '_params', num2str(nu));
    case "nonstatnongauss_exp"
        outputname = strcat( type, '_D', num2str(D), 'T21');
        T = 43;
    case "nonstatgauss_exp"
        outputname = strcat( type, '_D', num2str(D), 'T21');
        T = 43;
end

% load precomputed data
load(strcat( path_results, 'sim_ConfidenceBands_', outputname))

[mean(all(coverage_scb_HPE_true(:,:,1),1)),...
mean(all(coverage_scb_HPE_true(:,:,2),1)),...
mean(all(coverage_scb_HPE_true(:,:,3),1))]

[mean(all(coverage_scb_HPE(:,:,1),1)),...
mean(all(coverage_scb_HPE(:,:,2),1)),...
mean(all(coverage_scb_HPE(:,:,3),1))]


%--------------------------------------------------------------------------
%% Plot figures for the manuscript
for n = 1:length(Nsubj)
    nsubj = Nsubj(n);
    leg1 = strcat('$\bar\chi^{(', num2str(nsubj),')}$');
    leg2 = strcat('$\widehat{{\rm EEC}}^{(', num2str(nsubj),')}$');
    cov_ptw_HPE      = mean(coverage_ptw_HPE(:,:,n),2);
    cov_ptw_av       = mean(coverage_ptw_av(:,:,n),2);
    cov_ptw_HPE_true = mean(coverage_ptw_HPE_true(:,:,n),2);
    cov_ptw_av_true  = mean(coverage_ptw_av_true(:,:,n),2);

    % Coverage performance true variance
    figure(1), clf, hold on
    WidthFig = 400;
    HeightFig = 350;
        set(groot, 'defaultAxesTickLabelInterpreter','latex');
        set(groot, 'defaultLegendInterpreter','latex');
        set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
        set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig])
        plot(u, cov_ptw_av_true, 'Color', Vibrant(1,:), 'LineWidth', 2)
        plot(u, cov_ptw_HPE_true, 'Color', Vibrant(5,:), 'LineWidth', 2)
        plot([-8 8], [0.95 0.95], 'k--', 'LineWidth', 2)
        axis([u(1) u(end) 0 1])
%        xlabel('u')
        xticks([-5 -2.5 0 2.5 5])
        xticklabels( {'-5', '-2.5', '0', '2.5', '5'} )
        yticks([0 0.2 0.4 0.6 0.8 1])
        yticklabels( {'0' '0.2' '0.4' '0.6' '0.8' '1'} )
        hold off
%        h = xlabel('u');
%        set(h, 'Interpreter', 'latex');
        set(gca, 'FontSize', 20)
        h = legend(leg1, leg2, 'Location', 'south');
        legend('boxoff')
        set(h, 'FontSize', 25)
        set(h, 'Interpreter', 'latex');
        set(gcf,'papersize',[12 12])
        axis tight
        xlim([-6,6])
        ylim([0,1])
        hold off
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(strcat(path_pics, type, "EC_coverage_trueVar_nsim",num2str(Msim), "n", num2str(nsubj),".png"), '-dpng')
    
    
    % Coverage performance estimated variance
    figure(2), clf, hold on
        set(groot, 'defaultAxesTickLabelInterpreter','latex');
        set(groot, 'defaultLegendInterpreter','latex');
        set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
        set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig])
        plot(u, cov_ptw_av, 'Color', Vibrant(1,:), 'LineWidth', 2)
        plot(u, cov_ptw_HPE, 'Color', Vibrant(5,:), 'LineWidth', 2)
        plot([-8 8], [0.95 0.95], 'k--', 'LineWidth', 2)
        axis([u(1) u(end) 0 1])
%        xlabel('u')
        xticks([-5 -2.5 0 2.5 5])
        xticklabels( {'-5', '-2.5', '0', '2.5', '5'} )
        yticks([0 0.2 0.4 0.6 0.8 1])
        yticklabels( {'0' '0.2' '0.4' '0.6' '0.8' '1'} )
        hold off
%        h = xlabel('u');
%        set(h, 'Interpreter', 'latex');
        set(gca, 'FontSize', 20)
        h = legend(leg1, leg2, 'Location', 'south');
        legend('boxoff')
        set(h, 'FontSize', 25)
        set(h, 'Interpreter', 'latex');
        set(gcf,'papersize',[12 12])
        axis tight
        xlim([-6,6])
        ylim([0,1])
        hold off
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(strcat(path_pics, type, "EC_coverage_estVar_nsim",num2str(Msim), "n", num2str(nsubj),".png"), '-dpng')

end
