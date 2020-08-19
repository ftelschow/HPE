%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%  Simulation of EEC estimators and producing figures for the paper
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

%------ define parameters for this script
% simulation parameters
Nsim  = 1000;       % number of simulations
Nsubj = [10 50 100];% number of subjects/sample size

% Field parameters
TYPE = ["isotropic", "scale-space"];
L0   = 1;     % EC of domain
D    = 2;     % Dimension of domain (including scale)
T    = 50;    % size of domain

nu   = 5;         % parameter for isotropic
gamma = 4:.2:15;  % bandwiths for scale-space
sigma = 1;

% Thresholds
du = 0.05;
u = -6:du:6;

% Other parameters
FWERlevel = 0.05;

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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%---------- Simulation of EEC estimation using different estimators -------
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set seed for reproducibility
rng(23)
scale=0.9;

% loops over sample size and error fields
for type = TYPE
    for nsubj = Nsubj
        %------ generate data
        if strcmp(type, "isotropic")
            params = nu;
        elseif strcmp(type, "scale-space")
            params = gamma;
        end

        if nsubj ~=100
        tic
            [f, trueLKC] = generateField( nsubj, Nsim, [T T], type, params);
        toc
        else
            tic
            [f, ~] = generateField( nsubj/2, Nsim, [T T], type, params);
            toc
            tic
            [f1, trueLKC] = generateField( nsubj/2, Nsim, [T T], type, params);
            toc
            f = cat(3, f, f1);
            clear f1
            size(f)
        end


        %------ compute the empirical EC curves
        tic
        EC = EulerChar(f, u, D);
        toc
        % 1 minute per 100 fields for isotropic


        %------ estimate the LKCs and EEC with the HPE
        tic
        [LKC, EEC, rho] = HermiteEEC(EC, u, D, trueLKC, L0);
        toc
        Sigma = cov(LKC.hat1(:,:)');

        %------ evaluation: EEC
        % True accuracy
        C = rho(:, 1:D) * Sigma * rho(:, 1:D)';
        EEC_se = sqrt(diag(C) / nsubj);
        EEC_conf = cat(3, EEC.hatn - 1.96*repmat(EEC_se, 1, Nsim), EEC.hatn + 1.96*repmat(EEC_se, 1, Nsim));

        % Performance
        EEC_rmse = sqrt(mean(sum((EEC.hatn - repmat(EEC.true, 1, Nsim)).^2 .* repmat(exp(u'.^2/2), 1, Nsim)), 2))
        EEC_cover = mean((EEC_conf(:,:,1) < repmat(EEC.true, 1, Nsim)) & (repmat(EEC.true, 1, Nsim) < EEC_conf(:,:,2)), 2);
        EEC_cover_hat = mean((EEC.conf_hat(:,:,1) < repmat(EEC.true, 1, Nsim)) & (repmat(EEC.true, 1, Nsim) < EEC.conf_hat(:,:,2)), 2);


        %------ Nonparametric method
        EC_bar = squeeze(mean(EC, 2));
        EC_bar_se_hat = zeros(length(u), Nsim);
        for j = 1:Nsim
            EC_bar_se_hat(:,j) = sqrt(diag(cov(EC(:,:,j)')) / nsubj);
        end
        EC_bar_conf_hat = cat(3, EC_bar - 1.96*EC_bar_se_hat, EC_bar + 1.96*EC_bar_se_hat);

        % True accuracy
        EC_cov = cov(EC(:,:)');
        EC_bar_se = sqrt(diag(EC_cov) / nsubj);
        EC_bar_conf = cat(3, EC_bar - 1.96*repmat(EC_bar_se, 1, Nsim), EC_bar + 1.96*repmat(EC_bar_se, 1, Nsim));

        % Performance
        EC_bar_rmse = sqrt(mean(sum((EC_bar - repmat(EEC.true, 1, Nsim)).^2 .* repmat(exp(u'.^2/2), 1, Nsim)), 2))
        EC_bar_cover = mean((EC_bar_conf(:,:,1) < repmat(EEC.true, 1, Nsim)) & (repmat(EEC.true, 1, Nsim) < EC_bar_conf(:,:,2)), 2);
        EC_bar_cover_hat = mean((EC_bar_conf_hat(:,:,1) < repmat(EEC.true, 1, Nsim)) & (repmat(EEC.true, 1, Nsim) < EC_bar_conf_hat(:,:,2)), 2);


        %%------ Plot figures for the manuscript
        % Plot an instant of the random field
        figure(1), clf
        WidthFig = scale*500;
        HeightFig = scale*450;
            set(groot, 'defaultAxesTickLabelInterpreter','latex');
            set(groot, 'defaultLegendInterpreter','latex');
            set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
            set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig])
            i = 1; j = 1;
            imagesc(f(:,:,i,j)), axis equal tight
            colorbar;
            h = title('Single simulated field');
            set(h, 'Interpreter', 'latex');
            set(gca, 'FontSize', 20)
            set(gcf,'papersize',[12 12])
            xticks([10 20 30 40])
            xticklabels( {'10' '20' '30' '40'} )
            yticks([0 10 20 30 40])
            yticklabels( {'50' '40' '30' '20' '10'} )
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(strcat(path_pics, type, ".png"), '-dpng')

        % Plot single example of smoothing
        figure(2), clf, hold on,
        WidthFig = scale*500;
        HeightFig = scale*450;
        set(groot, 'defaultAxesTickLabelInterpreter','latex');
        set(groot, 'defaultLegendInterpreter','latex');
        set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
        set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig])
            h = title('Smoothing of single EC curve');
            set(h, 'Interpreter', 'latex');
            i = 1; j = 1;
            plot(u, EC(:,i,j), 'color', greyScale*[1 1 1])
            plot(u, EEC.hat1(:,i,j), 'Color', Vibrant(1,:), 'LineWidth', 2)
            h = xlabel('u');
            set(h, 'Interpreter', 'latex');
            h = ylabel('EC');
            set(h, 'Interpreter', 'latex');
            axis tight
            hold off
            xticks([-5 -2.5 0 2.5 5])
            xticklabels( {'-5', '-2.5', '0', '2.5', '5'} )
        set(gca, 'FontSize', 20)
        set(gcf,'papersize',[12 12])
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(strcat(path_pics, type, "EC1_nu5.png"), '-dpng')

        % Plot 10 examples of raw EC curves
        figure(3), clf, hold on,
        WidthFig = scale*500;
        HeightFig = scale*450;
        set(groot, 'defaultAxesTickLabelInterpreter','latex');
        set(groot, 'defaultLegendInterpreter','latex');
        set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
        set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig])
            h = title('EC curves');
            set(h, 'Interpreter', 'latex');
            j = 1;
            plot(u, EC(:,:,j), 'color', greyScale*[1 1 1])
            plot(u, EC_bar(:,j), 'Color', Vibrant(1,:), 'LineWidth', 2)
            plot(u, EC_bar_conf(:,j,1), 'Color', Vibrant(1,:), 'LineStyle', '--', 'LineWidth', 2)
            plot(u, EC_bar_conf(:,j,2), 'Color', Vibrant(1,:), 'LineStyle', '--', 'LineWidth', 2)
            plot(u, EEC.true, 'Color', Vibrant(5,:), 'LineWidth', 2)
            h = xlabel('u');
            set(h, 'Interpreter', 'latex');
            h = ylabel('EC');
            set(h, 'Interpreter', 'latex');
            axis tight
            hold off
            xticks([-5 -2.5 0 2.5 5])
            xticklabels( {'-5', '-2.5', '0', '2.5', '5'} )
        set(gca, 'FontSize', 20)
        set(gcf,'papersize',[12 12])
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(strcat(path_pics, type, "_average_nu5_n",num2str(nsubj),".png"), '-dpng')


        % Plot 10 examples of smoothed EC curves
        figure(4), clf, hold on,
        WidthFig = scale*500;
        HeightFig = scale*450;
        set(groot, 'defaultAxesTickLabelInterpreter','latex');
        set(groot, 'defaultLegendInterpreter','latex');
        set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
        set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig])
            h = title('Smoothed EC curves');
            set(h, 'Interpreter', 'latex');
            j = 1;
            plot(u, EEC.hat1(:,:,j), 'color', greyScale*[1 1 1])
            plot(u, EEC.hatn(:,j), 'Color', Vibrant(1,:), 'LineStyle', '-', 'LineWidth', 2)
            plot(u, EEC_conf(:,j,1), 'Color', Vibrant(1,:), 'LineStyle', '--', 'LineWidth', 2)
            plot(u, EEC_conf(:,j,2), 'Color', Vibrant(1,:), 'LineStyle', '--', 'LineWidth', 2)
            plot(u, EEC.true, 'Color', Vibrant(5,:), 'LineWidth', 2)
            h = xlabel('u');
            set(h, 'Interpreter', 'latex');
            h = ylabel('EC');
            set(h, 'Interpreter', 'latex');
            axis tight
            hold off
            xticks([-5 -2.5 0 2.5 5])
            xticklabels( {'-5', '-2.5', '0', '2.5', '5'} )
        set(gca, 'FontSize', 20)
        set(gcf,'papersize',[12 12])
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(strcat(path_pics, type, "_EEC_nu5_n",num2str(nsubj),".png"), '-dpng')

        % plot covariance of raw EC curves
        figure(5), clf
        WidthFig = scale*500;
        HeightFig = scale*450;
            set(groot, 'defaultAxesTickLabelInterpreter','latex');
            set(groot, 'defaultLegendInterpreter','latex');
            set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
            set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig])
            imagesc(u, u, flipud(EC_cov)), axis equal tight
            xticks([-5 -2.5 0 2.5 5])
            xticklabels( {'-5', '-2.5', '0', '2.5', '5'} )
            yticks([-5 -2.5 0 2.5 5])
            yticklabels( {'-5', '-2.5', '0', '2.5', '5'} )
            hax = gca;
            hax.YTickLabel = flipud(hax.YTickLabel);
            colorbar;
            h = title('Covariance of EC curves');
            set(h, 'Interpreter', 'latex');
            set(gca, 'FontSize', 20)
            set(gcf,'papersize',[12 12])
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(strcat(path_pics, type, "_EC-bar-cov",num2str(nsubj),".png"), '-dpng')

        % plot covariance of smoothed EC curves
        figure(6), clf
        WidthFig = scale*500;
        HeightFig = scale*450;
            set(groot, 'defaultAxesTickLabelInterpreter','latex');
            set(groot, 'defaultLegendInterpreter','latex');
            set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
            set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig])
            imagesc(u, u, flipud(C))
            axis equal tight
            xticks([-5 -2.5 0 2.5 5])
            xticklabels( {'-5', '-2.5', '0', '2.5', '5'} )
            yticks([-5 -2.5 0 2.5 5])
            yticklabels( {'-5', '-2.5', '0', '2.5', '5'} )
            hax = gca;
            hax.YTickLabel = flipud(hax.YTickLabel);
            colorbar;
            h = title('Covariance of HPE');
            set(h, 'Interpreter', 'latex');
            set(gca, 'FontSize', 20)
            set(gcf,'papersize',[12 12])
            axis square tight
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(strcat(path_pics, type, "_EEC-hat-cov",num2str(nsubj),".png"), '-dpng')

        % Comparison of variances
        figure(7), clf, hold on
        WidthFig = scale*600;
        HeightFig = scale*450;
            set(groot, 'defaultAxesTickLabelInterpreter','latex');
            set(groot, 'defaultLegendInterpreter','latex');
            set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
            set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig])
            h = title('Standard deviation');
            set(h, 'Interpreter', 'latex');
            plot(u, sqrt(diag(EC_cov)), 'Color', Vibrant(1,:), 'LineWidth', 2)
            plot(u, sqrt(diag(C)), 'Color', Vibrant(5,:), 'LineWidth', 2)
            h = xlabel('u');
            xticks([-5 -2.5 0 2.5 5])
            xticklabels( {'-5', '-2.5', '0', '2.5', '5'} )
            set(h, 'Interpreter', 'latex');
            set(gca, 'FontSize', 20)
            h = legend("$\bar\chi(u)$", "$\widehat{{\rm EEC}}(u)$", 'Location', 'northwest');
            legend('boxoff')
            set(h, 'Interpreter', 'latex');
            set(h, 'FontSize', 25)
            set(gcf,'papersize',[12 12])
            hold off
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(strcat(path_pics, type, "var",num2str(nsubj),".png"), '-dpng')


        % Coverage performance true variance
        figure(8), clf, hold on
        WidthFig = 400;
        HeightFig = 350;
            set(groot, 'defaultAxesTickLabelInterpreter','latex');
            set(groot, 'defaultLegendInterpreter','latex');
            set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
            set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig])
            plot(u, EC_bar_cover, 'Color', Vibrant(1,:), 'LineWidth', 2)
            plot(u, EEC_cover, 'Color', Vibrant(5,:), 'LineWidth', 2)
            plot([-8 8], [0.95 0.95], 'k--', 'LineWidth', 2)
            axis([u(1) u(end) 0 1])
            xlabel('u')
            xticks([-5 -2.5 0 2.5 5])
            xticklabels( {'-5', '-2.5', '0', '2.5', '5'} )
            yticks([0 0.2 0.4 0.6 0.8 1])
            yticklabels( {'0' '0.2' '0.4' '0.6' '0.8' '1'} )
            hold off
            h = xlabel('u');
            set(h, 'Interpreter', 'latex');
            set(gca, 'FontSize', 20)
            h = legend("$\bar\chi(u)$", "$\widehat{{\rm EEC}}(u)$", 'Location', 'south');
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
        print(strcat(path_pics, type, "EC_coverage_trueVar_nsim",num2str(Nsim), "n", num2str(nsubj),".png"), '-dpng')


        % Coverage performance estimated variance
        figure(9), clf, hold on
        WidthFig = 400;
        HeightFig = 350;
            set(groot, 'defaultAxesTickLabelInterpreter','latex');
            set(groot, 'defaultLegendInterpreter','latex');
            set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
            set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig])
            plot(u, EC_bar_cover_hat, 'Color', Vibrant(1,:), 'LineWidth', 2)
            plot(u, EEC_cover_hat, 'Color', Vibrant(5,:), 'LineWidth', 2)
            plot([-8 8], [0.95 0.95], 'k--', 'LineWidth', 2)
            axis([u(1) u(end) 0 1])
            xlabel('u')
            xticks([-5 -2.5 0 2.5 5])
            xticklabels( {'-5', '-2.5', '0', '2.5', '5'} )
            yticks([0 0.2 0.4 0.6 0.8 1])
            yticklabels( {'0' '0.2' '0.4' '0.6' '0.8' '1'} )
            hold off
            h = xlabel('u');
            set(h, 'Interpreter', 'latex');
            set(gca, 'FontSize', 20)
            h = legend("$\bar\chi(u)$", "$\widehat{{\rm EEC}}(u)$", 'Location', 'south');
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
        print(strcat(path_pics, type, "EC_coverage_estVar_nsim",num2str(Nsim), "n", num2str(nsubj),".png"), '-dpng')

    end
end