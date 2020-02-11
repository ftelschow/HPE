%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%    Plot the results from Script_Simulate_LKCestimates3D.m
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: This script plots the results from the simulation script
%              Script_Simulate_LKCestimates3D.m
%__________________________________________________________________________
% REFERENCES:
%
%__________________________________________________________________________
% AUTHOR: Fabian Telschow (ftelschow@ucsd.edu)
%         Wenyi Lin (wel316@ucsd.edu)
%__________________________________________________________________________
%------ prepare workspace
% clean history
clear
close all
% load generated profile containing the paths
load('paths.mat')

% create matlab paths to search for functions
addpath(genpath(path_main));

% %% Show results
% % standard color schemes 'https://personal.sron.nl/~pault/'
% BrightCol  = [[68 119 170];...    % blue
%               [102 204 238];...   % cyan
%               [34 136 51];...     % green
%               [204 187 68];...    % yellow
%               [238 102 119];...   % red
%               [170 51 119];...    % purple
%               [187 187 187]]/255; % grey
% 
% HighContr  = [[221, 170,  51];...   % yellow
%               [187,  85, 102];...   % red
%               [  0,  68, 136]]/255; % blue
% Vibrant    = [[0 119 187];... % blue
%               [51 187 238];...% cyan
%               [0 153 136];... % teal
%               [238 119 51];...% orange
%               [204 51 17];... % red
%               [238 51 119];...% magenta
%               [187 187 187]...% grey
%               ]/255;
% 
% colMat = Vibrant([1 3 4 5],:);
% 
% sfont = 20;
% addf = 2;
% 
% % load results
% load( strcat( "SimResults_Isotropic_T",...
%              num2str(T), "_D", num2str(D), "_Nfields",...
%              num2str(batch*Msim) ) )
% 
% 
% 
% tmp  = struct();
% tmp2 = LKC_HPE;
% tmp.hatn = tmp2;
% tmp.hatmean = mean( tmp.hatn, 4 );
% tmp.hatstd  = std( tmp.hatn, 0, 4 );
% LKC_HPE = tmp;
% 
% tmp  = struct();
% tmp2 = LKC_spm;
% tmp.hatn = tmp2;
% tmp.hatmean = mean( tmp.hatn, 4 );
% tmp.hatstd  = std( tmp.hatn, 0, 4 );
% LKC_spm = tmp;
% 
% scale = 0.9;
% WidthFig   = 1300;
% HeightFig  = WidthFig * scale;
% xvec       = (1:3) *15; %nvec;
% xtickcell  = {'50', '100', '200'};
% 
% set(groot, 'defaultAxesTickLabelInterpreter','latex');
% set(groot, 'defaultLegendInterpreter','latex');
% 
% % LKC plot
% for f = 1:length(FWHM)
% figure(f), clf, hold on;
% % Define size and location of the figure [xPos yPos WidthFig HeightFig]
% set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
% set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig])
% 
% 
% for i = 1:D
%         % generate subplot
%         subplot(1,D,i), hold on;
%         % Plot herm estimator
%         ErrorBar = errorbar(xvec-2, squeeze(LKC_HPE.hatmean(i,:,f)),...
%                                     squeeze(LKC_HPE.hatstd(i,:,f)), 'o');
%         set(ErrorBar,'Color', colMat(1,:),'LineWidth', 2)
%  
%         ErrorBar = errorbar(xvec+2, squeeze(LKC_spm.hatmean(i,:,f)),...
%                                     squeeze(LKC_spm.hatstd(i,:,f)), 'o');
%         set(ErrorBar,'Color', colMat(3,:),'LineWidth', 2)
% 
%         % Plot the true value
%         plot([xvec(1)-10 xvec(end)+10],[trueLKC(i,f) trueLKC(i,f)],'k', 'LineWidth', 2 ), hold on
% 
%         if i == 2
%         	% Add title
%             h = title( strcat('FWHM=',num2str(FWHM(f))), 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
%         end
% 
%         % Modify gloabal font size for this plot
%         set(gca,'FontSize', sfont)
% 
%         % Change axis style
%         xlim([xvec(1)-10 xvec(end)+10])
%         xticks(xvec)
%         xticklabels(xtickcell)
% 
%         h = xlabel('Sample Size [{\it N}\,]', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
%         h = ylabel( strcat('$\mathcal{L}_', num2str(i),'$'), 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
%         % Add legend
%         legend( 'HPE', 'spm12', 'Location', 'northeast' );
%  end
% %    saveas( gcf,strcat('pics/',outputname) )
%  set(gcf,'papersize',[12 12*scale])
%  fig = gcf;
%  fig.PaperPositionMode = 'auto';
%  fig_pos = fig.PaperPosition;
%  fig.PaperSize = [fig_pos(3) fig_pos(4)];
% hold off;
%  print( strcat( "Plots_LKC_FWHM_",num2str(FWHM(f)),"_SimResults_Isotropic_T",...
%              num2str(T), "_D", num2str(D), "_Nfields",...
%              num2str(batch*Msim) ), '-dpng')
% end
% 
% %% Show results
% % standard color schemes 'https://personal.sron.nl/~pault/'
% BrightCol  = [[68 119 170];...    % blue
%               [102 204 238];...   % cyan
%               [34 136 51];...     % green
%               [204 187 68];...    % yellow
%               [238 102 119];...   % red
%               [170 51 119];...    % purple
%               [187 187 187]]/255; % grey
% 
% HighContr  = [[221, 170,  51];...   % yellow
%               [187,  85, 102];...   % red
%               [  0,  68, 136]]/255; % blue
% Vibrant    = [[0 119 187];... % blue
%               [51 187 238];...% cyan
%               [0 153 136];... % teal
%               [238 119 51];...% orange
%               [204 51 17];... % red
%               [238 51 119];...% magenta
%               [187 187 187]...% grey
%               ]/255;
% 
% colMat = Vibrant([1 3 4 5],:);
% 
% sfont = 20;
% addf = 2;
% 
% % load results
% load( strcat( "SimResults_Isotropic_T",...
%              num2str(T), "_D", num2str(D), "_Nfields",...
%              num2str(batch*Msim) ) )
% 
% tmp  = struct();
% tmp2 = thr_HPE;
% tmp.hatn = squeeze(tmp2);
% tmp.hatmean = mean( tmp.hatn, 3 );
% tmp.hatstd  = std( tmp.hatn, 0, 3 );
% thr_HPE = tmp;
% 
% tmp  = struct();
% tmp2 = thr_spm;
% tmp.hatn = squeeze(tmp2);
% tmp.hatmean = mean( tmp.hatn, 3 );
% tmp.hatstd  = std( tmp.hatn, 0, 3 );
% thr_spm = tmp;
% 
% scale = 0.9;
% WidthFig   = 1300;
% HeightFig  = WidthFig * scale;
% xvec       = (1:3) *15; %nvec;
% xtickcell  = {'50', '100', '200'};
% 
% set(groot, 'defaultAxesTickLabelInterpreter','latex');
% set(groot, 'defaultLegendInterpreter','latex');
% figure(1), clf, hold on;
% % LKC plot
% for f = 1:length(FWHM)
% 
% % Define size and location of the figure [xPos yPos WidthFig HeightFig]
% set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
% set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig])
%         trueu   = get_EECthreshold(FWER, uvals, trueLKC(:,f));
%         % generate subplot
%         subplot(2,length(FWHM)/2,f), hold on;
%         % Plot herm estimator
%         ErrorBar = errorbar(xvec-2, squeeze(thr_HPE.hatmean(:,f)),...
%                                     squeeze(thr_HPE.hatstd(:,f)), 'o');
%         set( ErrorBar,'Color',colMat(1,:),'LineWidth', 2 )
%         ErrorBar = errorbar(xvec+2, squeeze(thr_spm.hatmean(:,f)),...
%                                     squeeze(thr_spm.hatstd(:,f)), 'o' );
%         set( ErrorBar,'Color',colMat(3,:),'LineWidth', 2)
% 
%         % Plot the true value
%         plot([xvec(1)-10 xvec(end)+10],[trueu trueu],'k', 'LineWidth', 2 ), hold on
%         % Add title
%         h = title( strcat('FWHM=',num2str(FWHM(f))), 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
% 
%         % Modify gloabal font size for this plot
%         set(gca,'FontSize', sfont)
% 
%         % Change axis style
%         xlim([xvec(1)-10 xvec(end)+10])
%         xticks(xvec)
%         xticklabels(xtickcell)
% 
%         h = xlabel( 'Sample Size [{\it N}\,]', 'fontsize', sfont+addf );  set(h, 'Interpreter', 'latex');
%         h = ylabel( strcat('FWER threshold'), 'fontsize', sfont+addf ); set(h, 'Interpreter', 'latex');
%         % Add legend
%         legend( 'HPE', 'spm12', 'Location', 'northeast' );
%  end
% %    saveas( gcf,strcat('pics/',outputname) )
%  set(gcf,'papersize',[12 12*scale])
%  fig = gcf;
%  fig.PaperPositionMode = 'auto';
%  fig_pos = fig.PaperPosition;
%  fig.PaperSize = [fig_pos(3) fig_pos(4)];
%  print( strcat( "Plots_FWERthresh_SimResults_Isotropic_T",...
%              num2str(T), "_D", num2str(D), "_Nfields",...
%              num2str(batch*Msim) ), '-dpng')
% hold off;