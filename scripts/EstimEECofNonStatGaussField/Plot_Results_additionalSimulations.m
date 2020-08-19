%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%  This script produces the simulations reported in
%%%%  Telschow, Fabian, et al.
%%%%  "Estimation of Expected Euler Characteristic Curves of
%%%%   Nonstationary Smooth Gaussian Random Fields."
%%%%  arXiv preprint arXiv:1908.02493 (2020).
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Fabian Telschow
%--------------------------------------------------------------------------
%------ prepare workspace
% clean history
clear
close all
% load generated profile containing the paths
load('paths.mat')

% create matlab paths to search for functions
addpath(genpath('/home/drtea/matlabToolboxes/HPE'));

% standard color schemes 'https://personal.sron.nl/~pault/'
BrightCol  = [[68 119 170];...    % blue
              [102 204 238];...   % cyan
              [34 136 51];...     % green
              [204 187 68];...    % yellow
              [238 102 119];...   % red
              [170 51 119];...    % purple
              [187 187 187]]/255; % grey

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

colMat = Vibrant([5 1 2 4],:);

% Global figure settings
sfont = 20;
addf  = 5;
scale = 3.5/12;

nsim    = 1e3;
nu      = 5;
nvec    = [10 50 100 200];
T       = 50;
D       = 2;

outputname = strcat( 'isotropic_D', num2str(D), 'T', num2str(T),...
                    '_params', num2str(nu) );

sim_identifier = strcat( 'maxN', num2str(max(nvec)), outputname );

% names for output graphics
outname = [ ...
            strcat(sim_identifier,"_LKCestim_Theory",".png");...
            strcat(sim_identifier,"_LKCestim_DemeanVar1",".png");...
            strcat(sim_identifier,"_LKCestim_Var1",".png");...
            strcat(sim_identifier,"_LKCestim_Demean",".png");...
            ];

WidthFig   = 1300;
HeightFig  = WidthFig * scale;
xvec       = (1:4) *15; %nvec;
xtickcell  = {'10', '50', '100', '200'};
yvec1      = [12 13 14 15 16 17 18 19 20 21 22];
yvec2      = [40 45 50 55 60 65 70 75 80 85];
ytickcell1 = {'12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22'};
ytickcell2 = {'40' '45' '50' '55' '60', '65', '70', '75' '80' '85'};


yvec12      = [12 13 14 15 16];
yvec22      = [40 45 50];
ytickcell12 = {'12' '13' '14' '15' '16'};
ytickcell22 = {'40' '45' '50'};

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');


%% Dependence on multipliers in the bootstrap

% files to be loaded
load( strcat( path_results, '/simLKChermBm_', sim_identifier,'.mat') )
% Global figure settings
trueLKC    = L(1:2);

for l = 1:2
%%%
    outputname = outname(l);
    figure(l), clf, hold on;
    % Define size and location of the figure [xPos yPos WidthFig HeightFig]
    set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
    set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig])

    % LKC plot
    for i = [1 2]
            % generate subplot
            subplot(1,2,i), hold on;
            % Plot gaussian bootstrap estimator
            ErrorBar = errorbar(xvec-3, LKChermB.hatmean(i,:,l), LKChermB.hatstd(i,:,l), 'o'); hold on;
            set(ErrorBar,'Color',colMat(1,:),'LineWidth',2)
            % Plot multinomial estimator
            ErrorBar = errorbar(xvec-1, LKChermBn.hatmean(i,:,l), LKChermBn.hatstd(i,:,l), 'o'); hold on ;
            set(ErrorBar,'Color',colMat(2,:),'LineWidth',2)
            % Plot corrected multinomial bootstrap estimator
            ErrorBar = errorbar(xvec+1, LKChermBm.hatmean(i,:,l) ,LKChermBm.hatstd(i,:,l), 'o'); hold on ;
            set(ErrorBar,'Color',colMat(3,:),'LineWidth',2)
            % Plot rademacher estimator
            ErrorBar = errorbar(xvec+3, LKChermBr.hatmean(i,:,l), LKChermBr.hatstd(i,:,l), 'o'); hold on ;
            set(ErrorBar,'Color',colMat(4,:),'LineWidth',2)
 
            % Plot the true value
            plot( [ xvec(1) - 10 xvec(end) + 10 ], [ trueLKC(i) trueLKC(i) ],...
                  'k', 'LineWidth', 2 ), hold on


            % Modify gloabal font size for this plot
            set( gca,'FontSize', sfont )

            % Change axis style
            xlim([xvec(1)-10 xvec(end)+10])
            xticks(xvec)
            xticklabels(xtickcell)

            h = xlabel('Sample Size [{\it N}]', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
            if i==1
                if l==1
                    ylim([yvec1(1) yvec1(end)+0.5]);
                    yticks(yvec1);
                    yticklabels(ytickcell1);
                else
                    ylim([yvec12(1) yvec12(end)+0.5]);
                    yticks(yvec12);
                    yticklabels(ytickcell12);
                end
                h = ylabel('$\mathcal{L}_1$', 'fontsize', sfont+addf);
                set(h, 'Interpreter', 'latex');
            else
                if l==1
                    ylim([yvec2(1)-3 yvec2(end)+3]);
                    yticks(yvec2);
                    yticklabels(ytickcell2);
                else
                    ylim([yvec22(1)-3 yvec22(end)+3]);
                    yticks(yvec22);
                    yticklabels(ytickcell22);
                end
                h = ylabel('$\mathcal{L}_2$', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
            end

            % Add legend
            if i==1 && l==2
               legend( 'bHPE (Gaussian)', 'bHPE (Multinomial)', 'bHPE (mod. Multinomial)', 'bHPE (Rademacher)',...
                            'Location', 'northeast' );
                set(legend, 'fontsize', sfont);
                legend boxoff
                %set(legend,'color','none')
            end
     end
%    saveas( gcf,strcat('pics/',outputname) )
     set(gcf,'papersize',[12 12*scale])
     fig = gcf;
     fig.PaperPositionMode = 'auto';
     fig_pos = fig.PaperPosition;
     fig.PaperSize = [fig_pos(3) fig_pos(4)];
     print(strcat(path_pics,outputname,'_additionalSimulations' ), '-dpng')
    hold off;
end


%% Dependence on connectivity
% y axis
yvec1      = [12 13 14 15 16];
yvec2      = [ 40 45 50 55];
ytickcell1 = {'12' '13' '14' '15' '16'};
ytickcell2 = {'40' '45' '50' '55'};

% files to be loaded
load( strcat( path_results, 'sim_varyingCC_LKCherm_', sim_identifier,'.mat') )
% Global figure settings
trueLKC    = L(1:2);

for l = 1
%%%
    outputname = outname(l);
    figure(l), clf, hold on;
    % Define size and location of the figure [xPos yPos WidthFig HeightFig]
    set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
    set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig])

    % LKC plot
    for i = [1 2]
            % generate subplot
            subplot(1,2,i), hold on;
            % Plot multinomial estimator
            ErrorBar = errorbar(xvec-1, LKCherm4.hatmean(i,:,l), LKCherm4.hatstd(i,:,l), 'o'); hold on ;
            set(ErrorBar,'Color',colMat(1,:),'LineWidth',2)
            % Plot corrected multinomial bootstrap estimator
            ErrorBar = errorbar(xvec+1, LKCherm8.hatmean(i,:,l) ,LKCherm8.hatstd(i,:,l), 'o'); hold on ;
            set(ErrorBar,'Color',colMat(3,:),'LineWidth',2)
 
            % Plot the true value
            plot( [ xvec(1) - 10 xvec(end) + 10 ], [ trueLKC(i) trueLKC(i) ],...
                  'k', 'LineWidth', 2 ), hold on


            % Modify gloabal font size for this plot
            set( gca,'FontSize', sfont )

            % Change axis style
            xlim([xvec(1)-10 xvec(end)+10])
            xticks(xvec)
            xticklabels(xtickcell)

            h = xlabel('Sample Size [{\it N}]', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
            if i==1
                ylim([yvec1(1) yvec1(end)+0.5]);
                yticks(yvec1);
                yticklabels(ytickcell1);
                h = ylabel('$\mathcal{L}_1$', 'fontsize', sfont+addf);
                set(h, 'Interpreter', 'latex');
            else
                ylim([yvec2(1)-3 yvec2(end)+3]);
                yticks(yvec2);
                yticklabels(ytickcell2);
                h = ylabel('$\mathcal{L}_2$', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
            end

            % Add legend
            if i==1 && l==1
               legend( 'HPE (4 cc)', 'HPE (8 cc)',...
                            'Location', 'northeast' );
                set(legend, 'fontsize', sfont);
                legend boxoff
                %set(legend,'color','none')
            end
     end
%    saveas( gcf,strcat('pics/',outputname) )
     set(gcf,'papersize',[12 12*scale])
     fig = gcf;
     fig.PaperPositionMode = 'auto';
     fig_pos = fig.PaperPosition;
     fig.PaperSize = [fig_pos(3) fig_pos(4)];
     print(strcat(path_pics,outputname,'_additionalSimulations_CCdependence' ), '-dpng')
    hold off;
end

%% Plot the dependence on resolution
% Change the axis labels
xvec       = (1:5) *15; %nvec;
xtickcell  = {'0.167', '0.25', '0.33', '0.5', '1'};
yvec1      = [12 13 14 15];
yvec2      = [ 40 45 50 55];
yvec2      = [ 42.5 45 47.5 50 ];
ytickcell1 = { '12' '13' '14' '15' };
ytickcell2 = { '42.5' '45' '47.5' '50' };

% files to be loaded
load( strcat( path_results, 'sim_varyingResolution_LKCherm.mat'))
% Global figure settings
trueLKC    = L(1:2);

figure, clf, hold on;
% Define size and location of the figure [xPos yPos WidthFig HeightFig]
set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig])

% LKC plot
for i = [1 2]
        % generate subplot
        subplot(1,2,i), hold on;
        % Plot multinomial estimator
        ErrorBar = errorbar(xvec-1, LKCherm4.hatmean(i,:), LKCherm4.hatstd(i,:), 'o'); hold on ;
        set(ErrorBar,'Color',colMat(1,:),'LineWidth',2)
        % Plot corrected multinomial bootstrap estimator
        ErrorBar = errorbar(xvec+1, LKCherm8.hatmean(i,:) ,LKCherm8.hatstd(i,:), 'o'); hold on ;
        set(ErrorBar,'Color',colMat(3,:),'LineWidth',2)

        % Plot the true value
        plot( [ xvec(1) - 10 xvec(end) + 10 ], [ trueLKC(i) trueLKC(i) ],...
              'k', 'LineWidth', 2 ), hold on


        % Modify gloabal font size for this plot
        set( gca,'FontSize', sfont )

        % Change axis style
        xlim([xvec(1)-10 xvec(end)+10])
        xticks(xvec)
        xticklabels(xtickcell)

        h = xlabel('Resolution', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
        if i==1
            ylim([yvec1(1) yvec1(end)+0.5]);
            yticks(yvec1);
            yticklabels(ytickcell1);
            h = ylabel('$\mathcal{L}_1$', 'fontsize', sfont+addf);
            set(h, 'Interpreter', 'latex');
        else
            ylim([yvec2(1)-0.5 yvec2(end)+0.5]);
            yticks(yvec2);
            yticklabels(ytickcell2);
            h = ylabel('$\mathcal{L}_2$', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
        end

        % Add legend
        if i==1
           legend( 'HPE (4 cc)', 'HPE (8 cc)',...
                        'Location', 'northeast' );
            set(legend, 'fontsize', sfont);
            legend boxoff
            %set(legend,'color','none')
        end
end

set(gcf,'papersize',[12 12*scale])
 fig = gcf;
 fig.PaperPositionMode = 'auto';
 fig_pos = fig.PaperPosition;
 fig.PaperSize = [fig_pos(3) fig_pos(4)];
 print(strcat(path_pics,'varyingResolution_LKCherm' ), '-dpng')
hold off;