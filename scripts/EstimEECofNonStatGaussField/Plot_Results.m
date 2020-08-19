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

colMat = Vibrant([1 3 4 5],:);

% Global figure settings
sfont = 20;
addf  = 5;
scale = 3.5/12;

%% visualize Simulation Results - Istropic
nsim    = 1e3;
nu      = 5;
nvec    = [10 20 50 75 100 150 200];
T       = 50;
D       = 2;

outputname = strcat( 'isotropic_D', num2str(D), 'T', num2str(T),...
                    '_params', num2str(nu));

sim_identifier = strcat( 'maxN', num2str(max(nvec)), outputname );

% files to be loaded
load( strcat( path_results, '/simLKCherm_', sim_identifier,'.mat') )
load( strcat( path_results, '/simLKCwarp','_',...
             sim_identifier,'.mat') )
load( strcat( path_results, '/simLKChermB_', sim_identifier,'.mat') )
load( strcat( path_results, '/simLKCiso_', sim_identifier,'.mat') )

% names for output graphics
outname = [ ...
            strcat(sim_identifier,"_LKCestim_Theory",".png");...
            strcat(sim_identifier,"_LKCestim_DemeanVar1",".png");...
            strcat(sim_identifier,"_LKCestim_Var1",".png");...
            strcat(sim_identifier,"_LKCestim_Demean",".png");...
            ];

% Global figure settings
trueLKC    = L(1:2);
trueLKC2   = [ L(1) 47.7 ]

WidthFig   = 1300;
HeightFig  = WidthFig * scale;
xvec       = (1:7) *15; %nvec;
xtickcell  = {'10', '30', '50', '75', '100', '150', '200'};
yvec1      = [12 13 14 15 16];
yvec2      = [40 45 50 55];
ytickcell1 = {'12' '13' '14' '15' '16'};
ytickcell2 = {'40' '45' '50' '55'};

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');


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
            % Plot herm estimator
            ErrorBar = errorbar(xvec-4, LKCherm.hatmean(i,:,l), LKCherm.hatstd(i,:,l), 'o'); hold on;
            set(ErrorBar,'Color',colMat(1,:),'LineWidth',2)
            % Plot herm bootstrap estimator
            ErrorBar = errorbar(xvec-2, LKChermB.hatmean(i,:,l) ,LKChermB.hatstd(i,:,l), 'o'); hold on ;
            set(ErrorBar,'Color',colMat(2,:),'LineWidth',2)
            % Plot warping estimator
            ErrorBar = errorbar(xvec, LKCwarp.hatmean(i,:,l), LKCwarp.hatstd(i,:,l), 'o'); hold on ;
            set(ErrorBar,'Color',colMat(3,:),'LineWidth',2)
            % Plot Worsleys isotropic estimator
            ErrorBar = errorbar(xvec+2, LKCiso.hatmean(i,:,l) ,LKCiso.hatstd(i,:,l), 'o'); hold on ;
            set(ErrorBar,'Color',colMat(4,:),'LineWidth',2)

            % Plot the true value
            plot([xvec(1)-10 xvec(end)+10],[trueLKC(i) trueLKC(i)],'k', 'LineWidth', 2 ), hold on
            plot([xvec(1)-10 xvec(end)+10],[trueLKC2(i) trueLKC2(i)],'k--', 'LineWidth', 2 ), hold on

            % Modify gloabal font size for this plot
            set(gca,'FontSize', sfont)

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
               legend( 'HPE', 'bHPE', 'WarpE', 'IsotE',...
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
     print(strcat(path_pics,outputname), '-dpng')
    hold off;
end

%% visualize Simulation Results - Scale
nsim    = 1e3;
nvec    = [10 20 50 75 100 150 200];
T       = 50;
D       = 2;
gamma   = [4 15];

outputname = strcat( 'scale-space_D', num2str(D), 'T', num2str(T),...
                    '_params', num2str(gamma(1)),...
                    '_', num2str(gamma(end)) );

sim_identifier = strcat( 'maxN', num2str(max(nvec)), outputname );

% files to be loaded
load( strcat( path_results, '/simLKCherm_', sim_identifier,'.mat') )
load( strcat( path_results, '/simLKCwarp','_',...
             sim_identifier,'.mat') )
load( strcat( path_results, '/simLKChermB_', sim_identifier,'.mat') )
load( strcat( path_results, '/simLKCiso_', sim_identifier,'.mat') )

% names for output graphics
outname = [ ...
            strcat(sim_identifier,"_LKCestim_Theory",".png");...
            strcat(sim_identifier,"_LKCestim_DemeanVar1",".png");...
            strcat(sim_identifier,"_LKCestim_Var1",".png");...
            strcat(sim_identifier,"_LKCestim_Demean",".png");...
            ];

% Global figure settings
trueLKC    = L(1:2);
WidthFig   = 1300;
HeightFig  = WidthFig * scale;
xvec       = (1:7)*15;% nvec;
xtickcell  = {'10', '30', '50', '75', '100', '150', '200'};
yvec1      = [5.5 6 6.5 7 7.5, 8];
ytickcell1 = {'5.5' '6' '6.5' '7' '7.5', '8'};
yvec2      = [1 2 3 4 5 6 7 8];
ytickcell2 = {'1' '2' '3' '4' '5' '6' '7' '8'};

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');


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
            % Plot herm estimator
            ErrorBar = errorbar(xvec-4, LKCherm.hatmean(i,:,l), LKCherm.hatstd(i,:,l), 'o'); hold on;
            set(ErrorBar,'Color',colMat(1,:),'LineWidth',2)
            % Plot herm bootstrap estimator
            ErrorBar = errorbar(xvec-2, LKChermB.hatmean(i,:,l) ,LKChermB.hatstd(i,:,l), 'o'); hold on ;
            set(ErrorBar,'Color',colMat(2,:),'LineWidth',2)
            % Plot warping estimator
            ErrorBar = errorbar(xvec, LKCwarp.hatmean(i,:,l), LKCwarp.hatstd(i,:,l), 'o'); hold on ;
            set(ErrorBar,'Color',colMat(3,:),'LineWidth',2)
            % Plot Worsleys isotropic estimator
                ErrorBar = errorbar(xvec+2, LKCiso.hatmean(i,:,l) ,LKCiso.hatstd(i,:,l), 'o'); hold on ;
                set(ErrorBar,'Color',colMat(4,:),'LineWidth',2)
            
            % Plot the true value
            plot([xvec(1)-10 xvec(end)+10],[trueLKC(i) trueLKC(i)],'k', 'LineWidth', 2 ), hold on


            % Modify gloabal font size for this plot
            set(gca,'FontSize', sfont)

            % Change axis style
            xlim([xvec(1)-10 xvec(end)+10])
            xticks(xvec)
            xticklabels(xtickcell)

            h = xlabel('Sample Size [{\it N}]', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
            if i==1
                ylim([yvec1(1) yvec1(end)]);
                yticks(yvec1);
                yticklabels(ytickcell1);
                h = ylabel('$\mathcal{L}_1$', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
            else
                ylim([yvec2(1) yvec2(end)]);
                yticks(yvec2);
                yticklabels(ytickcell2);
                h = ylabel('$\mathcal{L}_2$', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
            end
            
            if(i==1)
                if l==1
                    legend( 'HPE', 'bHPE', 'WarpE','IsotE',...
                            'Location', 'northeast' );
                    set(legend, 'fontsize', sfont);
                                    legend boxoff
                    set(legend,'color','none')
                end
                X = [0.375 0.375];
                Y = [0.33   0.22];
                h = annotation('textarrow',X,Y,'String','IsotE', 'Color', colMat(4,:), 'LineWidth', 2, 'fontsize', sfont-5);
                set(h, 'Interpreter', 'latex');
            end
     end
     set(gcf,'papersize',[12 12*scale])
     fig = gcf;
     fig.PaperPositionMode = 'auto';
     fig_pos = fig.PaperPosition;
     fig.PaperSize = [fig_pos(3) fig_pos(4)];
     print(strcat(path_pics,outputname), '-dpng')
    hold off;
end

%% visualize Simulation Results - Non gauss
nu      = 5;
nvec    = [10 20 50 75 100 150 200];
T       = 50;
D       = 2;

outputname = strcat( 'nongauss_D', num2str(D), 'T', num2str(T),...
                    '_params', num2str(nu));

sim_identifier = strcat( 'maxN', num2str(max(nvec)), outputname );

% files to be loaded
load( strcat( path_results, '/simLKCherm_', sim_identifier,'.mat') )
load( strcat( path_results, '/simLKCwarp','_',...
             sim_identifier,'.mat') )
load( strcat( path_results, '/simLKChermB_', sim_identifier,'.mat') )
load( strcat( path_results, '/simLKCiso_', sim_identifier,'.mat') )

% names for output graphics
outname = [ ...
            strcat(sim_identifier,"_LKCestim_Theory",".png");...
            strcat(sim_identifier,"_LKCestim_DemeanVar1",".png");...
            strcat(sim_identifier,"_LKCestim_Var1",".png");...
            strcat(sim_identifier,"_LKCestim_Demean",".png");...
            ];

% Global figure settings
trueLKC    = L(1:2);
WidthFig   = 1300;
HeightFig  = WidthFig * scale;
xvec       = (1:7)*15;
xtickcell  = {'10', '30', '50', '75', '100', '150', '200'};
yvec1      = [12 13 14 15 16];
yvec2      = [35 40 45 50 55];
ytickcell1 = {'12' '13' '14' '15' '16'};
ytickcell2 = {'35' '40' '45' '50' '55' '57.5'};

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');


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
            % Plot herm estimator
            ErrorBar = errorbar(xvec-4, LKCherm.hatmean(i,:,l), LKCherm.hatstd(i,:,l), 'o'); hold on;
            set(ErrorBar,'Color',colMat(1,:),'LineWidth',2)
            % Plot herm bootstrap estimator
            ErrorBar = errorbar(xvec-2, LKChermB.hatmean(i,:,l) ,LKChermB.hatstd(i,:,l), 'o'); hold on ;
            set(ErrorBar,'Color',colMat(2,:),'LineWidth',2)
            % Plot warping estimator
            ErrorBar = errorbar(xvec, LKCwarp.hatmean(i,:,l), LKCwarp.hatstd(i,:,l), 'o'); hold on ;
            set(ErrorBar,'Color',colMat(3,:),'LineWidth',2)
            % Plot Worsleys isotropic estimator
            ErrorBar = errorbar(xvec+2, LKCiso.hatmean(i,:,l) ,LKCiso.hatstd(i,:,l), 'o'); hold on ;
            set(ErrorBar,'Color',colMat(4,:),'LineWidth',2)

            % Plot the true value
            plot([xvec(1)-10 xvec(end)+10],[trueLKC(i) trueLKC(i)],'k', 'LineWidth', 2 ), hold on

            % Modify global font size for this plot
            set(gca,'FontSize', sfont)

            % Change axis style
            xlim([xvec(1)-10 xvec(end)+10])
            xticks(xvec)
            xticklabels(xtickcell)

            h = xlabel('Sample Size [{\it N}]', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
            if i==1
                ylim([yvec1(1) yvec1(end)]);
                yticks(yvec1);
                yticklabels(ytickcell1);
                h = ylabel('$\mathcal{L}_1$', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
            else
                ylim([yvec2(1)-1 yvec2(end)+2]);
                yticks(yvec2);
                yticklabels(ytickcell2);
                h = ylabel('$\mathcal{L}_2$', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
            end

            % Add legend
            if i==2 && l==1
               legend( 'HPE', 'bHPE', 'WarpE', 'IsotE',...
                            'Location', 'southeast' );
                set(legend, 'fontsize', sfont);
                legend boxoff
            end
     end
     set(gcf,'papersize',[12 12*scale])
     fig = gcf;
     fig.PaperPositionMode = 'auto';
     fig_pos = fig.PaperPosition;
     fig.PaperSize = [fig_pos(3) fig_pos(4)];
     print(strcat(path_pics,outputname), '-dpng')
    hold off;
end


%% visualize Simulation Results - Dependence of Smoothness
nvec    = [10 20 50 75 100];
T       = 50;
D       = 2;
nuVec   = [3 4 5 6 7]; % parameter for isotropic

outputname = strcat( 'D', num2str(D),'T', num2str(T),...
                     'Nsim', num2str(nsim),'_params', num2str(nuVec(1)),'_',num2str(nuVec(end)) );

sim_identifier = strcat('simSmoothnessLKC_maxN', num2str(max(nvec)), outputname)

% files to be loaded
load( strcat( 'simulations/', sim_identifier,'.mat') )
load( strcat( 'simulations/simLKCwarp_',...
             sim_identifier,'.mat') )
load( strcat( 'simulations/simLKChermB_', sim_identifier,'Weightsgaussian.mat') )
load( strcat( 'simulations/simLKCspm_', sim_identifier,'.mat') )

if simpleEst
    load( strcat( 'simulations/simLKCisoCorr_', sim_identifier,'.mat') )
end

% names for output graphics
outname = [ ...
            strcat(sim_identifier,"_LKCestim_TheoryBiasCor",num2str(biasCor),".png");...
            strcat(sim_identifier,"_LKCestim_DemeanVar1BiasCor",num2str(biasCor),".png");...
            strcat(sim_identifier,"_LKCestim_Var1BiasCor",num2str(biasCor),".png");...
            strcat(sim_identifier,"_LKCestim_DemeanBiasCor",num2str(biasCor),".png");...
            ];

% Global figure settings
trueLKC    = L(1:2);
scale      = 8 / 13;
WidthFig   = 1300;
HeightFig  = WidthFig * scale;
sfont      = 22;
xvec       = nvec;
xtickcell  = {'10', '30', '50', '75', '100'};
yvec1      = [12 12.5 13 13.5 14 14.5 15 15.5 16];
yvec2      = [40 42.5 45 47.5 50 52.5 55];
ytickcell1 = {'12' '12.5' '13' '13.5' '14' '14.5' '15' '15.5' '16'};
ytickcell2 = {'40' '42.5' '45' '47.5' '50' '52.5' '55'};

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');


for l = [1]
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
            if l==1 ||l==2
                o = 2;
            else
                o = 0;
            end
            ErrorBar = errorbar(nvec-2-o, LKCherm.hatmean(i,:,l), LKCherm.hatstd(i,:,l), 'o'); hold on;
            set(ErrorBar,'Color','red','LineWidth',2)
            if l==1 ||l==2
                ErrorBar = errorbar(nvec-2, LKChermB.hatmean(i,:,l) ,LKChermB.hatstd(i,:,l), 'o'); hold on ;
                set(ErrorBar,'Color','magenta','LineWidth',2)                
            end
            ErrorBar = errorbar(nvec, LKCwarp.hatmean(i,:,l), LKCwarp.hatstd(i,:,l), 'o'); hold on ;
            set(ErrorBar,'Color','blue','LineWidth',2)
            ErrorBar = errorbar(nvec+2, LKCspm.hatmean(i,:,l) ,LKCspm.hatstd(i,:,l), 'o'); hold on ;
            set(ErrorBar,'Color','cyan','LineWidth',2)
            if simpleEst
                ErrorBar = errorbar(nvec+4, LKCiso.hatmean(i,:,l) ,LKCiso.hatstd(i,:,l), 'o'); hold on ;
                set(ErrorBar,'Color','yellow','LineWidth',2)
            end

            % Plot the true value
            plot([nvec(1)-10 nvec(end)+10],[trueLKC(i) trueLKC(i)],'k', 'LineWidth', 2 ), hold on


            % Modify gloabal font size for this plot
            set(gca,'FontSize', sfont)

            % Change axis style
            xlim([nvec(1)-10 nvec(end)+10])
            xticks(xvec)
            xticklabels(xtickcell)

            h = xlabel('Sample Size [{\it N}]', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
            if i==1
                titlename = 'LKC1';
                ylim([yvec1(1) yvec1(end)]);
                yticks(yvec1);
                yticklabels(ytickcell1);
                h = ylabel('$\mathcal{L}_1$', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
            else
                titlename = 'LKC2';
                ylim([yvec2(1) yvec2(end)]);
                yticks(yvec2);
                yticklabels(ytickcell2);
                h = ylabel('$\mathcal{L}_2$', 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
            end

            % Add title
            h = title(titlename, 'fontsize', sfont+10); set(h, 'Interpreter', 'latex');
            hold on

            % Add legend
            if i==1
                if simpleEst
                    if l==1 ||l==2
                        legend( 'Herm Proj', 'Herm gBoot Proj', 'WarpE', 'SPM12', 'isotropic',...
                                'Location', 'northeast' );
                    else
                        legend( 'Hermite Proj', 'WarpE', 'SPM12', 'isotropic',...
                                'Location', 'northeast' );
                    end
                else
                    if l==1 ||l==2
                        legend( 'Herm Proj', 'Herm gBoot Proj', 'WarpE', 'SPM12',...
                                'Location', 'northeast' );
                    else
                        legend( 'Hermite Proj', 'WarpE', 'SPM12',...
                                'Location', 'northeast' );
                    end
                end
                set(legend, 'fontsize', sfont);
                set(legend,'color','none')
            end
    end
%    saveas( gcf,strcat('pics/',outputname) )
    set(gcf,'papersize',[12 12*scale])
     print(strcat('pics/',outputname), '-dpng', '-fillpage')
    hold off;
end