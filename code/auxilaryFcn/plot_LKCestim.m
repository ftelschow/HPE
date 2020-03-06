function plot_LKCestim( Y, figsetting, outputname )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots LCK estimator results
%
% Input: 
%   Y (structure of input data):
%           include name of method, estimation mean and standard deviation
%   fig_setting (structure): 
%           prespecified figure setting
%   outputname (string):
%           path for saving figure
% Output:
%        LKC_est (array of dimension Dx1): containing the estimated LKC
%                                 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods = Y.methods;
hatn    = Y.hatn;
hatstd  = Y.hatstd;
trueLKC = Y.trueLKC;

xvec  = figsetting.xvec;
sfont = figsetting.sfont;
addf  = figsetting.addf;
xtickcell = figsetting.xtickcell;
ylabels   = figsetting.ylabel;
titles    = figsetting.title;
legendon  = figsetting.legendon;
save = 1;

if ~exist( 'outputname', 'var' )
    save = 0;
end

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
for m = 1:length(methods)
    ErrorBar = errorbar(xvec + 2*(m-2), squeeze(hatn(m,:)),...
                            squeeze(hatstd(m,:)), 'o');
    set( ErrorBar,'Color',figsetting.colMat(m,:),'LineWidth', 2 )   
end
% Plot the true value
plot([xvec(1)-10 xvec(end)+10],[trueLKC trueLKC],'k', 'LineWidth', 2 ), hold on
% Add title
h = title( titles, 'fontsize', sfont+5); set(h, 'Interpreter', 'latex');
% Modify gloabal font size for this plot
set(gca,'FontSize', sfont)

% Change axis style
xlim([xvec(1)-10 xvec(end)+10])
xticks(xvec)
xticklabels(xtickcell)

h = xlabel( 'Sample Size [{\it N}\,]', 'fontsize', sfont+addf );  set(h, 'Interpreter', 'latex');
h = ylabel( ylabels, 'fontsize', sfont+addf ); set(h, 'Interpreter', 'latex');
% Add legend
if(legendon==1)
    h = legend(methods, 'Location', 'northwest' ); set(h, 'Interpreter', 'latex');
    legend boxoff
end

if(save==1)
    print( outputname, '-dpng')
end

 
