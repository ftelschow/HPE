%% Load data
clear all
close all
pic_path = '/Users/wenyilin/Documents/MATLAB/HPE/results/pics/';
res_path = '/Users/wenyilin/Documents/MATLAB/HPE/results/LKCestimation3D/';
load(strcat(res_path,'SimResults_Isotropic_T50_D3_Nfields4000_methods_Forman_Friston_HPE_spm_new.mat'));

%% Figure setup
Vibrant    = [[0 119 187];... % blue
              [51 187 238];...% cyan
              [0 153 136];... % teal
              [238 119 51];...% orange
              [204 51 17];... % red
              [238 51 119];...% magenta
              [187 187 187]...% grey
              ]/255;
          
scale = 0.9;
sfont = 20;
addf  = 5;
WidthFig   = 1300;
HeightFig  = WidthFig * scale;

figsetting = struct();
figsetting.xvec       = (1:length(Nsubj)) *15; %nvec;
figsetting.xtickcell  = cellstr(num2str(Nsubj'));
figsetting.colMat = Vibrant([1 3 4 5 6],:);
figsetting.sfont = sfont;
figsetting.addf = addf;

%% Plotting LKC estimator
D = 3;
T = 50;
methods = Isotropic.method;
FWHM = [3 6 12 15];
LKC_estim_hatmean = zeros(length(methods), D,length(Nsubj), length(FWHM));
LKC_estim_hatstd = zeros(length(methods),D, length(Nsubj), length(FWHM));
LKC_estim_trueLKC = Isotropic.trueLKC;

for m = 1:length(methods)
    LKC_estim_hatmean(m,:,:,:) = Isotropic.(methods{m}).LKChatmean;
    LKC_estim_hatstd(m,:,:,:) = Isotropic.(methods{m}).LKChatsd;
end

for f = 1 : length(FWHM)
    % generate subplot
    figure(f), clf, hold on;
    set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
    set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig]);
    outputname = strcat(pic_path, "Plots_SimResults_Isotropic_T",...
             num2str(T), "_D", num2str(D),"_FWHM",num2str(FWHM(f)));
    for d = 1:D
        subplot(1,D,d), hold on;
        if(d==1)
            figsetting.legendon = 1;
        else
            figsetting.legendon = 0;
        end
        figsetting.ylabel = strcat("$\mathcal{L}_", num2str(d), "$");
        figsetting.title = strcat('FWHM=',num2str(FWHM(f)));
        LKC_estim = struct();
        LKC_estim.methods = methods;
        LKC_estim.hatn = squeeze(LKC_estim_hatmean(:,d,:,f));
        LKC_estim.hatstd = squeeze(LKC_estim_hatstd(:,d,:,f));
        LKC_estim.trueLKC = LKC_estim_trueLKC(d,f);
        plot_LKCestim(LKC_estim,figsetting);
    end
    set(gcf,'papersize',[12 12*scale])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print( outputname, '-dpng')
    hold off;
end

%% Plotting threshold
D = 3;
T = 50;
thre_estim_hatmean = zeros(length(methods),length(Nsubj), length(FWHM));
thre_estim_hatstd = zeros(length(methods), length(Nsubj), length(FWHM));
thre_estim_trueu = Isotropic.trueu;

for m = 1:length(methods)
    thre_estim_hatmean(m,:,:) = Isotropic.(methods{m}).uhatmean;
    thre_estim_hatstd(m,:,:) = Isotropic.(methods{m}).uhatsd;
end

for f = 1 : length(FWHM)
    % generate subplot
    figure(f), clf, hold on;
    set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
    set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig]);
    outputname = strcat(pic_path, "Plots_threshold_SimResults_Isotropic_T",...
             num2str(T), "_D", num2str(D),"_FWHM",num2str(FWHM(f)));
    %for d = 1:D
%         subplot(1,D,d), hold on;
%         if(d==1)
%             figsetting.legendon = 1;
%         else
%             figsetting.legendon = 0;
%         end
    figsetting.legendon = 1;
    figsetting.ylabel = strcat("$\mathcal{Threshold}");
    figsetting.title = strcat('FWHM=',num2str(FWHM(f)));
    thre_estim = struct();
    thre_estim.methods = methods;
    thre_estim.hatn = squeeze(thre_estim_hatmean(:,:,f));
    thre_estim.hatstd = squeeze(thre_estim_hatstd(:,:,f));
    thre_estim.trueLKC = thre_estim_trueu(f);
    plot_LKCestim(thre_estim,figsetting);
%     end
    set(gcf,'papersize',[12 12*scale])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print( outputname, '-dpng')
    hold off;
end

%% Plotting FWHM
D = 3;
T = 50;
fwhm_estim_hatn = zeros(length(methods),length(Nsubj), length(FWHM),Msim);
fwhm_estim_trueu = FWHM;
for m = 1:length(methods)
    LKC_hatn = Isotropic.(methods{m}).LKChatn;
    for i = 1:length(Nsubj)
        for j = 1:length(FWHM)
            for k = 1:Msim
                [~,temp_f,~,~] = LKC_conversion(LKC_hatn(:,i,j,k),D,repmat(T, [1 D]),"LKC");
                fwhm_estim_hatn(m,i,j,k) = temp_f;
            end
        end
    end
end
fwhm_estim_hatmean = mean(fwhm_estim_hatn,4);
fwhm_estim_hatstd = std(fwhm_estim_hatn,0,4);

for f = 1 : length(FWHM)
    % generate subplot
    figure(f), clf, hold on;
    set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
    set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig]);
    outputname = strcat(pic_path, "Plots_fwhm_SimResults_Isotropic_T",...
             num2str(T), "_D", num2str(D),"_FWHM",num2str(FWHM(f)));
    figsetting.legendon = 1;
    figsetting.ylabel = strcat('$\mathcal{FWHM}$');
    figsetting.title = strcat('FWHM=',num2str(FWHM(f)));
    fwhm_estim = struct();
    fwhm_estim.methods = methods;
    fwhm_estim.hatn = squeeze(fwhm_estim_hatmean(:,:,f));
    fwhm_estim.hatstd = squeeze(fwhm_estim_hatstd(:,:,f));
    fwhm_estim.trueLKC = fwhm_estim_trueu(f);
    plot_LKCestim(fwhm_estim,figsetting);
%     end
    set(gcf,'papersize',[12 12*scale])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print( outputname, '-dpng')
    hold off;
end
