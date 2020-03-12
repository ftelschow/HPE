function [] = Simulate_LKCestim3D_isotropic( Msim , methods, errorproc, identifier )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%    Function to simulate LKC estimators and corresponding thresholds
%%%%    in 3D for isotropic random fields
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Function to simulate LKC estimators and corresponding
%              thresholds in 3D for isotropic random fields
%              For Msim = 1e3 methods = ["spm12", "HPE", "Forman", "Friston"]
%              it takes roughly 18hours
%              For Msim = 1e1 methods = ["bHPE"] it takes ~3hours
%__________________________________________________________________________
% REFERENCES:
%
%__________________________________________________________________________
% AUTHOR: Fabian Telschow (ftelschow@ucsd.edu)
%__________________________________________________________________________
%------ prepare workspace
% load generated profile containing the paths
load('paths.mat')
cd(path_main)

% create matlab paths to search for functions
addpath(genpath(path_main));
addpath(genpath(path_RFTtoolbox));
addpath(genpath(path_spm12));

%------ simulation parameters
Nfields = 4e3;
Nsubj   = [50 100 200];
FWER    = 0.05;

%------ field parameterscd ~
D    = 3;
T    = 50;
Dim  = ones([1 D])*T;
FWHM = [ 3 6 12 15 ];

%------ LKC estimator parameters
% general
uvals = -6:0.01:6; % evaluation of EEC curve for estimation of the thresholds
% number of bootstrap replicates for bHPE
Mboot = 3e3;

%------  general constants
indexD    = repmat( {':'}, 1, D );


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%------------------ Estimate LKCs and thresholds --------------------------
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
if strcmp(errorproc, "isotropic")
    %------ isotropic Gaussian fields
    % initialize buckets for estimators and threshold values
    Isotropic       = struct();
    Isotropic.methodnames = methods;
    Isotropic.FWER  = FWER;
    tic
    for method =  methods
        Isotropic.(method) = struct();
        Isotropic.(method).LKChatn    = zeros( D, length(Nsubj),...
                                                      length(FWHM), Msim);
        Isotropic.(method).LKChatmean = zeros( D, length(Nsubj),...
                                                      length(FWHM));
        Isotropic.(method).LKChatsd   = zeros( D, length(Nsubj),...
                                                      length(FWHM));
        Isotropic.(method).uhatn      = zeros( length(Nsubj),...
                                                      length(FWHM), Msim);
        Isotropic.(method).uhatmean   = zeros( length(Nsubj),...
                                                      length(FWHM));
        Isotropic.(method).uhatsd     = zeros( length(Nsubj),...
                                                      length(FWHM));
    end

    Isotropic.trueLKC = [];
    Isotropic.trueu = [];

    for fk = 1:length(FWHM)
        %------ load precomputed data
        load( fullfile( path_data, strcat("Isotropic_FWHM", num2str(FWHM(fk)), "_T",...
                             num2str(T), "_D", num2str(D), "_Nfields",...
                             num2str(Nfields) )), 'rfs', 'LKC' )
        %------ change name of the field and compute parameters
        Y = rfs;
        clear rfs;
        sY = size(Y);
        D  = length(sY)-1;
        N  = sY(end);

        %------ save true LKC and threshold
        Isotropic.trueLKC = [ Isotropic.trueLKC, LKC ];
        Isotropic.trueu   = [ Isotropic.trueu, get_EECthreshold(...
                                FWER, uvals, Isotropic.trueLKC(:, fk), D ) ];

        %------ prepare input for LKCestim_spm
        if ismember( "spm", methods )
            % clear tmp directory
            cd( path_tmp )
            delete *im.nii
            cd ..

            % save data and mask in .nii files
            path_mask = fullfile( path_tmp, 'mask.nii');
            nii_img   = make_nii(ones(Dim));
            save_nii( nii_img, path_mask );
            clear nii_img f

            % write fields into .nii
            for j = 1:N
                nii_img = make_nii( Y(indexD{:},j) );
                save_nii(nii_img, fullfile( path_tmp, strcat(num2str(j),'im.nii') ));
            end

            % save names from .nii files
            dinfo     = dir(fullfile( path_tmp, '*im.nii' ));
            Yspm.data = fullfile( path_tmp, {dinfo.name} );
            Yspm.path_mask = fullfile( path_tmp, 'mask.nii' );
            Yspm.path_tmp  = path_tmp;
            sY = size(Y);
            Yspm.dim = sY(1:end-1);
            Yspm.D   = length(sY)-1;
            clear sY
        end

         %------ apply different estimators
        if ismember( "spm", methods )
            method = struct();
            method.name = "spm";
            [LKCf, uf] = simulate_LKCThresh( Yspm, method, Nsubj, Msim,...
                                                     FWER, uvals );
            % fill the result structure with statistical descriptors
            Isotropic.spm.LKChatn(:,:,fk,:)  = LKCf;
            Isotropic.spm.LKChatmean(:,:,fk) = mean( LKCf, 3 );
            Isotropic.spm.LKChatsd(:,:,fk)   = std( LKCf, 0, 3 );
            Isotropic.spm.uhatn(:,fk,:)      = uf;
            Isotropic.spm.uhatmean(:,fk)     = mean( uf, 2 );
            Isotropic.spm.uhatsd(:,fk)       = std( uf, 0, 2 );
        end

        if ismember( "HPE", methods )
            method = struct();
            method.name = "HPE";
            [LKCf, uf] = simulate_LKCThresh( Y, method, Nsubj, Msim,...
                                                     FWER, uvals );
            % fill the result structure with statistical descriptors
            Isotropic.HPE.LKChatn(:,:,fk,:)  = LKCf;
            Isotropic.HPE.LKChatmean(:,:,fk) = mean( LKCf, 3 );
            Isotropic.HPE.LKChatsd(:,:,fk)   = std( LKCf, 0, 3 );
            Isotropic.HPE.uhatn(:,fk,:)      = uf;
            Isotropic.HPE.uhatmean(:,fk)     = mean( uf, 2 );
            Isotropic.HPE.uhatsd(:,fk)       = std( uf, 0, 2 );
        end

        if ismember( "bHPE", methods )
            method = struct();
            method.name  = "bHPE";
            method.Mboot = Mboot;
            [LKCf, uf] = simulate_LKCThresh( Y, method, Nsubj, Msim,...
                                                     FWER, uvals );
            % fill the result structure with statistical descriptors
            Isotropic.bHPE.LKChatn(:,:,fk,:)  = LKCf;
            Isotropic.bHPE.LKChatmean(:,:,fk) = mean( LKCf, 3 );
            Isotropic.bHPE.LKChatsd(:,:,fk)   = std( LKCf, 0, 3 );
            Isotropic.bHPE.uhatn(:,fk,:)      = uf;
            Isotropic.bHPE.uhatmean(:,fk)     = mean( uf, 2 );
            Isotropic.bHPE.uhatsd(:,fk)       = std( uf, 0, 2 );
        end

        if ismember( "Forman", methods )
            method = struct();
            method.name  = "Forman";
            [LKCf, uf] = simulate_LKCThresh( Y, method, Nsubj, Msim,...
                                                     FWER, uvals );
            % fill the result structure with statistical descriptors
            Isotropic.Forman.LKChatn(:,:,fk,:)  = LKCf;
            Isotropic.Forman.LKChatmean(:,:,fk) = mean( LKCf, 3 );
            Isotropic.Forman.LKChatsd(:,:,fk)   = std( LKCf, 0, 3 );
            Isotropic.Forman.uhatn(:,fk,:)      = uf;
            Isotropic.Forman.uhatmean(:,fk)     = mean( uf, 2 );
            Isotropic.Forman.uhatsd(:,fk)       = std( uf, 0, 2 );
        end

        if ismember( "Friston", methods )
            method = struct();
            method.name  = "Friston";
            [LKCf, uf] = simulate_LKCThresh( Y, method, Nsubj, Msim,...
                                                     FWER, uvals );
            % fill the result structure with statistical descriptors
            Isotropic.Friston.LKChatn(:,:,fk,:)  = LKCf;
            Isotropic.Friston.LKChatmean(:,:,fk) = mean( LKCf, 3 );
            Isotropic.Friston.LKChatsd(:,:,fk)   = std( LKCf, 0, 3 );
            Isotropic.Friston.uhatn(:,fk,:)      = uf;
            Isotropic.Friston.uhatmean(:,fk)     = mean( uf, 2 );
            Isotropic.Friston.uhatsd(:,fk)       = std( uf, 0, 2 );
        end
        clear Y
    end
    toc

    save( fullfile( path_results, strcat( "SimResults_Isotropic_T",...
                    num2str(T), "_D", num2str(D), "_Nfields",...
                    num2str(Nfields), "_methods_", strjoin(methods, '_'),...
                    identifier,".mat" ) ),'-v7.3', 'Isotropic', 'Nsubj', 'Msim' )
 
else
    %------ scale space Gaussian field
    % initialize buckets for estimators and threshold values
    scale       = struct();
    scale.methodnames = methods;
    scale.FWER  = FWER;
    tic
    for method =  methods
        scale.(method) = struct();
        scale.(method).LKChatn    = zeros( D, length( Nsubj ),...
                                                     Msim);
        scale.(method).LKChatmean = zeros( D, length( Nsubj ) );
        scale.(method).LKChatsd   = zeros( D, length( Nsubj ) );
        scale.(method).uhatn      = zeros( length( Nsubj ), Msim );
        scale.(method).uhatmean   = zeros( length( Nsubj ) );
        scale.(method).uhatsd     = zeros( length( Nsubj ) );
    end

    scale.trueLKC = [];
    scale.trueu = [];

    %------ load precomputed data
    load( fullfile( path_data, 'ScaleSpace_T50_D3_Nfields4000.mat'...
                                 ), 'rfs', 'LKC' );
    %------ change name of the field and compute parameters
    Y = rfs;
    clear rfs;
    sY = size(Y);
    D  = length(sY)-1;
    N  = sY(end);

    %------ save true LKC and threshold
    scale.trueLKC = LKC;
    scale.trueu   = get_EECthreshold(...
                            FWER, uvals, scale.trueLKC, D );

    %------ prepare input for LKCestim_spm
    if ismember( "spm", methods )
        % clear tmp directory
        cd( path_tmp )
        delete *im.nii
        cd ..

        % save data and mask in .nii files
        path_mask = fullfile( path_tmp, 'mask.nii');
        nii_img   = make_nii(ones(Dim));
        save_nii( nii_img, path_mask );
        clear nii_img f

        % write fields into .nii
        for j = 1:N
            nii_img = make_nii( Y(indexD{:},j) );
            save_nii(nii_img, fullfile( path_tmp, strcat(num2str(j),'im.nii') ));
        end

        % save names from .nii files
        dinfo     = dir(fullfile( path_tmp, '*im.nii' ));
        Yspm.data = fullfile( path_tmp, {dinfo.name} );
        Yspm.path_mask = fullfile( path_tmp, 'mask.nii' );
        Yspm.path_tmp  = path_tmp;
        sY = size(Y);
        Yspm.dim = sY(1:end-1);
        Yspm.D   = length(sY)-1;
        clear sY

         %------ apply different estimators
        if ismember( "spm", methods )
            method = struct();
            method.name = "spm";
            [LKCf, uf] = simulate_LKCThresh( Yspm, method, Nsubj, Msim,...
                                                     FWER, uvals );
            % fill the result structure with statistical descriptors
            scale.spm.LKChatn    = LKCf;
            scale.spm.LKChatmean = mean( LKCf, 3 );
            scale.spm.LKChatsd   = std( LKCf, 0, 3 );
            scale.spm.uhatn      = uf;
            scale.spm.uhatmean   = mean( uf, 2 );
            scale.spm.uhatsd     = std( uf, 0, 2 );
        end

        if ismember( "HPE", methods )
            method = struct();
            method.name = "HPE";
            [LKCf, uf] = simulate_LKCThresh( Y, method, Nsubj, Msim,...
                                                     FWER, uvals );
            % fill the result structure with statistical descriptors
            scale.HPE.LKChatn    = LKCf;
            scale.HPE.LKChatmean = mean( LKCf, 3 );
            scale.HPE.LKChatsd   = std( LKCf, 0, 3 );
            scale.HPE.uhatn      = uf;
            scale.HPE.uhatmean   = mean( uf, 2 );
            scale.HPE.uhatsd     = std( uf, 0, 2 );
        end

        if ismember( "bHPE", methods )
            method = struct();
            method.name  = "bHPE";
            method.Mboot = Mboot;
            [LKCf, uf] = simulate_LKCThresh( Y, method, Nsubj, Msim,...
                                                     FWER, uvals );
            % fill the result structure with statistical descriptors
            scale.bHPE.LKChatn    = LKCf;
            scale.bHPE.LKChatmean = mean( LKCf, 3 );
            scale.bHPE.LKChatsd   = std( LKCf, 0, 3 );
            scale.bHPE.uhatn      = uf;
            scale.bHPE.uhatmean   = mean( uf, 2 );
            scale.bHPE.uhatsd     = std( uf, 0, 2 );
        end

        if ismember( "Forman", methods )
            method = struct();
            method.name  = "Forman";
            [LKCf, uf] = simulate_LKCThresh( Y, method, Nsubj, Msim,...
                                                     FWER, uvals );
            % fill the result structure with statistical descriptors
            scale.Forman.LKChatn    = LKCf;
            scale.Forman.LKChatmean = mean( LKCf, 3 );
            scale.Forman.LKChatsd   = std( LKCf, 0, 3 );
            scale.Forman.uhatn      = uf;
            scale.Forman.uhatmean   = mean( uf, 2 );
            scale.Forman.uhatsd     = std( uf, 0, 2 );
        end

        if ismember( "Friston", methods )
            method = struct();
            method.name  = "Friston";
            [LKCf, uf] = simulate_LKCThresh( Y, method, Nsubj, Msim,...
                                                     FWER, uvals );
            % fill the result structure with statistical descriptors
            scale.Friston.LKChatn    = LKCf;
            scale.Friston.LKChatmean = mean( LKCf, 3 );
            scale.Friston.LKChatsd   = std( LKCf, 0, 3 );
            scale.Friston.uhatn      = uf;
            scale.Friston.uhatmean   = mean( uf, 2 );
            scale.Friston.uhatsd     = std( uf, 0, 2 );
        end
        clear Y
    end
    toc

    save( fullfile( path_results, strcat( "SimResults_scale_T",...
                    num2str(T), "_D", num2str(D), "_Nfields",...
                    num2str(Nfields), "_methods_", strjoin(methods, '_'),...
                    identifier,".mat" ) ),'-v7.3', 'scale', 'Nsubj', 'Msim' )
end
