function [LKC_est,thr_est] = simulate_LKCThresh( Y, method, Nsubj, Msim,...
                                                 FWER, uvals )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function simulates the LKCs and the thresholds
%
% Input: 
%   Y (structure array of input data):
%           N observations of residual fields over an L_1 x...x L_D-square.
%           Especially, it is assumed that the field has empirical mean 0.
%   method (structure array):
%           'name' contains the name of the method
%           'param1',...'paramN' are some parameters of the method
%   method_params (:
%           containing the parameters for the method        
%   D (integer):
%           dimension of the domain
%   Nsubj (integer):
%           number of subject
%   Msim (integer):
%           number of simulation
%   FWER (numeric):
%           error rate
%   uvals (numeric):
%           EEC thresholds
% Output:
%        LKC_est (array of dimension Dx1): containing the estimated LKC
%        thr_est (array of dimension Dx1): containing the estimated thresholds                              
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name = method.name;
% parameters obtained from the input field
if ~strcmp( name, 'spm' )
    sY   = size(Y);
    N    = sY(end);
    Dim  = sY(1:end-1);
    D    = length(sY) - 1;
    mask = ones( Dim ); 
else
    N       = length(Y.data);
    Dim     = Y.dim;
    D       = Y.D;
    m_nii   = make_nii(ones(Dim));
    save_nii( m_nii, Y.path_mask );
end

indexD   = repmat( {':'}, 1, D );

% output bucket
LKC_est  = zeros( D, length( Nsubj ), Msim );
thr_est  = zeros( length( Nsubj ), Msim );

if strcmp( name, 'HPE' )
    % simulate the estimator for different number of subjects
    for nsubj = 1:length( Nsubj ) % loop over number of subjects
        for m = 1:Msim % loop over number of simulations
            % get random indices to construct a sample
            index = randsample( 1:N, Nsubj(nsubj) );
            % estimate the LKC and threshold from the sample
                tmp = LKCestim_HermProjExact( Y( indexD{:}, index ),...
                                        D, mask, 1, [-1 1], 1, "C" );
                LKC_est(:,nsubj,m)  = tmp.hatn;
        end
    end
elseif strcmp( name, 'bHPE' )
    % simulate the estimator for different number of subjects
    for nsubj = 1:length( Nsubj ) % loop over number of subjects
        for m = 1:Msim % loop over number of simulations
            % get random indices to construct a sample
            index = randsample( 1:N, Nsubj(nsubj) );
            % estimate the LKC and threshold from the sample
                tmp = LKCestim_HermProjExact( Y( indexD{:}, index ),...
                                        D, mask, method.Mboot, [-1 1], 1, "C" );
                LKC_est(:,nsubj,m)  = tmp.hatn;
        end
    end
elseif strcmp( name, 'spm' )
    % simulate the estimator for different number of subjects
    for nsubj = 1:length( Nsubj ) % loop over number of subjects
        for m = 1:Msim % loop over number of simulations
            % get random indices to construct a sample
            index = randsample( 1:N, Nsubj(nsubj) );
            % estimate the LKC and threshold from the sample
                LKC_est(:,nsubj,m)  = LKCestim_spm( char(Y.data( index )), Nsubj(nsubj),...
                                                    Y.path_mask, 0, 0 );
        end
    end
elseif strcmp(name,'Friston') || strcmp(name,'Forman')
    % simulate the estimator for different number of subjects
    for nsubj = 1:length( Nsubj ) % loop over number of subjects
        for m = 1:Msim % loop over number of simulations
            % get random indices to construct a sample
            index = randsample( 1:N, Nsubj(nsubj) );
            % estimate the LKC and threshold from the sample
                LKC_est(:,nsubj,m) = LKCestim_SmoothEst( Y( indexD{:}, index ),...
                                                        D, mask, ones([1 D]), name);
        end
    end
else
        error('The method is not supported with this function.')
end
    % compute estimate of the FWER threshold
    thr_est(nsubj,m)  = get_EECthreshold( FWER, uvals,...
                                    LKC_est(:, nsubj, m) );
end