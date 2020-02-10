function LKCspm = LKCestim_spm(Y, df, path_mask, dir_data, del)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the Lipschitz Killing curvatures using spm.
%
% Input: 
%   Y (array of dimension L_1 x...x L_D x N or paths for saving nii files):
%           N observations of residual fields over an L_1 x...x L_D-square.
%   df (integer): 
%           degrees of freedom of the input ( usually df of the residuals
%                                             from the linear model)
%   path_mask (array of dimension L_1 x...x L_D):
%           folder path fto the mask.nii file
%   dir_data (char):
%           folder path for saving data as .nii files. Only required if Y
%           is an array.
%   del (integer)
%           delete all .nii files if del=1. Default del=1.      
% Output:
%        LKC (array of dimension Dx1): containing the estimated Lipschitz
%                                      killing curvatures LKC1 and LKC2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------ Check input and write out .nii files ------------------------
if isnumeric(Y) 
    % get size of the input array, dimension of domain and sample size
    sY = size(Y);
    N  = sY(end);
    D  = length(sY) - 1;
    index = repmat( {':'}, 1, D );
    % write data into nii files
    for j = 1:N
        nii_img = make_nii( Y( index{:}, j ) );
        save_nii( nii_img, strcat( dir_data, '/', num2str(j), 'im.nii' ) );
    end
    dinfo = dir(strcat( dir_data, '/', '*im.nii' ));
    names = {dinfo.name};
    if ~exist('del', 'var')
        del = 1;
    end
    Y = strcat(dir_data,'/',char(names));
elseif ischar(Y)
    N     = size(Y,1);
else
    error('The type of inputs should be either a data array or paths names')
end

%------------ Use spm to estimate the LKCs --------------------------------
[~,~,R] = spm_est_smoothness( Y, path_mask, [N df]);
LKCspm = ( R(2:end) .* sqrt( (4*log(2)).^(1:(size(R,2)-1)) ) )';

%------------ clear directory for nii files -------------------------------
if del == 1 
    delete( strcat( dir_data,'/','*im.nii' ) )
end