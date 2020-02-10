function LKC = LKCestim_iso( data, mask, type, df )
% EST_SMOOTH estimates the smoothness of a process in oned. NEED TO DO THE
% CROSS TERMS!!
%--------------------------------------------------------------------------
% ARGUMENTS
% data      Dim by nsubj, data which has nan where there is missing data.
% mask      A mask of the data which is made of 1s and 0s. 1s for where
%           there is data and nans for where there is no data. Default is
%           taken to be the mask with 1s everywhere.
%--------------------------------------------------------------------------
% OUTPUT
% fwhm_est      An estimate of the fwhm in each of the directions.
% Lambda_est    An estimate of the covariance matrix of the partial
%               derivatives.
% sigma_est     An estimate of the smoothness in terms of sigma.
%--------------------------------------------------------------------------
% EXAMPLES
% Dim = 160;
% nsubj = 20;
% noise = noisegen(Dim, nsubj, 6);
% est_smooth(noise)
%
% Dim = [250,250];
% nsubj = 20;
% noise = noisegen(Dim, nsubj, 10);
% est_smooth(noise)
%
% Dim = [250,250];
% nsubj = 1;
% noise = noisegen(Dim, nsubj, 10);
% est_smooth(noise)
%
% Dim = [250,250];
% nsubj = 100;
% noise = noisegen(Dim, nsubj, 6);
% est_smooth(noise)
%
% Dim = [91,109,91];
% nsubj = 100;
% noise = noisegen(Dim, nsubj, 2);
% est_smooth(noise) %Gets around 2.39
%
% Dim = [91,109,91];
% nsubj = 100;
% noise = noisegen(Dim, nsubj, 3);
% est_smooth(noise) %Gets around 3.28.
%--------------------------------------------------------------------------
% AUTHOR: Sam Davenport & Fabian Telschow.

%% Setup
size_of_data = size(data);
if length(size_of_data) < 5
    nDim = length(size_of_data) - 1;
else
    error('Est_smooth only works in 1,2 and 3 dimensions')
end
% nVox = prod(size_of_data(1:nDim));
nsubj = size_of_data(end);
% Index to remove cases for different dimensions
index = repmat( {':'}, 1, nDim );

%% Masking 
%Set the default mask to be all ones.
if nargin < 2
    if nDim == 1
        mask = ones(size_of_data(1), 1);
    else
        mask = ones(size_of_data(1:nDim));
    end
end

%Make the zero-entries of the mask nan.
mask = zero2nan(mask);

%% Standardize and Mask
mean_over_subjects = mean(data, length(size_of_data) );
nVox = sum(~isnan(mask(:)));

%Subtract the mean and multiply by the mask.
if nDim == 1 
    for I = 1:nsubj
        data(:, I) = (data(:, I) - mean_over_subjects).*mask;
    end
elseif nDim == 2
    if type(1)==1
        for I = 1:nsubj
            data(:,:, I) = (data(:,:, I) - mean_over_subjects).*mask;
        end
    else
         for I = 1:nsubj
            data(:,:, I) = data(:,:, I).*mask;
        end       
    end
elseif nDim == 3
    
    for I = 1:nsubj
        data(:,:,:, I) = (data(:,:,:, I) - mean_over_subjects).*mask;
    end
end

if type(2)==1
    var_est = sum(data(~isnan(data)).^2)/((nVox-1)*(nsubj - 1));
    data = data/sqrt(var_est);
end
%% Estimate Lambda Matrix and FWHMs
Lambda_est = zeros(nDim);
fwhm_est   = zeros(1,nDim);

Xderivmate = diff(data,1,1);
if type(3)==1
    tmp   = Xderivmate(index{:}, 1);
    tmp   = ~isnan(tmp);
    denom = sum(tmp(:))*df;
else
    denom = ((nVox-1)*(nsubj-1));
end

Lambda_est(1,1) = sum(Xderivmate(~isnan(Xderivmate)).^2)/denom;
fwhm_est(1)     = sqrt(4*log(2)/Lambda_est(1,1));
if nDim > 1
    Yderivmate      = diff(data,1,2);
    if type(3)==1
        tmp   = Yderivmate(index{:}, 1);
        tmp   = ~isnan(tmp);
        denom = sum(tmp(:))*df;
    else
        denom = ((nVox-1)*(nsubj-1));
    end

    Lambda_est(2,2) = sum(Yderivmate(~isnan(Yderivmate)).^2)/denom;
    fwhm_est(2)     = sqrt(4*log(2)/Lambda_est(2,2));
end
if nDim > 2
    Zderivmate = diff(data,1,3);
    if type(3)==1
        tmp   = Zderivmate(index{:}, 1);
        tmp   = ~isnan(tmp);
        denom = sum(tmp(:))*df;
    else
        denom = ((nVox-1)*(nsubj-1));
    end

    Lambda_est(3,3) = sum(Zderivmate(~isnan(Zderivmate)).^2)/denom;
    fwhm_est(3)     = sqrt(4*log(2)/Lambda_est(3,3));
end

mfwhm = mean(fwhm_est);
% Transform fwhm_est into LKC, using "Unified univariate and
% multivariate random field theory" keith Worsley (2004),
% Tab.2 and comparing it to
% "Detecting Sparse Signals in Random Fields, With an
% Application to Brain Mapping" Taylor 2007, eq. (10)
dim=size(mask);
LKC =  [ sum(dim-1)/mfwhm prod(dim-1)/mfwhm^2 ].* [(4*log(2))^(1/2) (4*log(2))];

end

%     nZvox = sum(~isnan(Zderivmate(:)));
%     Lambda_est(3,3) = sum(Zderivmate(:).^2)/((nVox-1)*(nsubj-1));
%     Lambda_est(3,3) = sum(Zderivmate(~isnan(Zderivmate)).^2)/((nVox-1)*(nsubj-1));
