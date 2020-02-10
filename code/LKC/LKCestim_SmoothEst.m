function L_hat = LKCestim_SmoothEst( Y, D, mask, voxdim, method )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the Lipschitz Killing curvature using the
% Forman or Friston estimator.
% It translates the SmoothEst function in R package AnalyzeFMRI into
% Matlab.
%
% Input: 
%   Y (array of dimension L_1 x...x L_D x N):
%           N observations of residual fields over an L_1 x...x L_D-square.
%           Especially, it is assumed that the field has empirical mean 0.
%   D (integer):
%           dimension of the domain
%   mask (boolean array of dimension L_1 x...x L_D):
%           indicates which voxels are in the mask, other values will be
%           set ti +oo
%   voxdim ?(Vector of length D): 
%            containing the voxel dimensions.
%   method (string):
%            Could be either forman or friston
% Output:
%        LKC (array of dimension 2x1): containing the estimated Lipschitz
%                                      killing curvatures LKC1 and LKC2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check that method is implemented for dimension D
if(D>4)
    error('The method is currently only implemented for field domains of dimension D<5.')
end

if(length(size(Y)) < D)
    error('The input array Y need to have at least D dimensions.')
elseif(length(size(Y)) > D+2)
    error('The input array Y need to have at most D+2 dimensions.')
end

sY     = size(Y);                   % get size of the input array
dim = sY(1:D);

tmp = zeros([D, sY(4)]);
%% set-up  
for i = 1:sY(4)    
    x = dim(1);
    y = dim(2);
    z = dim(3);
    b1 = zeros(dim+2);
    b2 = zeros(dim+2);
    m1 = zeros(dim+2);
    m2 = zeros(dim+2);
    % X
    b1(2:(x + 1), 2:(y + 1), 2:(z + 1)) = Y(:,:,:,i);
    b2(1:(x)    , 2:(y + 1), 2:(z + 1)) = Y(:,:,:,i);
    m1(2:(x + 1), 2:(y + 1), 2:(z + 1)) = mask;
    m2(1:(x)    , 2:(y + 1), 2:(z + 1)) = mask;
    m3 = (m1 + m2) == 2;
    if method == "Forman"        
      x_v0 = mean((b1.*m1).^2,'all');
      x_v1 = mean(((b1 - b2).* m3).^2,'all');
      xx = -(voxdim(1)^2) / (4 * log(1 - x_v1 / (2 * x_v0)));
    end
    if method == "Friston"
      xx = mean(((b1 - b2).*m3).^2,'all') / (voxdim(1)^2);
    end
    
    b1 = zeros(dim+2);
    b2 = zeros(dim+2);
    m1 = zeros(dim+2);
    m2 = zeros(dim+2);
    % Y
    b1(2:(x + 1), 2:(y + 1), 2:(z + 1)) = Y(:,:,:,i);
    b2(2:(x + 1), 1:(y)    , 2:(z + 1)) = Y(:,:,:,i);
    m1(2:(x + 1), 2:(y + 1), 2:(z + 1)) = mask;
    m2(2:(x + 1), 1:(y)    , 2:(z + 1)) = mask;
    m3 = (m1 + m2) == 2;
    if method == "Forman"        
      y_v0 = mean((b1.*m1).^2,'all');
      y_v1 = mean(((b1 - b2).* m3).^2,'all');
      yy = -(voxdim(2)^2) / (4 * log(1 - y_v1 / (2 * y_v0)));
    end
    if method == "Friston"
      yy = mean(((b1 - b2).*m3).^2,'all') / (voxdim(2)^2);
    end
    
    b1 = zeros(dim+2);
    b2 = zeros(dim+2);
    m1 = zeros(dim+2);
    m2 = zeros(dim+2);
    % Z
    b1(2:(x + 1), 2:(y + 1), 2:(z + 1)) = Y(:,:,:,i);
    b2(2:(x + 1), 2:(y + 1), 1:(z)) = Y(:,:,:,i);
    m1(2:(x + 1), 2:(y + 1), 2:(z + 1)) = mask;
    m2(2:(x + 1), 2:(y + 1), 1:(z)) = mask;
    m3 = (m1 + m2) == 2;
    if method == "Forman"        
      z_v0 = mean((b1.*m1).^2,'all');
      z_v1 = mean(((b1 - b2).* m3).^2,'all');
      zz = -(voxdim(3)^2) / (4 * log(1 - z_v1 / (2 * z_v0)));
    end
    if method == "Friston"
      zz = mean(((b1 - b2).*m3).^2,'all') / (voxdim(3)^2);
    end
    
    if method == "Forman" 
        sigma = [xx yy zz];
    end
        
    if method == "Friston"
        sigma = 1./(2 * [xx yy zz]);
    end
    
    fwhm =  mean((sqrt(8*log(2))) *sqrt(sigma(~isinf(sigma))));
    alpha = 1/(4*(fwhm/(2*sqrt(2*log(2))))^2);
    switch D
       case 1
           tmp(:,i) = (dim(1)-1) * sqrt(2*alpha);
       case 2
           tmp(:,i) = [2*(dim(1)-1); (dim(2)-1)^2] .* [sqrt(2*alpha); 2*alpha] ;
       case 3
           tmp(:,i) = [3*(dim(1)-1); 3*(dim(2)-1)^2; (dim(3)-1)^3] .* [sqrt(2*alpha); 2*alpha; (2*alpha)^(3/2)] ;
    end
end
L_hat = mean( tmp, 2 );
