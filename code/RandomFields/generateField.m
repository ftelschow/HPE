function [ f, LKC ] = generateField( n, nsim, dim, TYPE, params )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate random fields
% Input:
%   n:      Sample size
%   nsim:   Number of simulations
%   dim:    Dimension of domain (including scale)
%   TYPE:   Type of field ( 'isotropic' or 'scale-space' )
%   params: nu (scalar) for isotropic, params (vector) for scale-space
%
% Output:
%   f:   array of size dim x n x nsim
%   LKC: theoretical isotropic LKCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% constants from input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = length( dim );

%%%% Check input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp( TYPE, 'isotropic' ) && ~strcmp( TYPE, 'scale-space' )...
        && ~strcmp( TYPE, 'nongauss' )
    error( [ TYPE, ' not implemented' ] )
end

if strcmp( TYPE, 'isotropic' ) && ( D ~= 1 && D ~= 2 && D ~= 3 )
    error( 'dim must be of length smaller than 4!' );
end

if strcmp( TYPE, 'scale-space' ) && ( D ~= 2 )
    error( 'N must be 2' );
end

if ~exist( 'L', 'var' )
    L = 50;
end

if ~exist( 'params', 'var' ) && strcmp( TYPE, 'isotropic' )
    params = 5;
end

if ~exist( 'params', 'var' ) && strcmp( TYPE, 'scale-space' )
    params = 4:.2:40;
end

%%%% main functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h   = gaussFilter( params );
pad = ceil( 4 * params );

%%%% Generate fields
if ~strcmp( TYPE, "scale-space" )
switch D
    case 1
        % generate voxelwise noise
        if strcmp( TYPE, "isotropic" )
            w = randn( [ dim + repmat( 2 * pad, 1, D ), n, nsim ] );
        else
            df = 3;
            w = ( chi2rnd( df, [ dim + repmat( 2 * pad, 1, D ), n, nsim ] ) - df )...
                                    / sqrt( 2 * df );
        end

        % convolve voxelwise noise with the filter and cut to valid range
        f = convn( w, h, 'same' );
        f = f( pad+1:end-pad, :, : );
        
    case 2
        % generate voxelwise noise
        if strcmp( TYPE, "isotropic" )
            w = randn( [ dim + repmat( 2 * pad, 1, D ), n, nsim ] );
        else
            df = 3;
            w = ( chi2rnd( df, [ dim + repmat( 2 * pad, 1, D ), n, nsim ] ) - df )...
                                    / sqrt( 2 * df );
        end

        % convolve voxelwise noise with the filter and cut to valid range
        f = convn( w, h, 'same' );
        f = f( pad+1:end-pad, pad+1:end-pad, :, : );
    case 3
        f = NaN * zeros( [ dim n nsim ] );
        % loop over simulations to save memory
        for nn = 1:nsim
            % generate voxelwise noise
            if strcmp( TYPE, "isotropic" )
                w = randn( [ dim + repmat( 2 * pad, 1, D ) n ] );
            else
                df = 3;
                w = ( chi2rnd( df, [ dim + repmat( 2 * pad, 1, D ) n ] ) - df )...
                                        / sqrt( 2 * df );
            end
            
            tmp = convn( w, h, 'same' );
            f( :, :, :, :, nn) = tmp( pad+1:end-pad, pad+1:end-pad, pad+1:end-pad, : );
            clear tmp
        end
    
end
end

if D==2 && strcmp(TYPE, "scale-space")
    b = 5;
    % Generate fields
    whitenoise_large = randn(L + 2*b*params(end), n, nsim);
    f = zeros(L, length(params), n, nsim);
    for k = 1:length(params)
        whitenoise = whitenoise_large((1+round(b*params(end)-b*params(k))):(L+round(b*params(end)+b*params(k))), :, :);
        x = (-b*params(k):b*params(k))';
        w = params(k)^(-1/2)*pi^(-1/4)*exp(-x.^ 2 / (2 * params(k)^2));
        z_params = convn(whitenoise, w, 'valid');
        f(:,k,:,:) = z_params;
    end
end


%%%% compute the true isotropic LKC
if strcmp( TYPE, 'isotropic' ) || strcmp( TYPE, 'nongauss' )
    alpha = 1 / ( 4 * params^2 );
    switch D
        case 1
            LKC = ( L - 1 ) * sqrt( 2 * alpha );
        case 2
            LKC = [ 2 * ( L - 1 ); ( L - 1 )^2 ]...
                  .* [ sqrt( 2 * alpha ); 2 * alpha ] ;
        case 3
            LKC = [ 3 * ( L - 1 ); 3 * ( L - 1 )^2; ( L - 1 )^3 ]...
                  .* [ sqrt( 2 * alpha ); 2 * alpha; ( 2 * alpha )^( 3/2 ) ] ;
    end
else
    if D==2
        % LKCs (Siegmund 1995)
        NN = dim(1) - 1; % Dimension of domain (not including scale)
        kappa = 0.5;
        lambda = 0.5;

        params = params;
        L2 = (L-1) * (params(1)^(-1) - params(end)^(-1)) * sqrt(lambda*kappa);
        L1 = (L-1)/2 * (params(1)^(-1) + params(end)^(-1)) * sqrt(lambda) ...
                + 1 * log(params(end)/params(1)) * sqrt(kappa);
        LKC = [L1; L2] ;
    else
        LKC = NaN;
        disp("true LKC can not be computed for these fields")
    end
end