function LKC = LKCestim_HermProjExact( Y, D, mask, Mboot, ...
                                             L0, bootsmultiplier, cc, version )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the Lipschitz Killing curvatuversionre using the
% Hermite projection estimator proposed in Schwartzman et al (2020+).
% It uses a fast and exact way to compute the EC curves by looping through
% discrete critical values and computes there contribution to the change in
% Euler characteristic.
%
% Input: 
%   Y (array of dimension L_1 x...x L_D x N x nsim):
%           N observations of residual fields over an L_1 x...x L_D-square.
%           Especially, it is assumed that the field has empirical mean 0.
%   D (integer):
%           dimension of the domain
%   mask (boolean array of dimension L_1 x...x L_D):
%           indicates which voxels are in the mask, other values will be
%           set to +oo
%   Mboot (integer):
%           number of bootstrap processes used for estimation of LKC
%   L0 (integer):
%           Euler characteristic of the mask
%   bootsmultiplier (string)
%           Either "Gaussian", "Multinomial" or "Rademacher". Sets the
%           multipliers in the bootstrap to be either Gaussians or Multinomials.
%   version (string):
%           Either "C" or "matlab". Default is "C", which implments the EC
%           computation in C++ otherwise a matlab only version is used.
%
% Output:
%        LKC (array of dimension Dx1): containing the estimated Lipschitz
%                                      killing curvatures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Depends on:
%   
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check that method is implemented for dimension D
if D > 4
    error( 'The method is currently only implemented for field domains of dimension D<5.' )
end

if length( size( Y ) ) < D
    error( 'The input array Y need to have at least D dimensions.' )
elseif length( size( Y ) ) > D + 2 
    error( 'The input array Y need to have at most D+2 dimensions.' )
end


%%%%%%%%%%%%%%%%%%%%%%% Get constants from the input %%%%%%%%%%%%%%%%%%%%%%
sY     = size( Y );                 % get size of the input array
index  = repmat( {':'}, 1, D );     % get variable domain counter

% get parameter of the field
if length( sY ) > D + 1
    N    = sY( D + 1 ); % number of samples used to estimate the LKC
    nsim = sY( D + 2 ); % number of simulations
else
    N = sY( end );
    nsim = 1;
end

if nargin < 7
    version = "C";
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%% add default values %%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist( 'Mboot', 'var' )
   % default number of bootstrap replicates
   Mboot = 0;
end

if ~exist( 'bootsmultiplier', 'var' )
    bootsmultiplier = "Gaussian";
end

if ~exist( 'u', 'var' )
    % default sampling of R for output of covariances and variances
    u = -8 : 0.01 : 8;
end

% Catch wrong mask input or mask the data, if correct mask exists
if exist( 'mask', 'var' )
    if mask ~= -666
        L0 = EulerChar( mask, 0.5, D );
        if ~all( size( mask ) == sY( 1 : D ) )
            error( 'Incompatible input: The mask needs to have the same size as the first D dimensions of Y.' )
        else
            for i = 1 : sY( D + 1 )
                for j = 1 : nsim
                    tmp = Y( index{:}, i, j );
                    tmp( ~mask ) = -Inf;
                    Y( index{:}, i, j ) = tmp;
                end
            end
        end
    else
        mask = true( sY( 1:D ) );
    end
end
clear i j

if ~exist( 'L0', 'var' )
    L0 = 1;
end

if ~exist( 'cc', 'var' )
    switch D
        case 1
            cc = 2;
        case 2
            cc = 4;
        case 3
            cc = 6;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the LKC output
if Mboot > 1
    L_hat = zeros( [ D, Mboot, nsim ] );
else
    L_hat = zeros( [ D, N, nsim ] );    
end
% get size of residuals
sR = sY( 1 : end-1 );

% scaling vector for integral
p = [ sqrt( 2 * pi ); pi; ( 2 * pi )^( 3 / 2 ) / factorial( 3 ); ...
      ( 2 * pi )^( 4 / 2 ) / factorial( 4 ) ];

% Compute LKCs depending on method
for j = 1 : nsim % start loop over simulation
    % get the residuals for the j-th simulation
    R = Y( index{:}, :, j );
    
    if( Mboot > 1 )
        % Get weights for the bootstrap process. Defalut is Gaussian
        % weights
        if strcmp( bootsmultiplier, "Multinomial" )
            multiplier = mnrnd( N, ones( [1, N] ) / N, Mboot )' - 1;
        elseif strcmp( bootsmultiplier, "Naive" )
            multiplier = mnrnd( N, ones( [1, N] ) / N, Mboot )';

        elseif strcmp( bootsmultiplier, "Rademacher" )
            multiplier = randi(2,[N,Mboot])*2-3;
        else
            % Get weights for the multiplier bootstrap
            multiplier = normrnd( 0, 1, [ N, Mboot ] );            
        end
        
        % reshape and and standardize the field, such that it has variance
        % 1
        R = reshape( R, prod( sR(1:D) ), N );
        % normalize the residuals
        R = R ./ sqrt(sum( R.^2, 2 ));
        
          
        for i = 1:Mboot
            % get the bootstrapped process
            mR = reshape( R * multiplier( :, i ), sR );

            % Get the EC stepfunctions
            EC = EulerCharCrit(mR, mask, L0, cc, version);
            EC = EC{ 1 };

            % Get LKC by integrating the Euler Char curves against the Hermite
            % polynomials
            v =  EC( 2:end-1, 1 )';
            b = -diff( EC( 2:end, 2 ) );
            H =  [ v; ( v.^2 - 1 ); ( v.^3 - 3*v ); ( v.^4 - 6*v.^2 + 3 ) ];    % Hermite polynomials

            L_hat( :, i, j ) = p( 1:D ) .* ( H( 1:D, : ) * b );
        end
    else
        % Get the EC stepfunctions
        ECall = EulerCharCrit( R, mask, L0, cc, version );

        for i = 1:N  % start loop over realisations (bootstrap or normal)
            % Get LKC by integrating the Euler Char curves against the Hermite
            % polynomials
            v =  ECall{ i }( 2:end-1, 1 )';
            b = -diff( ECall{ i }( 2:end, 2 ) );
            H =  [ v; ( v.^2 - 1 ); ( v.^3 - 3 * v ); ( v.^4 - 6 * v.^2 + 3 ) ];    % Hermite polynomials

            L_hat( :, i, j ) = p( 1:D ) .* ( H( 1:D, : ) * b );
        end % end loop over realisations
    end
end % end loop over simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stat summary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Summary estimators: LKCs
L_hatn    = permute( mean( L_hat, 2 ), [ 1 3 2 ] );
Sigma_hat = zeros( D, D, nsim );
L_se_hat  = zeros( D, nsim );

for j = 1:nsim
    Sigma_hat( :, :, j ) = cov( L_hat( :, :, j )' );
    L_se_hat( :, j ) = sqrt( diag( Sigma_hat( :, :, j ) ) / size( L_hat, 2 ) );
end
L_conf_hat = cat( 3, L_hatn - 1.96 * L_se_hat, L_hatn + 1.96 * L_se_hat );

% Summarize output
LKC  = struct( 'hat1', L_hat, 'hatn', L_hatn, 'Sigma_hat', Sigma_hat, ...
               'se_hat', L_se_hat, 'conf_hat', L_conf_hat );
return