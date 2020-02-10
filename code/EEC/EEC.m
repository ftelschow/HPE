function [ EEC_vals, EEC_se_hat, C ] = EEC( uvals, LKC, LKC0, type, se )
% EEC(u, LKC, type)
% calculates the expected Euler characteristic at all values uvals
%--------------------------------------------------------------------------
% ARGUMENTS
% uvals	an Nuvals length vector whose entries are the real numbers, where
%       the EEC is evaluated
% LKC   an D+1 length vector containing the Lipschitz Killing curvatures.
% type	string specifying the type of random field the EEC is computed from.
%       Default option is "gaussian". (later "t","F") will follow.
%--------------------------------------------------------------------------
% OUTPUT
% EEC_vals    an Nuvals length vector whose entries are the EEC evaluated
%             at the locations specified by uvals
%--------------------------------------------------------------------------
% EXAMPLES  
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow
if nargin < 3
    LKC0   = 1;
    type = "gaussian";
end
if nargin < 4
    type = "gaussian";
end

D = length(LKC);

if strcmp( type, "gaussian" )
    % matrix containing the Hermite polynomials evaluated on uvals
    H      = [ ones( 1, length(uvals) ); uvals; ( uvals.^2 - 1 ) ]';
    % matrix containing the  EC densities evaluated on uvals
    rho    = H .* [ exp( -uvals.^2/2 ) / (2*pi)^(2/2);...
                    exp( -uvals.^2/2 ) / (2*pi)^(3/2); ...
                    exp( -uvals.^2/2 ) / (2*pi)^(4/2) ]';
         
    EEC_vals = ( 1 - normcdf( uvals ) )' * LKC0 + rho(:, 1:D) * LKC(1:D);
    
    % compute estimated se of EEC estimator and standard error
    if isfield( se, 'LKCcov' ) && isfield( se, 'N' )
        C = rho( :, 1:D ) * se.LKCcov * rho( :, 1:D )';
        EEC_se_hat = sqrt( diag(C) / se.N );
    else
        C = NaN;
        EEC_se_hat = NaN;
    end
end