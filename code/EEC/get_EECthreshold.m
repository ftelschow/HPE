function [u_FWER, u_FWER_se] = get_EECthreshold( FWER, uvals, LKC, LKC0,...
                                                 type, se)
% get_EECthreshold(u, LKC, type)
% calculates the EECthreshold.
%--------------------------------------------------------------------------
% ARGUMENTS
% FWER  an number between [0 1] controlling the FWER
% uvals	an Nuvals length vector whose entries are the real numbers, where
%       the EEC is evaluated
% LKC   an D+1 length vector containing the Lipschitz Killing curvatures.
% se    structure containing covariance matrix of LKC estimators and sample
%       size, i.e., fields 'LKCcov' and 'N'. 
% type	string specifying the type of random field the EEC is computed from.
%       Default option is "gaussian". (later "t","F") will follow.
%--------------------------------------------------------------------------
% OUTPUT
% u_FWER  threshold controlling the FWER using the EEC as approximation for
%         the max distribution
%--------------------------------------------------------------------------
% EXAMPLES  
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow
if nargin < 4
    LKC0   = 1;
    type = "gaussian";
end
if nargin < 5
    type = "gaussian";
end

if strcmp(type, "gaussian")
    [ hatEEC, ~, C ] = EEC( uvals, LKC, LKC0, type, se );
    [ind, u0] = crossing( hatEEC, uvals, FWER, 'linear' );
    u_FWER  = max(u0);
end

% compute the variability of the threshold, if covariance matrix of LKCs is
% provided
if isfield(se,'LKCcov')
    tau       = [ 0; diff( hatEEC ) ] ./ diff(uvals);
    u_FWER_se = sqrt( C( max(ind), max(ind ) ) ...
                            / ( se.N * tau(max(ind))^2 ) );
else
    u_FWER_se = NaN;
end