function [Y, mY, stdY] = standardize(Y, demean, destd)
%__________________________________________________________________________
% Standardizes a sample
% Input:
%   Y       an array of size d1 x ... dl x N x M
%   demean  boolean; 1 if sample mean should be subtracted
%   destd   boolean; 1 if sample std should be divided
%
% Output:
%   Y     array  of size d1 x ... dl x N x M containing the normalized
%         data
%   mY    an array of size d1 x ... dl x M containing the mean of the
%         data
%   stdY  an array of size d1 x ... dl x M containing the standard deviation
%         of the data
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
%__________________________________________________________________________
%
% Start of function
%
switch nargin
    case 1
        demean = 1;
        destd  = 1;
    case 2
        destd  = 0;
end

dim  = size(Y);
l    = length(dim); 
N    = dim(end);

% subtract mean and make sure that variance is still 1
if( demean==1 )
    mY = squeeze( mean(Y,l) );
    Y  =  sqrt( N/(N-1) ) * (Y - mY ) ;
end

% make empirical variance equal to 1
if( destd==1 )
    stdY = squeeze( std(Y,0,l) );
    Y    = Y ./ stdY;
end