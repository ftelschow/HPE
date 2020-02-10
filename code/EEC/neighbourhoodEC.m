function dEC = neighbourhoodEC( Z, pos, D )

% finds the EC of a neighbourhood excursion. 
% Input:
%  Z:       field over a domain in R^2/3
%  pos:     coordinates of the point
% Output:
%  values of i-th smallest value of the points in mask
%
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 10/05/2018
%__________________________________________________________________________
if nargin < 3
    D   = length(size(Z));
end

NumberCrits = size(pos,1);
dEC = zeros(NumberCrits,2);

% Z = padarray(Z, ones([1,D]), -Inf);
pos = pos+1;

% find the change of EC at the crits
switch D
    case 2
        for k = 1:NumberCrits
            dEC(k,1) = Z(pos(k,1),pos(k,2));
            % cutting the circle at Z(pos(k,1)-1,pos(k,2), i.e.
            %
            % + - +       
            % - 0 -   =>  - + - + - - - +
            % + - -
            circle = [ -Inf; Z(pos(k,1)-1,pos(k,2)); Z((pos(k,1)-1):(pos(k,1)+1), pos(k,2)+1);...
                       Z(pos(k,1)+1,pos(k,2)); Z((pos(k,1)+1):-1:(pos(k,1)-1), pos(k,2)-1); -Inf] > dEC(k,1);              
            % compute EC of circle using gluing and 1D line EC (Note
            % that we subtract also an additional -1, since we assume
            % in each of the circles appears thecritical point appears the pattern - + -              
            dEC(k,2) = sum(diff(circle)==1)-circle(2)*circle(end-1) -1;
        end
    case 3
    
end