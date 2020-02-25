function [sigma,fwhm,LKC,resel] = LKC_conversion(x, D, dim, TYPE)
% A conversion function between sigma and fwhm, generating LKC and resel in isotropic GRF
% Input:
%   x:      input value
%   D:      dimensions of the domain
%   dim:    size of domain
%   TYPE:   input variable type

if ~exist('D'), D = 1; end
if (length(dim) ~= D) error('The dimension of the domain should equal D'); end
if ~exist('TYPE'), TYPE = 'fwhm'; end

switch TYPE
    case 'sigma'
       sigma = x;
       fwhm =  (sqrt(8*log(2))) *sigma;
       alpha = 1/(4*(fwhm/(2*sqrt(2*log(2))))^2);
       switch D
           case 1
               LKC = (dim(1)-1) * sqrt(2*alpha);
           case 2
               LKC = [2*(dim(1)-1); (dim(2)-1)^2] .* [sqrt(2*alpha); 2*alpha] ;
           case 3
               LKC = [3*(dim(1)-1); 3*(dim(2)-1)^2; (dim(3)-1)^3] .* [sqrt(2*alpha); 2*alpha; (2*alpha)^(3/2)] ;
       end
       resel = (LKC' ./ sqrt((4*log(2)).^(1:D)))';
    case 'fwhm'
        fwhm = x;
        sigma = fwhm / (sqrt(8*log(2)));
        alpha = 1/(4*(fwhm/(2*sqrt(2*log(2))))^2);
        switch D
           case 1
               LKC = (dim(1)-1) * sqrt(2*alpha);
           case 2
               LKC = [2*(dim(1)-1); (dim(2)-1)^2] .* [sqrt(2*alpha); 2*alpha] ;
           case 3
               LKC = [3*(dim(1)-1); 3*(dim(2)-1)^2; (dim(3)-1)^3] .* [sqrt(2*alpha); 2*alpha; (2*alpha)^(3/2)] ;
       end
        resel = (LKC' ./ sqrt((4*log(2)).^(1:D)))';
    case 'LKC'
        LKC = x;
        resel = LKC ./ sqrt((4*log(2)).^(1:D));
        switch D
           case 1
               alpha = (LKC / (dim(1)-1))^2/2;
           case 2
               alpha = [[LKC(1); LKC(2)]./[2*(dim(1)-1); (dim(2)-1)^2]].^[2; 1]/2;
           case 3
               alpha = [[LKC(1); LKC(2);LKC(3)]./[3*(dim(1)-1); 3*(dim(2)-1)^2; (dim(3)-1)^3]].^[2; 1;2/3]/2;
        end
        fwhm = mean((2*sqrt(2*log(2)))./sqrt(4*alpha));
        sigma = fwhm / (sqrt(8*log(2)));
    case 'resel'
        resel = x;
        LKC = resel .* sqrt((4*log(2)).^(1:D));
        switch D
           case 1
               alpha = (LKC / (dim(1)-1))^2/2;
           case 2
               alpha = [[LKC(1); LKC(2)]./[2*(dim(1)-1); (dim(2)-1)^2]].^[2; 1]/2;
           case 3
               alpha = [[LKC(1); LKC(2);LKC(3)]./[3*(dim(1)-1); 3*(dim(2)-1)^2; (dim(3)-1)^3]].^[2; 1;2/3]/2;
        end
        fwhm = mean((2*sqrt(2*log(2)))./sqrt(4*alpha));
        sigma = fwhm / (sqrt(8*log(2)));        
end

        

