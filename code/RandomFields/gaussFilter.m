function h = gaussFilter( nu, D, delta )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalized Gauss filter tp produce variance 1 noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Input check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist( 'D', 'var' )
    D = 2;
end
if ~exist( 'delta', 'var' )
    delta = ones( [ 1, D ] );
end

if all( size( nu ) == [ 1 1 ] )
    nu = repmat( nu, [ 1 D ] );
end

%%%% Main function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fix threshold to declare the value to be zero for numerical stability
eps = 1e-10;

% numbers of zeros to be padded to each side (heuristic 4*std)
siz   = ceil( 4 * nu );

switch D
    case 1
        x = -siz:delta( 1 ):siz;
        arg   = -( x.*x / nu( 1 )^2 ) ./ 2;

        h            = exp( arg ) / sqrt( ( 2 * pi )^D * prod( nu ) );
        h( h < eps ) = 0;
        
    case 2
        [x,y] = meshgrid( -siz:delta(1):siz, -siz:delta(2):siz );
        arg   = -( x.*x / nu( 1 )^2 + y.*y / nu( 2 )^2 ) ./ 2;

        h        = exp( arg ) / sqrt( ( 2 * pi )^D * prod( nu ) );
        h(h<eps) = 0;

    case 3
        [x, y, z] = meshgrid( -siz : siz, -siz : siz, -siz : siz );
        arg   = -( x.*x / nu( 1 )^2 + y.*y / nu( 2 )^2 + z.*z / nu( 3 )^2 ) / 2;

        h     = exp(arg) / sqrt( (2*pi)^D * prod( nu ) );
        h(h<eps) = 0;
end
% normalize the kernel to square integrate to 1!
h = h / sqrt( sum( h(:).^2 ) );