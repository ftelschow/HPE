%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%     Estimates L1 and L2 for N realisations of a Gaussian field over a
%%%     square
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code estimates the LKC according to Taylor Worsley 2007.
% Input: 
% Res (array of dimension KxLxN): N observations of residual fields
% over an KxL-square.
%
% center (boolean): If 'true' the sample mean of the input field 'Res' will be
% subtracted and the variance will be corrected by sqrt( N/(N-1) ).
%
% normalize (boolean): If 'true' estimate the standard deviation of the
% field 'Res' and divide by it.
%
%
% Output:
%        LKC (array of dimension 2x1): containing the estimated Lipschitz
%                                      killing curvatures LKC1 and LKC2
function LKC = LKCestim_warp( f, center, normalize )

sf = size( f );
D  = length( sf )-1;
N  = sf( D+1 );

% Initiate output
LKC = zeros( 1, D );

% If necessary center and normalize the input fields
if( center )
    f = f - mean( f, D+1 );
end

if( normalize )
    f = f ./ sqrt( sum( f.^2, D+1 ) );
else
    f = f / sqrt(N);
end

% Calculate the LKCs
switch D
    case 1
        % note that this implementation does not take disconnected domains
        % into account yet!
        f = f(mask, :);
        LKC(1) = sum( sqrt( sum( diff( f, 1, 1).^2, 2 ) ) );
    case 2
        % Compute length of horizontal/vertical edges, the
        % triangulation of the mask is given by the following pattern for
        % 3 x 3 square. This is a special choice and could be replaced by
        % any other triangulation
        %
        % o---o---o---o---o
        % | / | / | / | / |
        % o---o---o---o---o
        % | / | / | / | / |
        % o---o---o---o---o
        edges_vert = diff( f, 1, 1 );
        edges_horz = diff( f, 1, 2 );
        l_edges_vert = sqrt(sum( edges_vert.^2, 3 ));
        l_edges_horz = sqrt(sum( edges_horz.^2, 3 ));

        % Compute estimator L1 = sum_edges mu_1(edges) - sum_triangles mu_1(triangles)
        % Note that the interior cancels out!
        LKC(1) = 0.5*( sum( l_edges_horz(1,:) ) + sum( l_edges_horz(end,:) ) + ...
                       sum( l_edges_vert(:,1) ) + sum( l_edges_vert(:,end) ) ) ;

        % Compute estimator L2 = sum_triangles mu_2(triangles)
        % Taylor (2007) "Detecting Sparse Signals..." p920 eq. (21)
        LKC(2) = sum(sum(...
                    sqrt(...
                            l_edges_vert(:,1:(end-1)).^2 .* l_edges_horz(1:(end-1),:).^2 ...
                            - ( sum( edges_vert(:,1:(end-1),:) .* edges_horz(1:(end-1),:,:), 3 ) ).^2 ...
                         )...
                 )) / 2 + ...
                 sum(sum(...
                    sqrt(...
                            l_edges_vert(:,2:end).^2 .* l_edges_horz(2:end,:).^2 ...
                            - ( sum( edges_vert(:,2:end,:) .* edges_horz(2:end,:,:), 3 ) ).^2 ...
                         )...
                 )) / 2;
end
 