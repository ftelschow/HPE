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
% unbiased (boolean): If 'true' use biased corrected estimator of the
% standard deviation.
% avLKC (boolean): If 'true' the average over the two different meshes is
% taken. Default 'true'.
%
% Output:
%        LKC (array of dimension 2x1): containing the estimated Lipschitz
%                                      killing curvatures LKC1 and LKC2
function LKC = LKCestim_warp( Res, center, normalize, unbiased, avLKC )
    % Initiate output
    LKC = zeros( 2, 1 ) ;
    
    % Find size of the residuals
    sZ   = size( Res ) ;
    N    = sZ( end ) ;
    
    % If necessary center and normalize the input fields
    if( center )
        Res = sqrt(N/(N-1)) * (Res - repmat(mean(Res, 3), [1 1 N])) ; 
    end
    
    if( normalize )
        Res = Res ./ reshape( repmat( sqrt( sum(Res.^2,  3) ), 1, N),...
              [size(Res(:,:,1)), N] ) ;
    else
        % If not normalized divide by degrees of freedom in order to
        % approximate the variance correctly
        if( center )
            Res = Res / sqrt( N-1 ) ;
        else
            Res = Res / sqrt( N ) ;
        end
        if( unbiased )
            % Make the standard deviation estimation unbiased, use
            % exp(gammaln(...)) since the Gamma function gamma(...) is
            % manually redefined by Armins code
            Res = Res * ( sqrt( (N-1)/2 ) * exp(gammaln((N-1)/2)) / ...
                  exp(gammaln(N/2)) ) ;
        end
    end
    
    % Compute length of horizontal/vertical edges, the
    % triangulation of the rectangle is given by the following pattern for
    % 3 x 3 square
    %
    % o---o---o---o---o
    % | / | / | / | / |
    % o---o---o---o---o
    % | / | / | / | / |
    % o---o---o---o---o
    edges_vert = sqrt(sum( (Res(2:end,:,:) - Res(1:(end-1),:,:)).^2, 3 )) ;
    edges_horz = sqrt(sum( (Res(:,2:end,:) - Res(:,1:(end-1),:)).^2, 3 )) ;
    % edges_diag = sqrt(sum( (Z(1:(end-1),1:(end-1),:) - Z(2:end,2:end,:)).^2, 3 )) ;
    
    % Compute estimator L1 = sum_edges mu_1(edges) - sum_triangles mu_1(triangles)
    % Note that the interior cancels out!
    LKC(1) = 0.5*( sum(edges_horz(1,:)) + sum(edges_horz(end,:)) + ...
               sum(edges_vert(:,1)) + sum(edges_vert(:,end)) ) ;
    
    if avLKC == 0
        % Compute estimator L2 = sum_triangles mu_2(triangles)
        % Taylor (2007) "Detecting Sparse Signals..." p920 eq. (21)

        LKC(2) = sum(sum( edges_vert(:,1:(end-1)) .* edges_horz(1:(end-1),:) )) / 2 + ...
        sum(sum( edges_vert(:,1:end-1) .* edges_horz(2:end,:) )) / 2 ;
    else
        % Compute estimator L2 = sum_triangles mu_2(triangles)
        LKC(2) = 0.5*( ...
            sum(sum( edges_vert(:,1:(end-1)) .* edges_horz(1:(end-1),:) )) / 2 + ...
            sum(sum( edges_vert(:,2:end)     .* edges_horz(2:end,:)     )) / 2 + ...
            sum(sum( edges_vert(:,1:(end-1)) .* edges_horz(2:end,:)     )) / 2 + ...
            sum(sum( edges_vert(:,2:end)     .* edges_horz(1:(end-1),:) )) / 2 ...
        );
    end
 