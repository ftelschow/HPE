function EC = EulerCharCrit(f, mask, L0, cc, version)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the Euler characteristic of excursion sets
%
% Mandatory inputs:
%   f (array K_1, x ... x K_D x N): field the Euler characteristic curves
%                                   needs to be computed for
%
% Optional inputs:
%   mask (array K_1, x ... x K_D): boolean array having TRUE, for the values
%                                  voxels to be considered belonging to the
%                                  mask. (default ones([1, d]))
%   L0 (integer): Euler characteristic of mask. (default 1)
%   version (string): default "C" uses a C implementation for speed up.
%                     "matlab" uses a matlab only implementation.
%
% Output:
%   EC curve is computed for all N fields f. The 2 x N_crit array contains
%   in the first column the heights of the critical values and in the
%   second column the Euler characterisitc of the upper level set.
%   4-connectivity for 2D fields and 6-connectivity for 3D fields is used.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Preliminary constants, input check etc.
% Get constants from the input field
sf = size( f );
if  ~exist( 'mask', 'var' )
    D = length( sf( 1:end-1 ) );
else
    D = length( size( mask ) );
end

% introduce index for accessing fields of different dimensions
index  = repmat( {':'}, 1, D );

% Initialize the EC curve cells, one for each realisation of the field
if length( sf ) > D
    nEC = sf( end );
    EC = cell( [ 1 nEC ] );
else
    nEC = 1;
    EC  = cell( 1 );
end

% Set default parameter, if not provided
if ~exist( 'mask', 'var' )
    % masks which voxels should be considered
    mask = -666;
elseif ~all( size( mask ) == sf( 1:D ) )
    error( "Please, specify an input mask, which has the same dimension as the domain of the data" )
end

if ~exist( 'L0', 'var' )
    if mask ~= -666
        L0 = EulerChar( mask, 0.5, D );
    else
        L0 = 1;
    end
end

if ~exist( 'version', 'var' )
    version = "C";
end

if D > 3
    error( "The method is unfortunately until now only implemeted for D < 3." )
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


%%%%%%%%% Beginning of the algorithm
% set all values outside the mask to -Inf
if mask(1) ~= -666
    for kk = 1:size( f, D+1 )
        tmp = f( index{:}, kk );
        tmp( ~mask ) = -Inf;
        f( index{:}, kk ) = tmp;
    end
end

% Pad large negative values to the array, which is required for 
% EulerCharCrit_c input
f_tmp = f;
if length( sf ) == D
    f     = -Inf * ones( sf + [ repmat( 2, [ 1 D ] ) ] );
else
    f     = -Inf * ones( sf + [ repmat( 2, [ 1 D ] ) 0 ] );
end

switch D
    case 1
        f( 2:end-1, : ) = f_tmp;
    case 2
        f( 2:end-1, 2:end-1, : ) = f_tmp;
    case 3
        f( 2:end-1, 2:end-1, 2:end-1, : ) = f_tmp;
end
clear f_tmp

%%%% Compute the EC curves
if strcmp( version, "C" )
    % C based implementation
    for n = 1:nEC
        ECn = EulerCharCrit_c( f( index{:}, n ), cc )';
        ECn = ECn( ECn( :, 2 ) ~= 0, : );
        ECn = ECn( ~isnan( ECn( :, 1 ) ), : );
        ECn = ECn( ECn( :, 1 )~=-Inf, : );

        [ ~, I ]   = sort( ECn( :, 1 ), 'ascend' );
        ECn     = ECn( I, : );
        EC{ n } = [ [ -Inf; ECn( :, 1 ); Inf ], [ 1; 1; 1 + ...
                    cumsum( ECn( :, 2 ) ) ] ];
    end
else
    % treat cases split by dimension
    switch D
        case 1
            % Compute signs of horizontal and vertical gradient
            df     = sign( f( 2:end, : ) - f( 1:end-1, : ) );
            minMax = df( 1:end-1, : ) - df( 2:end, : );
            f = f( 2:end-1, : );

            % Get locations of local minima and maxima and there numbe rin each
            % function
            mins = minMax == -2;
            maxs = minMax == 2;

            lmins = sum( mins );
            lmaxs = sum( maxs );

            % Compute the EC curves
            for k = 1:nEC
                dEC = [ [ f( mins ), ones( [ lmins( k ), 1 ] ) ];...
                        [ f( maxs ), -ones( [lmaxs( k ), 1 ] ) ] ];

                % sort crits according to its critical height
                [ ~, I ]    = sort( dEC( :, 1 ), 'ascend' );
                dEC         = dEC( I, : );
                dEC( 1, 2 ) = dEC( 1, 2 ) + L0;
                EC{ k }     = [ [ -Inf; dEC( :, 1 ); Inf ],...
                                [ 1; 1; cumsum( dEC( :, 2 ) ) ] ];
            end

        case 2

           for l = 1 : nEC
                %%%%%%% Find critical points
                CritMat = [ ];
                for i = 2:( sf( 1 ) + 1 )
                    for j = 2:( sf( 2 ) + 1 )
                            dEC = ECchange( f( i-1:i+1, j-1:j+1, l ), cc );
                            if dEC ~= 0
                                CritMat = [ CritMat; [ f(i, j, l) dEC ] ];
                            end
                    end
                end
                [ ~, I ]   = sort( CritMat( :, 1 ), 'ascend' );

                CritMat = CritMat( I, : );
                clear I
                EC{l}   = [ [ -Inf; CritMat( :, 1 ); Inf ], ...
                            [ 1; 1; 1 + cumsum( CritMat( :, 2 ) ) ] ];
           end
        case 3
            for l = 1 : nEC
                %%%%%%% Find critical points
                CritMat = [ ];
                for i = 2:( sf( 1 ) + 1 )
                    for j = 2:( sf( 2 ) + 1 )
                        for k = 2:( sf( 3 ) + 1 )
                            dEC = ECchange( f( i-1:i+1, j-1:j+1, k-1:k+1, l ), cc );
                            if dEC ~= 0
                                CritMat = [ CritMat; [ f(i, j, k, l) dEC ] ];
                            end
                        end
                    end
                end
                [ ~, I ]   = sort( CritMat( :, 1 ), 'ascend' );

                CritMat = CritMat( I, : );
                clear I
                EC{l}   = [ [ -Inf; CritMat( :, 1 ); Inf ], ...
                            [ 1; 1; 1 + cumsum( CritMat( :, 2 ) ) ] ];
            end
    end
end
%%%% End of the algorithm