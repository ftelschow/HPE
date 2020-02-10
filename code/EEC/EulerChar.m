function EC = EulerChar(f, u, N, dir)

% Compute Euler characteristic of excursion set
% f is observed process
% u is vector of thresholds
% N is field dimension (1, 2, or 3)
% The first N dimensions represent the field
% EC is computed for all available fields of dimension N
% Uses 4-connectivity for 2D fields and 6-connectivity for 3D fields
% Uses continuouity in u to gain 20% computational efficiency

% Check input
if(N >= 4), error('N >= 4 not implemented'), end
sz = size(f);
if(length(sz) < N), error('N is larger than array size f'), end

% Size of EC array
if(length(sz) == N), sz = [length(u) 1];
else sz = [length(u) sz(N+1:end)];
end

if ~exist('dir','var')
    dir = 1;
end

if dir ~= 1
    u = -u;
end

EC = zeros(length(u)+1, prod(sz(2:end)));   % Extra element deleted later
A1 = false(size(f));
for j = length(u):-1:1
    EC(j,:) = EC(j+1,:);
    if dir == 1
        A = (f > u(j));
    else
         A = (f < u(j));
    end
    D = A & ~A1;    % set added
    switch N
        case 1
            I = find(any(D, 1));
            vertices = sum(A(:,I), 1);
            edges = sum(A(1:end-1,I) & A(2:end,I), 1);
            EC(j,I) = vertices - edges;
        case 2
            I = find(any(any(D, 1), 2));
            % 4-connectivity
            vertices = sum(sum(A(:,:,I), 1), 2);
            edges    = sum(sum(A(1:end-1,:,I) & A(2:end,:,I), 1), 2) + ...
                       sum(sum(A(:,1:end-1,I) & A(:,2:end,I), 1), 2);
            faces    = sum(sum(A(1:end-1,1:end-1,I) & A(2:end,1:end-1,I) & ...
                            A(1:end-1, 2:end, I) & A(2:end, 2:end, I), 1), 2);
            EC(j,I) = vertices - edges + faces;
        case 3
            I = find(any(any(any(D, 1), 2), 3));
            % 6-connectivity
            vertices = sum(sum(sum(A(:,:,:,I), 1), 2), 3);
            edges = sum(sum(sum(A(1:end-1,:,:,I) & A(2:end,:,:,I), 1), 2), 3) + ...
                    sum(sum(sum(A(:,1:end-1,:,I) & A(:,2:end,:,I), 1), 2), 3) + ...
                    sum(sum(sum(A(:,:,1:end-1,I) & A(:,:,2:end,I), 1), 2), 3);
            faces = sum(sum(sum(A(1:end-1,1:end-1,:,I) & A(2:end,1:end-1,:,I) & ...
                                A(1:end-1, 2:end, :,I) & A(2:end, 2:end, :,I), 1), 2), 3) + ...
                    sum(sum(sum(A(1:end-1,:,1:end-1,I) & A(2:end,:,1:end-1,I) & ...
                                A(1:end-1,:, 2:end, I) & A(2:end,:, 2:end ,I), 1), 2), 3) + ...
                    sum(sum(sum(A(:,1:end-1,1:end-1,I) & A(:,2:end,1:end-1,I) & ...
                                A(:,1:end-1, 2:end, I) & A(:,2:end, 2:end ,I), 1), 2), 3);
            cells = sum(sum(sum(A(1:end-1,1:end-1,1:end-1,I) & A(2:end,1:end-1,1:end-1,I) & ...
                                A(1:end-1, 2:end, 1:end-1,I) & A(2:end, 2:end, 1:end-1,I) & ...
                                A(1:end-1,1:end-1, 2:end, I) & A(2:end,1:end-1, 2:end, I) & ...
                                A(1:end-1, 2:end,  2:end, I) & A(2:end, 2:end,  2:end, I), 1), 2), 3);
            EC(j,I) = vertices - edges + faces - cells;
    end
    A1 = A; % new set becomes old
end
EC = reshape(EC(1:end-1,:), sz);
