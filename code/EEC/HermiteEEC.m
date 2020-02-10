function [LKC, EEC, rho] = HermiteEEC(EC, u, N, trueLKC, L0)
% Estimate LKC and EEC by the Hermite functional method
% [LKC, EEC] = HermiteEEC(EC, u, N, LKC, L0)
%
% Input:
%   EC      an array of observed EC curves of size length(u) x n x nsim
%   u       row vector of thresholds
%   N       dimension of domain (including scale)
%   trueLKC true LKC values (default = zeros(N, 1))
%   L0      EC of domain (default = 1)
%
% Output:
%   LKC     array with fields:
%               true        true LKC values
%               hat1        estimates from one curve at a time
%               hatn        estimate from all n curves
%               Sigma_hat   estimate of LKC covariance Sigma
%               se_hat      estimate of LKC estimation standard error
%               conf_hat    normal confidence intervals based on se_hat
%   EEC     array with fields:
%               true        true EEC values
%               hat1        estimates from one curve at a time
%               hatn        estimate from all n curves
%               se_hat      estimate of EEC estimation standard error
%               conf_hat    normal confidence intervals based on se_hat
%   rho     the EC densities

% Check input
sz = size(EC);
if(sz(1) ~= length(u))
    error('Dimensions of EC should be length(u) x n x nsim')
end
n = sz(2);
if(length(sz)>=3), nsim = sz(3); else nsim = 1; end
if ~exist('trueLKC'), trueLKC = zeros(N,1); end
if ~exist('L0'), L0 = 1; end

% Other parameters
du = u(2)-u(1);     % threshold increment


%% Hermite functional method

% Theoretical EEC
H = [ones(1,length(u)); u; (u.^2 - 1)]';    % Hermite polynomials
rho = H .* [exp(-u.^2/2)/(2*pi)^(2/2); exp(-u.^2/2)/(2*pi)^(3/2); ...
    exp(-u.^2/2)/(2*pi)^(4/2)]';            % EC densities
HH = H * diag([(2*pi)^(1/2); (2*pi)^(2/2); (2*pi)^(3/2)/2]);    % Operator
EEC = (1-normcdf(u))' * L0 + rho(:, 1:N) * trueLKC(1:N);

% Fit one curve at a time
L_hat = zeros(N, n, nsim);
EC_hat = zeros(length(u), n, nsim);
for j = 1:nsim
    for i = 1:n
        y = EC(:, i, j) - (1-normcdf(u))' * L0;
        L_hat(:, i, j) = HH(:, 1:N)' * y * du;
        EC_hat(:, i, j) = (1-normcdf(u))' * L0 + rho(:, 1:N) * L_hat(1:N, i, j);
    end
end

% Summary estimators: LKCs
L_hatn = permute(mean(L_hat, 2), [1 3 2]);
Sigma_hat = zeros(N, N, nsim);
L_se_hat = zeros(N, nsim);
for j = 1:nsim,
    Sigma_hat(:,:,j) = cov(L_hat(:,:,j)');
    L_se_hat(:, j) = sqrt(diag(Sigma_hat(:,:,j)) / n);
end
L_conf_hat = cat(3, L_hatn - 1.96*L_se_hat, L_hatn + 1.96*L_se_hat);

% Summary estimator: EEC
EEC_hatn = squeeze(mean(EC_hat, 2));
C = zeros(length(u), length(u), nsim);
EEC_se_hat = zeros(length(u), nsim);
for j = 1:nsim
    C(:,:,j) = rho(:, 1:N) * Sigma_hat(:,:,j) * rho(:, 1:N)';
    EEC_se_hat(:,j) = sqrt(diag(C(:,:,j)) / n);
end
EEC_conf_hat = cat(3, EEC_hatn - 1.96*EEC_se_hat, EEC_hatn + 1.96*EEC_se_hat);

% Summarize output
LKC  = struct('true', trueLKC, 'hat1', L_hat, 'hatn', L_hatn, 'Sigma_hat', Sigma_hat, ...
    'se_hat', L_se_hat, 'conf_hat', L_conf_hat);
EEC = struct('true', EEC, 'hat1', EC_hat, 'hatn', EEC_hatn, ...
    'cov', C, 'se_hat', EEC_se_hat, 'conf_hat', EEC_conf_hat);

return