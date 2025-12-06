function samples = quasi_mvnrnd(mu, Sigma, n)
%QUASI_MVNRND Quasi-random sampling from a multivariate normal distribution
%   samples = quasi_mvnrnd(mu, Sigma, n) returns n samples from a
%   multivariate normal distribution with mean mu and covariance Sigma,
%   using Sobol quasi-random sequences.

    % Validate inputs
    d = length(mu);
    if ~isequal(size(Sigma), [d, d])
        error('Covariance matrix Sigma must be of size length(mu) x length(mu)');
    end

    % Create Sobol sequence generator
    p = sobolset(d, 'Skip', 1e3, 'Leap', 1e2); 
    p = scramble(p, 'MatousekAffineOwen'); % scrambling for better uniformity

    % Generate quasi-random samples from unit hypercube
    u = net(p, n); % n-by-d matrix with values in (0,1)

    % Convert uniform samples to standard normal using inverse CDF
    z = norminv(u, 0, 1); % n-by-d matrix of standard normal samples

    % Transform to multivariate normal using Cholesky decomposition
    L = chol(Sigma, 'lower');
    samples = bsxfun(@plus, z * L', mu(:)');
end
