function [pi_k, mu_k, Sigma_k, nu_k, responsibilities] = t_distribution_em(X, K, max_iter, tol)
% T_DISTRIBUTION_EM Performs EM algorithm for a mixture of t-distributions
%   X         : Data matrix (N x D), where N is the number of samples and D is the dimensionality
%   K         : Number of clusters/components
%   max_iter  : Maximum number of iterations
%   tol       : Tolerance for convergence
%
% Outputs:
%   pi_k           : Mixing coefficients (1 x K)
%   mu_k           : Means of the components (K x D)
%   Sigma_k        : Covariance matrices of the components (D x D x K)
%   nu_k           : Degrees of freedom for each component (1 x K)
%   responsibilities: Responsibility matrix (N x K)

% Get the size of the data
[N, D] = size(X);

% Initialize the parameters
[pi_k, mu_k, Sigma_k, nu_k] = initialize_parameters(X, K);

% Initialize variables
responsibilities = zeros(N, K);
log_likelihood = -inf;
converged = false;
iter = 0;

while ~converged && iter < max_iter
    iter = iter + 1;
    prev_log_likelihood = log_likelihood;
    
    % E-step: Compute responsibilities and auxiliary variables
    [responsibilities, w_ik, delta_ik] = expectation_step(X, pi_k, mu_k, Sigma_k, nu_k);
    
    % M-step: Update parameters
    [pi_k, mu_k, Sigma_k] = maximization_step(X, responsibilities, w_ik);
    
    % Update degrees of freedom (nu_k) using Newton-Raphson method
    nu_k = update_degrees_of_freedom(nu_k, responsibilities, w_ik, delta_ik, D);
    
    % Compute the log-likelihood
    log_likelihood = compute_log_likelihood(X, pi_k, mu_k, Sigma_k, nu_k);
    
    % Check for convergence
    if abs(log_likelihood - prev_log_likelihood) < tol
        converged = true;
    end
end

end

function [pi_k, mu_k, Sigma_k, nu_k] = initialize_parameters(X, K)
% Random initialization of parameters

[N, D] = size(X);

% Randomly initialize mixing coefficients
pi_k = ones(1, K) / K;

% Randomly initialize means by selecting random data points
rand_indices = randperm(N, K);
mu_k = X(rand_indices, :);

% Initialize covariance matrices as identity matrices
Sigma_k = repmat(eye(D), [1, 1, K]);

% Initialize degrees of freedom (can be set to a fixed value or randomly initialized)
nu_k = repmat(5, 1, K); % Starting with 5 degrees of freedom

end

function [responsibilities, w_ik, delta_ik] = expectation_step(X, pi_k, mu_k, Sigma_k, nu_k)
% E-step of the EM algorithm

[N, D] = size(X);
K = length(pi_k);

% Preallocate arrays
responsibilities = zeros(N, K);
w_ik = zeros(N, K);
delta_ik = zeros(N, K);

for k = 1:K
    % Compute delta_ik for each data point
    diff = X - mu_k(k, :);
    inv_Sigma = inv(Sigma_k(:, :, k));
    delta_ik(:, k) = sum((diff * inv_Sigma) .* diff, 2);
    
    % Compute the t-distribution pdf
    coef = exp(gammaln((nu_k(k) + D) / 2) - gammaln(nu_k(k) / 2));
    coef = coef / ((nu_k(k) * pi)^(D / 2) * sqrt(det(Sigma_k(:, :, k))));
    pdf = coef * (1 + delta_ik(:, k) / nu_k(k)).^(-(nu_k(k) + D) / 2);
    
    % Compute responsibilities (unnormalized)
    responsibilities(:, k) = pi_k(k) * pdf;
    
    % Compute auxiliary variable w_ik
    w_ik(:, k) = (nu_k(k) + D) ./ (nu_k(k) + delta_ik(:, k));
end

% Normalize responsibilities
sum_responsibilities = sum(responsibilities, 2);
responsibilities = responsibilities ./ sum_responsibilities;

end

function [pi_k, mu_k, Sigma_k] = maximization_step(X, responsibilities, w_ik)
% M-step of the EM algorithm

[N, D] = size(X);
K = size(responsibilities, 2);

% Update mixing coefficients
Nk = sum(responsibilities);
pi_k = Nk / N;

% Update means
mu_k = zeros(K, D);
for k = 1:K
    weights = responsibilities(:, k) .* w_ik(:, k);
    mu_k(k, :) = sum(weights .* X) / sum(weights);
end

% Update covariance matrices
Sigma_k = zeros(D, D, K);
for k = 1:K
    diff = X - mu_k(k, :);
    weights = responsibilities(:, k) .* w_ik(:, k);
    Sigma_k(:, :, k) = (diff' * (diff .* weights)) / sum(responsibilities(:, k));
end

end

function nu_k = update_degrees_of_freedom(nu_k, responsibilities, w_ik, delta_ik, D)
% Update degrees of freedom using Newton-Raphson method

K = length(nu_k);
epsilon = 1e-6; % Convergence threshold for Newton-Raphson
max_iter = 100; % Maximum iterations for Newton-Raphson

for k = 1:K
    nu = nu_k(k);
    for iter = 1:max_iter
        Nk = sum(responsibilities(:, k));
        w = w_ik(:, k);
        delta = delta_ik(:, k);
        
        % Compute the first derivative (dL/dnu)
        psi_nu_half = psi(nu / 2);
        psi_nu_D_half = psi((nu + D) / 2);
        term1 = 0.5 * Nk * (log(nu / 2) - psi_nu_half + 1);
        term2 = -0.5 * sum(responsibilities(:, k) .* (log((nu + delta) / 2) - psi_nu_D_half + w - 1));
        dL = term1 + term2;
        
        % Compute the second derivative (d^2L/dnu^2)
        trigamma_nu_half = psi(1, nu / 2);
        trigamma_nu_D_half = psi(1, (nu + D) / 2);
        term1 = 0.25 * Nk * (1 / (nu / 2) - trigamma_nu_half);
        term2 = 0.25 * sum(responsibilities(:, k) .* (1 / ((nu + delta) / 2) - trigamma_nu_D_half));
        d2L = term1 + term2;
        
        % Update nu using Newton-Raphson
        nu_new = nu - dL / d2L;
        
        % Ensure nu stays positive
        if nu_new <= 0
            nu_new = epsilon;
        end
        
        % Check for convergence
        if abs(nu_new - nu) < epsilon
            break;
        end
        
        nu = nu_new;
    end
    nu_k(k) = nu;
end

end

function log_likelihood = compute_log_likelihood(X, pi_k, mu_k, Sigma_k, nu_k)
% Compute the log-likelihood of the data given the current parameters

[N, D] = size(X);
K = length(pi_k);
log_likelihood = 0;

for n = 1:N
    prob = 0;
    for k = 1:K
        diff = X(n, :) - mu_k(k, :);
        delta = diff * (Sigma_k(:, :, k) \ diff');
        coef = exp(gammaln((nu_k(k) + D) / 2) - gammaln(nu_k(k) / 2));
        coef = coef / ((nu_k(k) * pi)^(D / 2) * sqrt(det(Sigma_k(:, :, k))));
        pdf = coef * (1 + delta / nu_k(k))^(-(nu_k(k) + D) / 2);
        prob = prob + pi_k(k) * pdf;
    end
    log_likelihood = log_likelihood + log(prob);
end

end
