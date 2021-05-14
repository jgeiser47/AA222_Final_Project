function [x_hist, samples_hist, elite_samples_hist] = cross_entropy(x0)
% Purpose: run Cross-Entropy algorithm 

% Hyperparameters
N_ITERS = 3;
N_SAMPLE = 500;
N_ELITE = 20;
dim = 2;

% Initial covariance size
cov_0 = 200^2 .* eye(2);

% Set initial conditions
mu = x0;
cov = cov_0;

% Prepare history arrays
x_hist = [x0];
samples_hist = [];
elite_samples_hist = [];

for i=1:N_ITERS
    
    % Get samples
    samples = sqrtm(cov) * randn(dim, N_SAMPLE) + mu;
    
    % Calculate objective function value for each sample
    f_samples = zeros(1,N_SAMPLE);
    for j = 1:N_SAMPLE
        f_samples(j) = f(samples(1,j), samples(2,j));
    end
    
    % Get indices of elite samples
    [~, indices_sorted] = sort(f_samples);
    indices_elite = indices_sorted(1:N_ELITE);
    
    % Fit distribution to elite samples
    elite_samples = samples(:, indices_elite);
    mu = mean(elite_samples, 2);
    cov = zeros(dim, dim);
    for j = 1:N_ELITE
        cov = cov + (elite_samples(:,j) - mu) * (elite_samples(:,j) - mu)';
    end
    cov = cov / N_ELITE;
    
    % Append to history arrays for plotting 
    x_hist = [x_hist, mu];
    samples_hist(:,:,i) = samples;
    elite_samples_hist(:,:,i) = elite_samples;
end
end