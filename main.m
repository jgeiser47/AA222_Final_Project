% Joshua Geiser
% AA 222 Final Project
clear all; close all; clc;

%% Plot of Lambert Arcs and Transfer Orbit Geometry
% cal_1 = [11, 05, 2020];
% MJD_1 = cal_to_MJD(cal_1);
% cal_2 = [03, 24, 2021];
% MJD_2 = cal_to_MJD(cal_2);
% data = solve_lambert(MJD_1, MJD_2);
% plot_lambert_arc(data);

%% Plot of Optimal Type I and Type II Transfers
% MJD_1 = 59056.573659816;
% MJD_2 = 59263.3863702012;
% data = solve_lambert(MJD_1, MJD_2);
% plot_lambert_arc(data);
% title('Type I Trajectory');
% 
% MJD_1 = 59041.0301348091;
% MJD_2 = 59365.698726598;
% data = solve_lambert(MJD_1, MJD_2);
% plot_lambert_arc(data);
% title('Type II Trajectory');

%% Blank Porchop Plot
%plot_contour_blank()

%% Comparison Plot of Gradient Descent and Adam
% % Get initial conditions
% rng(5);
% N = 3;
% x0s = get_x0s(N);
% 
% % Get N runs of gradient descent and adam for comparison
% x_hist_grad = [];
% x_hist_adam = [];
% for i = 1:N
%     x_hist_grad(:,:,i) = grad_descent(x0s(:,i));
%     x_hist_adam(:,:,i) = adam(x0s(:,i));
% end
% 
% % Plot comparison of first-order methods (gradient descent and adam)
% plot_contour_first_order_methods(x_hist_grad, x_hist_adam);

%% Cross-Entropy Method

% Run Cross-Entropy
rng(1);
N = 1;
x0 = get_x0s(N);
[x_hist, samples_hist, elite_samples_hist] = cross_entropy(x0);

% Plot contours for successive iterations of Cross-Entropy
plot_contour_cross_entropy(x_hist, samples_hist, elite_samples_hist);


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

function x_hist = adam(x0)
% Purpose: run Adam algorithm 

% Hyperparameters
alpha = 5;
gammav = 0.5;
gammas = 0.5;
eps = 1e-8; 
N_MAX = 1000;

% Initial conditions
x_hist = [x0];
x = x0;
k = 0; 
v_k = [0; 0];
s_k = [0; 0];

% Iterate until max iterations
while k < N_MAX
    
    % Update step
    g_k = g(x(1), x(2)); 
    v_k = gammav*v_k + (1-gammav)*g_k;
    s_k = gammas*s_k + (1-gammas)*(g_k.^2);
    k = k + 1; 
    v_hat = v_k ./ (1 - gammav^k);
    s_hat = s_k ./ (1 - gammas^k);
    x = x - ((alpha.*v_hat) ./ (eps + sqrt(s_hat)));
    
    % Append to x_hist array 
    x_hist = [x_hist, x];
end
end

function x_hist = grad_descent(x0)
% Purpose: run simple gradient descent algorithm (with normalized descent
%          direction and decaying learning rate)

% Hyperparameters and max iterations
alpha = 10;
gamma = 0.99;
N_MAX = 1000;

% Initial conditions
x_hist = [x0];
x = x0;

% Iterate until max iterations
i = 1;
while i <= N_MAX
    
    % Update step
    gx = g(x(1), x(2));
    x = x - alpha * (gx ./ norm(gx)); 
    alpha = alpha * gamma;
    
    % Append to x_hist array and increment iteration
    x_hist = [x_hist, x];
    i = i+1;
end
end

function x0s = get_x0s(N)
% Purpose: Get an array of N random initial conditions x0 for use in 
%          multiple optimization test cases

[x_beg, x_end, y_beg, y_end] = get_datelims();

x = (x_end - x_beg) .* rand(N,1) + x_beg;
y = (y_end - y_beg) .* rand(N,1) + y_beg;
x0s = [x, y]';
end

function plot_contour_first_order_methods(x_hist_grad, x_hist_adam)
% Purpose: Create porkchop plots comparing performance of the two first
%          order methods used: gradient descent and ADAM

% Get contour lines data by running Lambert's algorithm across search space
s = get_contour_data();

% Create contour plot
figure(); 

% Number of times each algorithm was run
N = size(x_hist_grad, 3);

% Gradient Descent subplot
subplot(1,2,1); hold on; grid on; axis equal;
contour_helper(s);
for i = 1:N
    plot(toDateNum(x_hist_grad(1,:,i)), toDateNum(x_hist_grad(2,:,i)), 'k');
end
title('Gradient Descent');

% Adam subplot
subplot(1,2,2); hold on; grid on; axis equal;
contour_helper(s);
for i = 1:N
    plot(toDateNum(x_hist_adam(1,:,i)), toDateNum(x_hist_adam(2,:,i)), 'k');
end
title('Adam');
hcb = colorbar; 
hcb.Title.String = '\DeltaV';
hcb.Title.FontWeight = 'bold';
end

function plot_contour_cross_entropy(x_hist, samples_hist, elite_samples_hist)
% Purpose: Create porkchop plot of cross-entropy iterations

% Get contour lines data by running Lambert's algorithm across search space
s = get_contour_data();

% Create contour plot
figure(); 

% Create a new subplot for each iteration
N_iters = 3;
for i = 1:N_iters
    subaxis(1,3,i, 'SpacingHorizontal', 0.1, 'SpacingVertical', 0.1); 
    hold on; grid on; axis equal;
    contour_helper(s);
    scatter(toDateNum(samples_hist(1,:,i)), toDateNum(samples_hist(2,:,i)), ...
            5, [0.5 0.5 0.5], 'filled');
    plot(toDateNum(x_hist(1,1:i+1)), toDateNum(x_hist(2,1:i+1)), 'k');
    scatter(toDateNum(elite_samples_hist(1,:,i)), toDateNum(elite_samples_hist(2,:,i)), ...
            10, 'r', 'filled');  
    scatter(toDateNum(x_hist(1,i+1)), toDateNum(x_hist(2,i+1)), ...
            100, 'g', 'filled', 'p');
    title(sprintf('Iteration %d', i));
end

% Add colorbar and adjust position
hcb = colorbar; 
hcb.Title.String = '\DeltaV';
hcb.Title.FontWeight = 'bold';
hcb.Position = [0.93 0.2 0.02 0.6];

% Adjust some other settings of printing/saving of figure as PDF
set(gcf, 'PaperPositionMode', 'Auto');
set(gcf, 'PaperOrientation', 'landscape');
end