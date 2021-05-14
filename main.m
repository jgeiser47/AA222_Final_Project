% Joshua Geiser
% AA 222 Final Project
clear all; close all; clc;

%% Plot of Aribtrary Lambert Arc
cal_1 = [11, 05, 2020];
MJD_1 = cal_to_MJD(cal_1);
cal_2 = [03, 24, 2021];
MJD_2 = cal_to_MJD(cal_2);
data = solve_lambert(MJD_1, MJD_2);
plot_lambert_arc(data);

%% Plot of Optimal Type I and Type II Transfers
MJD_1 = 59056.573659816;
MJD_2 = 59263.3863702012;
data = solve_lambert(MJD_1, MJD_2);
plot_lambert_arc(data);
title('Type I Trajectory');

MJD_1 = 59041.0301348091;
MJD_2 = 59365.698726598;
data = solve_lambert(MJD_1, MJD_2);
plot_lambert_arc(data);
title('Type II Trajectory');

%% Blank Porchop Plot
plot_contour_blank()

%% Comparison Plot of Gradient Descent and Adam
% Get initial conditions
rng(5);
N = 3;
x0s = get_x0s(N);

% Get N runs of gradient descent and adam for comparison
x_hist_grad = [];
x_hist_adam = [];
for i = 1:N
    x_hist_grad(:,:,i) = grad_descent(x0s(:,i));
    x_hist_adam(:,:,i) = adam(x0s(:,i));
end

% Plot comparison of first-order methods (gradient descent and adam)
plot_contour_first_order_methods(x_hist_grad, x_hist_adam);

%% Cross-Entropy Method
% Get initial conditions
rng(1);
N = 1;
x0 = get_x0s(N);

% Run 1 run of Cross-Entropy algorithm
[x_hist, samples_hist, elite_samples_hist] = cross_entropy(x0);

% Plot contours for successive iterations of Cross-Entropy
plot_contour_cross_entropy(x_hist, samples_hist, elite_samples_hist);
