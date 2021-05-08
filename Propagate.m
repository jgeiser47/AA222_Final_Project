% Joshua Geiser
% AA 222 Final Project
clear all; close all; clc;

% cal_1 = [11, 05, 2020];
% MJD_1 = cal_to_MJD(cal_1);
% cal_2 = [03, 24, 2021];
% MJD_2 = cal_to_MJD(cal_2);
% data = solve_lambert(MJD_1, MJD_2);
% plot_lambert_arc(data);

%x_hist = grad_descent();
x_hist = adam();

plot_contour_blank()
plot_contour_first_order_methods(x_hist)

%x0s = get_x0s(3)

function x_hist = adam()
% Purpose: run ADAM algorithm 

% Hyperparameters
alpha = 5;
gammav = 0.5;
gammas = 0.5;
eps = 1e-8; 
N_MAX = 100;

% Initial condition
x0 = [58900; 59600]; 
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

function x_hist = grad_descent()
% Purpose: run simple gradient descent algorithm

% Hyperparameters and max iterations
alpha = 10;
gamma = 0.99;
N_MAX = 1000;

% Initial condition
x0 = [59100; 59800]; %[59000; 59100]; [59200; 59800];
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
x0s = [x, y];
end

function plot_contour_first_order_methods(x_hist)
% Purpose: Create porkchop plots comparing performance of the two first
%          order methods used: gradient descent and ADAM

% Get contour lines data by running Lambert's algorithm across search space
s = get_contour_data();

% Create contour plot
figure(); 

subplot(1,2,1); hold on; grid on; axis equal;
contour_helper(s);
title('Gradient Descent');

subplot(1,2,2); hold on; grid on; axis equal;
contour_helper(s);
title('Adam');
colorbar;

if (x_hist)
    plot(toDateNum(x_hist(1,:)), toDateNum(x_hist(2,:)), 'k')
end
end