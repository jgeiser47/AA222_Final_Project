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