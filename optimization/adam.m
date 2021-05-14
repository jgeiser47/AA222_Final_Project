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