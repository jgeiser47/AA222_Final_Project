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