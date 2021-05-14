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