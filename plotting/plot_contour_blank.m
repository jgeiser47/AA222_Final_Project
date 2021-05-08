function plot_contour_blank()
% Purpose: Iterate over 2D array of departure dates and arrival dates to
%          create a porkchop plot. 

% Get contour lines data by running Lambert's algorithm across search space
s = get_contour_data();

% Create contour plot
figure(); hold on; grid on; axis equal;
contour_helper(s);
title('\Delta V vs Departure/Arrival Date');
hcb = colorbar; 
hcb.Title.String = '\DeltaV';
hcb.Title.FontWeight = 'bold';
end