function x0s = get_x0s(N)
% Purpose: Get an array of N random initial conditions x0 for use in 
%          multiple optimization test cases

[x_beg, x_end, y_beg, y_end] = get_datelims();

x = (x_end - x_beg) .* rand(N,1) + x_beg;
y = (y_end - y_beg) .* rand(N,1) + y_beg;
x0s = [x, y]';
end