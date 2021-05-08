function s = get_contour_data()
% Purpose: Get data for contour lines by running Lambert's algorithm across
%          search space. This function helps prevent running this over and
%          over again for each contour plot. 

% Output struct
s = struct(); 

% How fine the mesh grid is
N = 100;

% Get axis limits
[x_beg, x_end, y_beg, y_end] = get_datelims();

% Get arrays of MJDs
x = linspace(x_beg, x_end, N); 
y = linspace(y_beg, y_end, N); 

% Get mesh grid and setup Z (aka f(x)) variable
[X,Y] = meshgrid(x,y);
Z = zeros(size(X));

% Populate f(x) values given inputs x1 and x2
for i = 1:size(X,1)
    for j = 1:size(X,2)
        Z(i,j) = f(X(i,j), Y(i,j), 1);
    end
end

% Constraint lines
c1_x = linspace(x_beg, x_end+100, N);
c1_y = c1_x;
c2_y = c1_x + 365*1.5;

% Add values to struct
s.N = N;
s.x_beg = x_beg; s.x_end = x_end; 
s.y_beg = y_beg; s.y_end = y_end;
s.X = X; s.Y = Y; s.Z = Z;
s.c1_x = c1_x; s.c1_y = c1_y; s.c2_y = c2_y;
s.levels = [5,5.5,6,6.5,7,7.5,8,8.5,9,10,11,12,13,14,15];
end