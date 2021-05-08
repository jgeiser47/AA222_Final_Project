function gx = g(x1, x2)
% Purpose: Get gradient of function using central difference method

% Step size
h = 0.1;

% Partial derivatives wrt x1 and x2
fprime1 = (f(x1+(h/2),x2) - f(x1-(h/2),x2)) / h;
fprime2 = (f(x1,x2+(h/2)) - f(x1,x2-(h/2))) / h;

% Return gradient
gx = [fprime1; fprime2];
end