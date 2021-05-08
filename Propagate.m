% Joshua Geiser
% AA 222 Final Project
clear all; close all; clc;

cal_1 = [11, 05, 2020];
MJD_1 = cal_to_MJD(cal_1);
cal_2 = [03, 24, 2021];
MJD_2 = cal_to_MJD(cal_2);
data = solve_lambert(MJD_1, MJD_2);
%plot_lambert_arc(data);

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

function [x_beg, x_end, y_beg, y_end] = get_datelims()
% Purpose: Get beginning/end MJDs for X and Y axes of search space

% 2020-2021 launch window (Mars 2020 dates: July 30, 2020 - Feb 18, 2021)
x_beg = cal_to_MJD([01, 01, 2020]);
x_end = cal_to_MJD([01, 05, 2021]); 
y_beg = cal_to_MJD([08, 01, 2020]);
y_end = cal_to_MJD([09, 01, 2022]);
end

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

function contour_helper(s)
% Purpose: Helper function for formatting of the contour plots so that
%          these lines of code don't have to be repeated for each subplot

contour(toDateNum(s.X), toDateNum(s.Y), s.Z, s.levels);
patch(toDateNum([s.c1_x s.c1_x(end)]), toDateNum([s.c1_y s.c1_y(1)]), 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
patch(toDateNum([s.c1_x s.c1_x(1)]), toDateNum([s.c2_y s.c2_y(end)]), 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
xlabel('Earth Departure Date');
ylabel('Mars Arrival Date');
dateformat = 2;
datetick('x', dateformat);
datetick('y', dateformat);
xtickangle(45);
xlim([toDateNum(s.x_beg), toDateNum(s.x_end)]);
ylim([toDateNum(s.y_beg), toDateNum(s.y_end)]);
xticks(toDateNum([cal_to_MJD([1,  1, 2020]), ...
                  cal_to_MJD([4,  1, 2020]), ...
                  cal_to_MJD([7,  1, 2020]), ...
                  cal_to_MJD([10, 1, 2020]), ...
                  cal_to_MJD([1,  1, 2021])]));
end

function plot_contour_blank()
% Purpose: Iterate over 2D array of departure dates and arrival dates to
%          create a porkchop plot. 

% Get contour lines data by running Lambert's algorithm across search space
s = get_contour_data();

% Create contour plot
figure(); hold on; grid on; axis equal;
contour_helper(s);
colorbar;
title('\Delta V vs Departure/Arrival Date');

end

function plot_contour_first_order_methods(x_hist)
% Purpose: Iterate over 2D array of departure dates and arrival dates to
%          create a porkchop plot. 

% Get contour lines data by running Lambert's algorithm across search space
s = get_contour_data();

% Create contour plot
figure(); 

subplot(1,2,1); hold on; grid on; axis equal;
contour_helper(s);

subplot(1,2,2); hold on; grid on; axis equal;
contour_helper(s);
colorbar;

if (x_hist)
    plot(toDateNum(x_hist(1,:)), toDateNum(x_hist(2,:)), 'k')
end
end

function val_out = toDateNum(MJD)
% Purpose: Helper function for converting MJD to datenum for plotting
val_out = datenum(datetime(MJD, 'convertfrom', 'modifiedjuliandate'));
end

function fx = f(x1, x2, no_penalty)
% Purpose: Wrapper for optimization problem. Returns f(x) given inputs x1
%          and x2, where x1 is the Earth departure MJD and x2 is the Mars 
%          arrival MJD

% Flag for adding penalties or not
if ~exist('no_penalty', 'var')
    no_penalty = 0;
end

% Calculate total dV first without penalties
data = solve_lambert(x1, x2);
fx = data.dv_tot;

% Add penalties to f(x) (if running an algorithm, not for plotting)
if ~no_penalty
    rho_count = 5;
    rho_quadr = 0.05;
    cx = c(x1, x2);
    fx = fx + (rho_count * sum(max(cx,0)>0)) + (rho_quadr * sum(max(cx,0).^2));
end
end

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

function cx = c(x1, x2)
% Purpose: Calculate constraint violations (constraints violated for values
%          greater than 0)

% Trajectory most have TOF between 0 and 1.5 years
cx = [x1 - x2; ...
      x2 - x1 - 365*1.5];
end

function data = solve_lambert(MJD_1, MJD_2)
% Purpose: Main function for setting up dynamics and solving Lambert's
%          problem for two input MJDs. Returns 'data' struct containing
%          relevant values for this solution

% Initial orbital elements of Earth and Mars
s = get_ICs();

% Get pos/vel of Earth at desired Earth departure epoch
rv_Earth = get_rv_Earth(s, MJD_1);
r_Earth = rv_Earth(1:3);
v_Earth = rv_Earth(4:6);

% Get pos/vel of Mars at desired Mars arrival epoch
rv_Mars = get_rv_Mars(s, MJD_2);
r_Mars = rv_Mars(1:3);
v_Mars = rv_Mars(4:6);

% Time of flight (seconds)
dt = (MJD_2 - MJD_1) * 86400;

% Solve Lambert's problem for both Type I and Type II trajectories
[v1_out_l,v2_out_l,errorl] = AA279lambert_vallado_u(s.mu,r_Earth,r_Mars,'l',0,dt);
[v1_out_s,v2_out_s,errorl] = AA279lambert_vallado_u(s.mu,r_Earth,r_Mars,'s',0,dt);

% Calculate departure C3 for both solution
C3_l = norm(v_Earth - v1_out_l)^2;
C3_s = norm(v_Earth - v1_out_s)^2;

% Determine which solution is better
if C3_l < C3_s
    v1_out = v1_out_l;
    v2_out = v2_out_l;
else
    v1_out = v1_out_s;
    v2_out = v2_out_s;
end

% Format output variables into 'data' structure
data = struct();
data.MJD_1 = MJD_1;
data.MJD_2 = MJD_2;
data.dt = dt;
data.r_Earth = r_Earth; 
data.v_Earth = v_Earth;
data.r_Mars = r_Mars;
data.v_Mars = v_Mars;
data.v1_out = v1_out;
data.v2_out = v2_out;
data.C3_1 = norm(v_Earth - v1_out)^2;
data.C3_2 = norm(v_Mars - v2_out)^2;
data.C3_tot = data.C3_1 + data.C3_2;
data.dv_1 = get_Earth_dv(data.C3_1);
data.dv_2 = get_Mars_dv(data.C3_2);
data.dv_tot = data.dv_1 + data.dv_2;
end

function plot_lambert_arc(data)
% Purpose: Plot a single lambert arc solution in the XY plane (since
%          solutions are primarily in the XY plane

% Initial pos/vel of Earth and Mars
s = get_ICs();

% Get full orbits of Earth and Mars in r_Earth_arr and r_Mars_arr
N = 100;
M_arr = linspace(-pi,pi,N);
r_Earth_arr = zeros(3,N);
r_Mars_arr  = zeros(3,N);
for i = 1:N
    M = M_arr(i);
    
    % Earth position
    kep_Earth = s.kep_E0;
    kep_Earth(6) = M;
    rv_Earth = oe_to_rv(kep_Earth);
    r_Earth_arr(:,i) = rv_Earth(1:3);
    
    % Mars position
    kep_Mars = s.kep_M0;
    kep_Mars(6) = M;
    rv_Mars = oe_to_rv(kep_Mars);
    r_Mars_arr(:,i) = rv_Mars(1:3);
end

% Get start and end mean anomalies of transfer arc
kep_transfer0 = rv_to_oe([data.r_Earth; data.v1_out]);
n_transfer = sqrt(s.mu / kep_transfer0(1)^3);
M0 = kep_transfer0(6);
Mf = M0 + n_transfer*(data.dt);

% Get full transfer arc history in r_transfer 
M_arr = linspace(M0, Mf, N);
r_transfer = zeros(3,N);
for i = 1:N
    M = M_arr(i);
    kep_transfer = kep_transfer0;
    kep_transfer(6) = M;
    rv_transfer = oe_to_rv(kep_transfer); 
    r_transfer(:,i) = rv_transfer(1:3);
end

% Generate trajectory plot
figure(); hold on; grid on; axis equal;
plot(r_Earth_arr(1,:), r_Earth_arr(2,:), 'b--');
plot(r_Mars_arr(1,:), r_Mars_arr(2,:), 'r--');
plot(r_transfer(1,:), r_transfer(2,:), 'm--');
plot(0,0, 'y*');
plot(data.r_Earth(1), data.r_Earth(2), 'bo', 'MarkerFaceColor', 'b');
plot(data.r_Mars(1), data.r_Mars(2), 'ro', 'MarkerFaceColor', 'r');
xlabel('X (km)'); ylabel('Y (km)');
end
