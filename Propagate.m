% Joshua Geiser
% AA 222 Final Project
clear all; close all; clc;

cal_1 = [11, 05, 2013];
MJD_1 = cal_to_MJD(cal_1);
cal_2 = [09, 24, 2014];
MJD_2 = cal_to_MJD(cal_2);
data = solve_lambert(MJD_1, MJD_2);
plot_lambert_arc(data);

plot_contour()

function plot_contour()
% Purpose: Iterate over 2D array of departure dates and arrival dates to
%          create a porkchop plot. 

% 2013-2014
x = linspace(56550, 56800, 50);
y = linspace(56750, 57250, 50);

% 2020-2021 (Mars 2020: July 30, 2020 - Feb 18, 2021)
x = linspace(58950, 59200, 50);
y = linspace(59140, 59700, 50);

% Get mesh grid and setup Z (aka f(x)) variable
[X,Y] = meshgrid(x,y);
Z = zeros(size(X));

% Populate f(x) values given inputs x1 and x2
for i = 1:size(X,1)
    for j = 1:size(X,2)
        Z(i,j) = f(X(i,j), Y(i,j));
    end
end

% Create contour plot
figure(); hold on; grid on; axis equal;
levels = [7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39];
contour(X,Y,Z,levels);
colorbar;
xlabel('Earth Departure Date (MJD)');
ylabel('Mars Arrival Date (MJD)');
end


function fx = f(x1, x2)
% Purpose: Wrapper for optimization problem. Returns f(x) given inputs x1
%          and x2, where x1 is the Earth departure MJD and x2 is the Mars 
%          arrival MJD

data = solve_lambert(x1, x2);
fx = data.C3_1;
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
    
function rv_out = get_rv_Earth(s, MJD_f)
% Purpose: Get pos/vel of Earth in heliocentric frame at desired MJD

% Elapsed seconds
dt = (MJD_f-s.MJD_0) * 86400;

% Get mean motion 
a = s.kep_E0(1);
n = sqrt(s.mu / (a^3));

% Calculate final Mean Anomaly
M0 = s.kep_E0(6);
Mf = wrapToPi(M0 + (n*dt));

% Get pos/vel at final time
kep_Ef = s.kep_E0; 
kep_Ef(6) = Mf;
rv_out = oe_to_rv(kep_Ef);
end

function rv_out = get_rv_Mars(s, MJD_f)
% Purpose: Get pos/vel of Mars in heliocentric frame at desired MJD

% Elapsed seconds
dt = (MJD_f-s.MJD_0) * 86400;

% Get mean motion 
a = s.kep_M0(1);
n = sqrt(s.mu / (a^3));

% Calculate final Mean Anomaly
M0 = s.kep_M0(6);
Mf = wrapToPi(M0 + (n*dt));

% Get pos/vel at final time
kep_Mf = s.kep_M0; 
kep_Mf(6) = Mf;
rv_out = oe_to_rv(kep_Mf);
end

function s = get_ICs()
% Purpose: Get initial orbital elements and pos/vel of Earth and Mars at
%          J2000 epoch in heliocentric frame. These will be used to
%          propagate Earth/Mars to desired dates

% Structure for returning data
s = struct(); 

% Gravitational parameter of the Sun
s.mu = 1.327122e11;

% Initial time - 01/01/2000 (J2000 epoch)
s.cal_0 = [01, 01, 2000];
s.MJD_0 = cal_to_MJD(s.cal_0); 

% Earth elements at J2000 epoch
s.kep_E0 = [1; 0.01671123; -0.00001531; 100.46457166; 102.93768193; 0];
s.kep_E0 = jpl2kep(s.kep_E0);
s.rv_E0 = oe_to_rv(s.kep_E0);

% Mars elements at J2000 epoch
s.kep_M0 = [1.52371034; 0.09339410; 1.84969142; -4.55343205; -23.94362959; 49.55953891];
s.kep_M0 = jpl2kep(s.kep_M0);
s.rv_M0 = oe_to_rv(s.kep_M0);
end

function val_km = au2km(val_au)
% Purpose: Convert AU to km

val_km = val_au * 149597870.700; 
end

function val_kep = jpl2kep(val_jpl)
% Purpose: Given initial Earth and Mars orbital elements at J2000 epoch 
%          (from JPL website), convert to standard Keplerian elements of
%          the form [a, e, i, RAAN, AOP, M)

% Get desired elements
a = au2km(val_jpl(1));
e = val_jpl(2);
i = val_jpl(3);
RAAN = val_jpl(6); 
AOP = val_jpl(5) - RAAN;
M = val_jpl(4) - val_jpl(5);

% Put into array, convert angles to rad and wrap between -pi and pi
val_kep = [a; e; i; RAAN; AOP; M];
val_kep(3:6) = wrapToPi(deg2rad(val_kep(3:6)));
end