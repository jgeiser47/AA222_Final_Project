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
plot(0,0, 'y*');
plot(data.r_Earth(1), data.r_Earth(2), 'bo', 'MarkerFaceColor', 'b');
plot(data.r_Mars(1), data.r_Mars(2), 'ro', 'MarkerFaceColor', 'r');
plot(r_transfer(1,:), r_transfer(2,:), 'm--');
plot(r_Earth_arr(1,:), r_Earth_arr(2,:), 'b--');
plot(r_Mars_arr(1,:), r_Mars_arr(2,:), 'r--');
xlabel('X (km)'); ylabel('Y (km)');
legend('Sun', 'Earth', 'Mars', 'Transfer');
end