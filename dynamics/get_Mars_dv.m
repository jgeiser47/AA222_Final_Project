function Mars_dv = get_Mars_dv(C3)
% Purpose: Calculate Mars arrival delta-V required for given Lambert arc
%          solution. Assumes no plane change and final 400 km altitude
%          circular orbit. 

% Mars gravitational parameter
mu = 42828.3719; 

% Initial velocity for 400 km altitude circular orbit
r1 = 3396.19 + 400;
v1 = sqrt(mu/r1); 

% Periapsis velocity of hyperbolic trajectory 
energy = C3 / 2;
v2 = sqrt(2 * (energy + (mu/r1)));

% Calculate delta-V required
Mars_dv = abs(v2 - v1);
end