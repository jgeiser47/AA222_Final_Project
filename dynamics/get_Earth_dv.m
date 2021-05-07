function Earth_dv = get_Earth_dv(C3)
% Purpose: Calculate Earth departure delta-V required for given Lambert arc
%          solution. Assumes no plane change and initial 400 km altitude
%          circular orbit. 

% Earth gravitational parameter
mu = 3.986e5; 

% Initial velocity for 400 km altitude circular orbit
r1 = 6378.137 + 400;
v1 = sqrt(mu/r1); 

% Periapsis velocity of hyperbolic trajectory 
energy = C3 / 2;
v2 = sqrt(2 * (energy + (mu/r1)));

% Calculate delta-V required
Earth_dv = abs(v2 - v1);
end