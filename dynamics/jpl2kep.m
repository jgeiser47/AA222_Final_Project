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