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