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