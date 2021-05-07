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