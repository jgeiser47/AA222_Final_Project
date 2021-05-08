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