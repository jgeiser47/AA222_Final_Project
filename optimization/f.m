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