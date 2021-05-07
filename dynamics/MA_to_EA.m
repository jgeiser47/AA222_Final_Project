function EA = MA_to_EA(M, e, tol)
    % Author: Joshua Geiser
    % Inputs: Mean anomaly (rad), eccentricity (--), tolerance (rad)
    % Output: Eccentric anomaly (rad)
    
    % Wrap M between 0 and 2 pi just to be safe
    M = wrapTo2Pi(M);
    
    % Safe initial guess
    E0 = pi;
    
    % Initialize some parameters
    EA = E0;
    delta_i = 1e10;
    
    % Newton-Raphson method
    while (abs(delta_i) > tol)
        delta_i = -1 * (EA - e*sin(EA) - M)/(1 - e*cos(EA));
        EA = EA + delta_i;
    end
end