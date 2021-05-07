function oe = rv_to_oe(rv_ECI)
    % Author: Joshua Geiser
    % Inputs: Position vector in ECI (km), Velocity vector in ECI (km/s)
    % Outputs: Keplerian elements (km, --, rad, rad, rad, rad)

    % Assumes Sun orbit
    mu = 1.327122e11; 
    
    % Separate pos/vel vectors
    r_ECI = rv_ECI(1:3);
    v_ECI = rv_ECI(4:6);

    % Angular momentum vector and unit vector
    h_vec = cross(r_ECI,v_ECI);
    h_hat = h_vec / norm(h_vec);
    
    % SMA, mean motion, semi-latus rectum, and eccentricty
    a = 1 / ( (2/norm(r_ECI)) - (norm(v_ECI)^2/mu) );
    n = sqrt(mu / (a^3) );
    p = norm(h_vec)^2 / mu; 
    e = sqrt(1 - (p/a) );
    
    % Calc inclination based on angular momentum vector
    i = atan2(sqrt(h_hat(1)^2 + h_hat(2)^2), h_hat(3));
    
    % Calc RAAN based on angular momentum vector
    RAAN = atan2(h_hat(1), -h_hat(2));

    % Calculate Eccentric Anomaly, True Anomaly, and Mean Anomaly
    E = atan2(dot(r_ECI, v_ECI)/(a^2*n), (1-(norm(r_ECI)/a)));
    TA = 2*atan( tan(E/2) * sqrt( (1+e)/(1-e) ) );
    MA = E - e*sin(E);

    % Calcuate argument of latitude and AOP
    u = atan2(r_ECI(3)/sin(i), (r_ECI(1)*cos(RAAN) + r_ECI(2)*sin(RAAN)));
    AOP = u - TA;

    % Make sure to wrap angles between 0-2pi
    RAAN = wrapToPi(RAAN);
    AOP = wrapToPi(AOP);
    TA = wrapToPi(TA);
    MA = wrapToPi(MA);
    
    oe = [a;e;i;RAAN;AOP;MA];
end