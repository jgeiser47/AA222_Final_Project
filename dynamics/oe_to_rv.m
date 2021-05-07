function rv_ECI = oe_to_rv(oe)
    % Author: Joshua Geiser
    % Inputs: Keplerian elements (km, --, rad, rad, rad, rad)
    % Outputs: 3x2 matrix containing pos/vel vectors in inertial frame

    % Get orbital elements
    a = oe(1);
    e = oe(2);
    i = oe(3); 
    RAAN = oe(4); 
    AOP = oe(5); 
    MA = oe(6);
    
    % Assumes Sun orbit
    mu = 1.327122e11; 

    % Calculate a couple necessary parameters
    tol = 10e-10;
    E = MA_to_EA(MA, e, tol);
    n = sqrt(mu / a^3);
    
    % Rotation matrix from perifocal to inertial (ECI) coordinates
    R_PQW_to_IJK = PQW_to_IJK(i, RAAN, AOP);
    
    % Position vector
    r_PQW = [a*(cos(E)-e); a*sqrt(1-e^2)*sin(E); 0];
    r_IJK = R_PQW_to_IJK * r_PQW;
    
    % Velocity vector
    v_PQW = ((a*n)/(1-e*cos(E))) .* [-sin(E); sqrt(1-e^2)*cos(E); 0];
    v_IJK = R_PQW_to_IJK * v_PQW; 
    
    rv_ECI = [r_IJK; v_IJK];
end