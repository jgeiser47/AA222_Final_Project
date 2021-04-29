function MJD = cal_to_MJD(UT1_cal, varargin)
    % Author: Joshua Geiser
    % Inputs: Required: UT1 calendar date in array form [MM,DD,YYYY]
    %         Optional: Time of calendar day in array form [HH,MM,SS]
    % Outputs: UT1 epoch expressed in Modified Julian Days (MJD)
    
    % Set month, day, and year variables
    MM = UT1_cal(1);
    DD = UT1_cal(2);
    YYYY = UT1_cal(3); 
    
    % Year runs from March 1 to end of February in computing leap years
    if MM <= 2
        y = YYYY-1;
        m = MM + 12;
    else
        y = YYYY;
        m = MM;
    end
    
    % Auxillary quantity to account for leap days
    B = (y/400) - (y/100) + (y/4);
    
    % Calculate MJD according to Eq (A.6) in Montenbruck
    MJD = 365*y - 679004 + floor(B) + floor(30.6001 * (m+1)) + DD;
    
    % Calculate fraction of day (if HH:MM:SS input) and add to MJD
    if nargin == 2
        HHMMSS_arr = varargin{1};
        HH = HHMMSS_arr(1);
        MM = HHMMSS_arr(2);
        SS = HHMMSS_arr(3);
        
        frac_MJD = (HH + (MM/60) + (SS/3600)) / 24;
        MJD = MJD + frac_MJD;
    end
end