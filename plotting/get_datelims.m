function [x_beg, x_end, y_beg, y_end] = get_datelims()
% Purpose: Get beginning/end MJDs for X and Y axes of search space

% 2020-2021 launch window (Mars 2020 dates: July 30, 2020 - Feb 18, 2021)
x_beg = cal_to_MJD([01, 01, 2020]);
x_end = cal_to_MJD([01, 05, 2021]); 
y_beg = cal_to_MJD([08, 01, 2020]);
y_end = cal_to_MJD([09, 01, 2022]);
end