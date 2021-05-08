function val_out = toDateNum(MJD)
% Purpose: Helper function for converting MJD to datenum for plotting
val_out = datenum(datetime(MJD, 'convertfrom', 'modifiedjuliandate'));
end