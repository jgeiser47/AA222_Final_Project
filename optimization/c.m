function cx = c(x1, x2)
% Purpose: Calculate constraint violations (constraints violated for values
%          greater than 0)

% Trajectory most have TOF between 0 and 1.5 years
cx = [x1 - x2; ...
      x2 - x1 - 365*1.5];
end