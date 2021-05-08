function contour_helper(s)
% Purpose: Helper function for formatting of the contour plots so that
%          these lines of code don't have to be repeated for each subplot

contour(toDateNum(s.X), toDateNum(s.Y), s.Z, s.levels);
patch(toDateNum([s.c1_x s.c1_x(end)]), toDateNum([s.c1_y s.c1_y(1)]), 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
patch(toDateNum([s.c1_x s.c1_x(1)]), toDateNum([s.c2_y s.c2_y(end)]), 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
xlabel('Earth Departure Date');
ylabel('Mars Arrival Date');
dateformat = 2;
datetick('x', dateformat);
datetick('y', dateformat);
xtickangle(45);
xlim([toDateNum(s.x_beg), toDateNum(s.x_end)]);
ylim([toDateNum(s.y_beg), toDateNum(s.y_end)]);
xticks(toDateNum([cal_to_MJD([1,  1, 2020]), ...
                  cal_to_MJD([4,  1, 2020]), ...
                  cal_to_MJD([7,  1, 2020]), ...
                  cal_to_MJD([10, 1, 2020]), ...
                  cal_to_MJD([1,  1, 2021])]));
end