
function grid_area = calc_grid_area_header(header)

% Dec 4 2020  Fandy

error('This function has been deprecated! use ''calc_grid_area'' instead!!');

r_earth=6.371e6;

%/ Since the meridional length of the grid does not vary with latitude
grid_y = r_earth*header.dyout/180*pi;
grid_y_array = repmat(grid_y, 360/header.dxout, 1); %/ column vector

grid_center_lat = header.latp;
R = r_earth*cosd(grid_center_lat);
R = R(1:end-1);                                     %/ discard the first and last element as it is related to the north pole. No area for points at 90N or 90S.
grid_x_array = R*header.dxout/180*pi;               %/ row vector

grid_area = grid_y_array*grid_x_array;


end
