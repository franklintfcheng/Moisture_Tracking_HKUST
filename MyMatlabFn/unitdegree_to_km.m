function [grid_y, grid_x] = unitdegree_to_km(lon, lat)


    lon_range = [lon-0.5, lon+0.5];
    lat_range = [lat-0.5, lat+0.5];

    dx = abs(diff(lon_range));
    dy = abs(diff(lat_range));
    
    r_earth = 6.371e6;
    grid_y  = r_earth*dy/180*pi/1e5;         %/ in  100km
    R       = r_earth*cosd(mean(lat_range));
    grid_x  = R*dx/180*pi/1e5;               %/ in  100km
end