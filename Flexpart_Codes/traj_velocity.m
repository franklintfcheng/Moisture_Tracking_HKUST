%%
function [u, v, w] = traj_velocity(varargin)

% create a set of valid parameters and their default value
pnames = {'lon1', 'lat1', 'z1', 'lon2', 'lat2', 'z2', 'dt'};
dflts  = repmat([], 1, length(pnames));
         [  lon1,   lat1,   z1,   lon2,   lat2,   z2,  dt] = ...
                internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments

%/ NOTE: lon, lon, z can be arrays 
%/
%/       assume z1 and z2 are in m.
%/       assume dt is in hr.

%/ Compute rhumb line between two points (this assumes the propagation direction remains constant)
[arclen,az] = distance('rh', lat1,lon1,lat2,lon2); %/ rhumb line distance. arclen == theta between two points (in degree)
                                                   %/ IMPORTANT: This function comes from 'map_toolbox'!!

r = 6371;                           %/ Earth's radius in km
D = arclen/180*pi*r;                %/ distance over a sphere in km

%/ 90 - az will always obtain the correct angle  
speed = D/dt/3.6;                   %/ from km/h to m/s (assume dt is in hr)

%/ obtain the uv components
u = speed.*cosd(90 - az);            %/ in m/s. Here, a rhumb line distance is based on a constant direction.
v = speed.*sind(90 - az);  

%/ vertical velocity
w = (z2 - z1)/(dt*3600)*100;         %/ ** NOTE: in cm/s **

end