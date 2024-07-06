%%
function T2mSST = T2mSST(varargin)

pnames = {'Res', 'lon', 'lat', 'slct_field', 'T2m', 'SST'};

dflts  = {1, [] [] [] [] []};

[Res, lon, lat, slct_field, T2m, SST] = internal.stats.parseArgs(pnames, dflts, varargin{:});

stlon = min(lon);
edlon = max(lon);
stlat = min(lat);
edlat = max(lat);

path(path, 'landmask');

[lon_2D, lat_2D] = meshgrid(stlon:Res:edlon,edlat:-Res:stlat);
[bigy, bigx] = size(lat_2D);

% [lat_2D_wrap,slon_2D_wrap] = meshgrid(stlat:Res:edlat, ...
%                                      [stlon:Res:(180-Res),-180:Res:(edlon-360)]);
% [lon_2D_wrap,lat_2D_wrap] = meshgrid([stlon:Res:(180-Res),-179:Res:(edlon-360+1)],...
%                                        edlat:-Res:stlat);
x = stlon:Res:edlon;
xind = find(x > 179);
x(xind) = x(xind) - 360;

% [lon_2D_wrap,lat_2D_wrap] = meshgrid(x,  linspace(edlat,stlat, 106));
[lon_2D_wrap,lat_2D_wrap] = meshgrid(x,  edlat:-Res:stlat);

ocean_ind = ~landmask(lat_2D_wrap, lon_2D_wrap, 100); % 100% quality
ocean_ind = ocean_ind'; % transpose

% slct_field = {'daily'};
% slct_field = {'daily_clim'};
% fld = fieldnames(SST);
% for i = 1:length(fld)
%     if ~ismember(fld(i), slct_field) == 1
%         T2mSST.(fld{i}) = SST.(fld{i});
%     end
% end

disp(size(SST))
T2mSST = nan(size(SST));
disp(size(T2mSST))

 
% disp(size(ocean_ind));
% disp(ocean_ind)
% disp(bigy)
% filling oceanic grids with SST
for j = 1:bigy
%     disp(j)
%     disp(size(ocean_ind(:,j)))

    T2mSST(ocean_ind(:,j), j, :) = SST(ocean_ind(:,j), j,:); 
end

% filling land grids with T2m
[ind_land_lon, ind_land_lat] = find(isnan(T2mSST(:,:,1)));

for i = 1:length(ind_land_lon)
    T2mSST(ind_land_lon(i),ind_land_lat(i),:)...
                = T2m(ind_land_lon(i),ind_land_lat(i),:);
end

%/ convert to degree C if it is not an anomaly field
if ~ismember(slct_field, {'dailyAnom_PMC', 'monthlyAnom'})
    T2mSST = T2mSST - 273.15;
end

end