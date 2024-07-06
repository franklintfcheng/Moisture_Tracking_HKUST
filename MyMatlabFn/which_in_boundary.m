%%
function [map_data_in, area_in, ind_in] = which_in_boundary(varargin)
    
    % create a set of valid parameters and their default value
    pnames = {'bndry_data_inpoly2', 'map_data', 'lon', 'lat', 'output_1D_or_2D'};
    dflts  = cell(length(pnames), 1);
    [bndry_data_inpoly2, map_data, lon, lat, output_1D_or_2D] ...
                        = internal.stats.parseArgs(pnames, dflts, varargin{:});

    %/ Check dims
    if (size(map_data,1) ~= length(lon)) || (size(map_data,2) ~= length(lat))
        error('Dimension of the input map_data does not fit lon & lat!')
    end
    ndate = size(map_data, 3);  %/ ndate == 1 if input map_data is 2D.
    
    
    %/ [IMPORTANT]: Always convert from [0 360) to (-180, 180]
    lon(lon > 180)             = lon(lon > 180) - 360;
    bndry_lon                  = bndry_data_inpoly2(:,1);
    bndry_lon(bndry_lon > 180) = bndry_lon(bndry_lon > 180) - 360;
    bndry_data_inpoly2(:, 1)   = bndry_lon;
    
    
    [lon_2D, lat_2D] = meshgrid(lon, lat);
    lon_2D = lon_2D';  lat_2D = lat_2D';
    lon_2Dto1D  = reshape(lon_2D, [], 1);
    lat_2Dto1D  = reshape(lat_2D, [], 1);
    area        = calc_grid_area('lon', lon,  'lat', lat);     %/ m^2
    

    [in,~] = inpoly2([lon_2Dto1D, lat_2Dto1D], bndry_data_inpoly2);        %/ inpoly2 is 600xx faster than inpolygon!!
    ind_in = find(in == 1);

    map_data_2Dto1D = reshape(map_data, [], ndate);
    area_2Dto1D     = reshape(area, [], 1);                                %/ be careful to check if the size is correct

    
    if isequal(output_1D_or_2D, '1D')
        map_data_in = map_data_2Dto1D(ind_in,:);
        area_in     = area_2Dto1D(ind_in);
        
    elseif isequal(output_1D_or_2D, '2D')
        map_data_in                  = nan(size(map_data));
        map_data_in_2Dto1D           = reshape(map_data_in, [], ndate);                                   
        map_data_in_2Dto1D(ind_in,:) = map_data_2Dto1D(ind_in,:);      
        map_data_in                  = reshape(map_data_in_2Dto1D, length(lon), length(lat), ndate);
        
        area_in         = nan(length(lon), length(lat));
        area_in         = reshape(area_in, [], 1);                                    
        area_in(ind_in) = area_2Dto1D(ind_in); 
        
        area_in         = reshape(area_in, length(lon), length(lat)); %/ reshape back to 2D
        
        
        
    else
        error('wrong input of ''output_1D_or_2D''!');
    end
end