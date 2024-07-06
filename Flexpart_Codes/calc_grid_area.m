function grid_area = calc_grid_area(varargin)

    pnames = {'lon', 'lat', 'debug_mode'};

    dflts = cell(1,length(pnames));

    [         lon,    lat,   debug_mode] = internal.stats.parseArgs(pnames, dflts, varargin{:});
    %%
    %==================================================================================================
    %/       Author: Franklin Cheng (fandycheng@ust.hk)
    %/  Last Update: Mar 21, 2024
    %/
    %/       Output: area (m^2)
    %==================================================================================================
    % 
    % %/ Determine grid sizes dlat and dlon: 
    % [dlat1,dlat2] = gradient(lat); 
    % [dlon1,dlon2] = gradient(lon); 
    % 
    % % We don't know if lat and lon were created by 
    % % [lat,lon] = meshgrid(lat_array,lon_array) 
    % % or [lon,lat] = meshgrid(lon_array,lat_array),
    % % but we can find out: 
    % 
    % if isequal(dlat1,zeros(size(lat)))
    %    dlat = dlat2; 
    %    dlon = dlon1; 
    %    assert(isequal(dlon2,zeros(size(lon)))==1,'Error: lat and lon must be monotonic grids, as if created by meshgrid.') 
    % else
    %    dlat = dlat1; 
    %    dlon = dlon2; 
    %    assert(isequal(dlon1,dlat2,zeros(size(lon)))==1,'Error: lat and lon must be monotonic grids, as if created by meshgrid.') 
    % end
    % 
    % dy = dlat.*R*pi/180;
    % dx = (dlon/180).*pi.*R.*cosd(lat); 
    % 
    % A = double(abs(dx.*dy)); 

    r_earth = 6371e3;
    if isvector(lon) && isvector(lat)
        if size(lat, 2) == 1         %/ column to row vector
            lat = lat';         
        end   

        dx = abs(diff(lon(1:2)));
        dy = abs(diff(lat(1:2)));

        %/ Since the meridional length of the grid (almost) does not vary with latitude
        grid_y = r_earth*dy/180*pi;
        grid_y_array = repmat(grid_y, length(lon), 1); %/ column vector, round() is to make the number integer.

        R = r_earth*cosd(lat);
        grid_x_array = R.*dx/180*pi;               %/ row vector
        grid_area    = grid_y_array*grid_x_array;
    else
        %/=================================================================
        %/ NOTE: Since the input lon, lat are likely non-monotonic matrices
        %/       it is almost impossible to shift the grids such that we
        %/       could compute the area centered at (i,j).
        %/       Hence, area of the last lon and lat will be zero.
        %/=================================================================

        % Convert lon and lat to radians
        lon_rad = deg2rad(lon);
        lat_rad = deg2rad(lat);

        % Initialize grid area matrix
        grid_area = zeros(size(lon));

        % Compute grid cell area
        for i = 1:size(lon, 1)-1
            for j = 1:size(lon, 2)-1
        % for i = 97
            % for j = 291
                % Define polygon vertices
                lon1 = lon_rad(i, j);
                lat1 = lat_rad(i, j);
                lon2 = lon_rad(i+1, j);
                lat2 = lat_rad(i+1, j);
                lon3 = lon_rad(i+1, j+1);
                lat3 = lat_rad(i+1, j+1);
                lon4 = lon_rad(i, j+1);
                lat4 = lat_rad(i, j+1);
                
                % fprintf('lon1 = %.3f%s, radian = %.3f\n', lon(i,  j),   char(176), lon1);
                % fprintf('lon2 = %.3f%s, radian = %.3f\n', lon(i+1,j),   char(176), lon2);
                % fprintf('lon3 = %.3f%s, radian = %.3f\n', lon(i+1,j+1), char(176), lon3);
                % fprintf('lon4 = %.3f%s, radian = %.3f\n', lon(i,  j+1), char(176), lon4);

                %/ Handle grid cells that cross the prime meridian
                thres = 1.5*pi;
                if (lon1 - lon2) > thres
                    lon2 = lon2 + 2*pi;
                elseif (lon1 - lon2) < -thres
                    lon2 = lon2 - 2*pi;
                end
                
                if (lon1 - lon3) > thres
                    lon3 = lon3 + 2*pi;
                elseif (lon1 - lon3) < -thres
                    lon3 = lon3 - 2*pi;
                end

                if (lon1 - lon4) > thres
                    lon4 = lon4 + 2*pi;
                elseif (lon1 - lon4) < -thres
                    lon4 = lon4 - 2*pi;
                end

                lon_poly = [lon1 lon2 lon3 lon4 lon1];
                lat_poly = [lat1 lat2 lat3 lat4 lat1];

                % Compute polygon area using polyarea function
                R = r_earth * cosd(mean([lat1, lat2, lat3, lat4]));
                area = polyarea(lon_poly, lat_poly) * R^2;

                % Assign area to grid cell
                grid_area(i, j) = area;
            end
        end

        %/==== OLD CODE, fail to handle non-monotonic grid cells ===
        %/ A bit complicated to compute the area of uneven grid cells
        % dx      = diff(lon,[],1);
        % dx_half = reshape([dx'/2;dx'/2], size(lon,2), [])';
        % dx_half = [dx_half(1,:); dx_half; dx_half(end,:)];
        % dx_adj  = reshape(sum(reshape(dx_half,2,[]),1), size(lon,1), []);
        % R       = r_earth*cosd(lat);
        % grid_x  = R.*dx_adj/180*pi;               %/ row vector
        % 
        % dy      = diff(lat,[],2);
        % dy_half = reshape([dy/2; dy/2], size(lat,1), []);
        % dy_half = [dy_half(:,1), dy_half, dy_half(:,end)];
        % dy_adj  = reshape(sum(reshape(dy_half',2,[]),1), size(lat,1), []);
        % 
        % grid_y    = r_earth*dy_adj/180*pi;
        % grid_area = grid_x.*grid_y;
    end

    if debug_mode   %/ IMPORTANT!
        close all
        figure
        imagesc(grid_area)
        colorbar
        min(grid_area, [], 'all', 'omitnan')
        max(grid_area, [], 'all', 'omitnan')
        [ind_max_row, ind_max_col] = find(grid_area == max(grid_area, [], 'all', 'omitnan'));
        fprintf('ind_max_row(1) = %d,  ind_max_col(1) = %d\n', ind_max_row(1), ind_max_col(1));

        if ~isempty(find(grid_area < 0, 1))
            [ind_nve_row, ind_nve_col] = find(grid_area < 0);
            fprintf('ind_nve_row(1) = %d,  ind_nve_col(1) = %d\n', ind_nve_row(1), ind_nve_col(1));
            error('grid_area contains negative values! Check your code or input!')
        end
    end

end
