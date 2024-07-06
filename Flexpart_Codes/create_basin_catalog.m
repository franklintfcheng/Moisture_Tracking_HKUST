function [basin_catalog, str_domain_trajfile, str_domain] = create_basin_catalog(varargin)

    pnames = {'slct_reg', 'slct_var', 'X_tavg_filter', 'X_daily', 'lon', 'lat', 'date',...
              'data_folder', 'save_catalog', 'recompute_catalog'};

    dflts  = cell(1, length(pnames));

              [slct_reg,   slct_var,   X_tavg_filter,   X_daily,   lon,   lat,   date,...
               data_folder,   save_catalog,   recompute_catalog] = ...
                    internal.stats.parseArgs(pnames, dflts, varargin{:});
%%
    %======================================================================
    %/      Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: 28 Jun 2024
    %======================================================================
    if ischar(slct_reg)
        slct_reg = {slct_reg};
    end

    [~, bndry_data, name, id, ~] = reg_extractor('lon', lon, 'lat', lat, 'slct_reg', slct_reg,...
                                                 'data_folder', data_folder, 'savemat', 1, 'recompute', 0);
    
    catalog_name     = strjoin(slct_reg,'_');
    catalog_filename = strcat(data_folder, 'catalog_', catalog_name, '.mat');
    if isfile(catalog_filename) && recompute_catalog == 0
        fprintf('*** Loading: %s *** \n', catalog_filename)
        load(catalog_filename, 'basin_catalog');
    else
        lon(lon > 180) = lon(lon > 180) - 360;
        [lon_2D, lat_2D] = meshgrid(lon, lat);                                 %/ since Pm map does not contain two poles.
        lon_2D = lon_2D';  
        lat_2D = lat_2D';
        lon_2Dto1D = reshape(lon_2D, [], 1);                                   %/ be careful to check if the size is correct
        lat_2Dto1D = reshape(lat_2D, [], 1);                                   %/ be careful to check if the size is correct
        area        = calc_grid_area('lon', lon,  'lat', lat)*1e-6;            %/ m^2 to km^2
        area_2Dto1D = reshape(area,   [], 1);                                  %/ be careful to check if the size is correct
    
        if ~isempty(X_daily) && ~isempty(X_tavg_filter)
            X_tavg_filter_2Dto1D  = reshape(X_tavg_filter, [], 1);             %/ filtered 2D to 1D
        end
        
        basin_catalog = []; cnt = 0;
        for i = 1:length(bndry_data)
            if length(bndry_data{i}(:,1)) <= 3 
                continue;                                                      %/ skip it if the bndry data contains only 3 vertices or less!
            end 
    
            try
                [in, ~] = inpoly2([lon_2Dto1D, lat_2Dto1D], bndry_data{i});    %/ inpoly2 is 600xx faster than inpolygon!! 
            catch E
                warning(E.identifier, 'Error caught!: %s. Proceed with inpolygon().', E.message)
                [in, ~] = inpolygon(lon_2Dto1D, lat_2Dto1D, bndry_data{i}(:,1), bndry_data{i}(:,2));    %/ sometimes inpoly2 could run into a bug, then we use inpolygon.
            end
            ind_in = find(in == 1);
            
            hotspot_in_lon = lon_2Dto1D(ind_in);
            hotspot_in_lat = lat_2Dto1D(ind_in);
    
            ind_grids = find(ismember([lon_2Dto1D, lat_2Dto1D], [hotspot_in_lon, hotspot_in_lat], 'rows'));
    %         disp(length(ind_grids))
            
            if isempty(ind_grids)
                fprintf('skipped since no grid being assigned \n')
                continue;                                                      
            end 
    
            %/ Put cnt here after all if-statements involving continue
            cnt = cnt + 1;
            if ~isempty(id)        basin_catalog(cnt).id     = id(i);           else   basin_catalog(cnt).id     = '';       end
            if ~isempty(name)      basin_catalog(cnt).name   = name(i);         else   basin_catalog(cnt).name   = '';       end
            
            %/ Find the centroid of the hotspot region
            polyin = polyshape({bndry_data{i}(:,1)}, {bndry_data{i}(:,2)});
            [x,y] = polycenter(polyin);
    
            basin_catalog(cnt).area      = sum(area_2Dto1D(ind_grids));
            basin_catalog(cnt).bndry     = bndry_data{i};
            basin_catalog(cnt).centroid  = [x, y];
            
            %/ Optional
            if ~isempty(X_tavg_filter)
                basin_catalog(cnt).(strcat(slct_var, '_sum'))       = sum(X_tavg_filter_2Dto1D(ind_grids));
                basin_catalog(cnt).(strcat(slct_var, '_avg'))       = mean(X_tavg_filter_2Dto1D(ind_grids));
                basin_catalog(cnt).(strcat(slct_var, '_area_sum'))  = sum(X_tavg_filter_2Dto1D(ind_grids).*area_2Dto1D(ind_grids)); %/ consider area when comparing low-lat hotspots with high-lat ones.
                basin_catalog(cnt).(strcat(slct_var, '_area_wavg')) = sum(X_tavg_filter_2Dto1D(ind_grids).*area_2Dto1D(ind_grids))/sum(area_2Dto1D(ind_grids));
            end
            if ~isempty(date)   
                ndate = length(date);   
                basin_catalog(cnt).(strcat(slct_var, '_date')) = date;
            end
            
            %/ After finding hotspots, we restore unfiltered bwd contr. in the hotspots to compute its time series!!
            if ~isempty(X_daily)
                hotspot_var_daily_3Dto2D = reshape(X_daily, [], ndate);        %/ reshape to [grids time]
    
                basin_catalog(cnt).(strcat(slct_var, '_sum_ts')) = sum(hotspot_var_daily_3Dto2D(ind_grids, :), 1)';
                basin_catalog(cnt).(strcat(slct_var, '_avg_ts')) = mean(hotspot_var_daily_3Dto2D(ind_grids, :), 1)';
                basin_catalog(cnt).(strcat(slct_var, '_area_sum_ts'))  = sum(hotspot_var_daily_3Dto2D(ind_grids,:).*area_2Dto1D(ind_grids));
                basin_catalog(cnt).(strcat(slct_var, '_area_wavg_ts')) = sum(hotspot_var_daily_3Dto2D(ind_grids,:).*area_2Dto1D(ind_grids))/sum(area_2Dto1D(ind_grids));
            end
        end
        
        % %/ rank by strength. (??)
        % if ~isempty(slct_var) && ~isempty(criteria)
        %     T       = struct2table(basin_catalog);                             %/ convert to a table matrix first
        %     sortedT = sortrows(T, strcat(slct_var, '_', criteria), 'descend'); %/ sort by the requested criteria (i.e., BL_Pm_area_sum)
        %     basin_catalog = table2struct(sortedT);                             %/ change it back to struct array
        % end
    
        %/ Save the catalog
        if save_catalog
            fprintf('*** Saving basin_catalog into %s *** \n', catalog_filename)
            save(catalog_filename, 'basin_catalog', '-v7.3');
        end
    end

    %/ Some hard-coded strings just for my recent projects 
    if isequal(slct_reg, {'TP_basins'})
        str_domain_trajfile = sprintf('TP_basins%d', length(basin_catalog));
        str_domain      = sprintf('TP_basins%d', length(basin_catalog));

    elseif isequal(slct_reg, {'TP_Grass', 'TP_Semidesert', 'TP_Tundra'})
        str_domain_trajfile = 'TP_basins15';   %/ adopt the traj file saved in from_basin == 4 to save time and space!
        str_domain      = 'TP_main_LoVeg';

    elseif isequal(slct_reg, {'TP_grids_1x1'})
        str_domain_trajfile = 'TP_basins15';   %/ Since a grid box may cover areas outside TP boundary, 
                                               %/ we use the traj file saved in from_basin == 4 for consistency and convenience.
        str_domain      = catalog_name;

    else  %/ By default
        str_domain_trajfile = catalog_name;
        str_domain      = catalog_name;
    end

end