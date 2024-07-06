%%
function [bndry_data, hotspot] = saveload_hs(varargin)

%/ create a set of valid parameters and their default value
pnames = {'recompute', 'hs_data_filename', 'top', 'WSV', 'slct_var', 'prct', 'criteria', 'topx', 'include_oceans'};  
dflts  = cell(1, length(pnames));

[recompute, hs_data_filename, top, WSV, slct_var, prct, criteria, topx, include_oceans] ...
               = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments

% flag_exist = 0;
% %/ check if data exists
% if ismember(save_or_load, {'load'})
%     
%     if isfile(hs_data_filename)
%         fprintf('*** The queried data are found. *** \n*** Loading: %s *** \n', hs_data_filename)
%         load(hs_data_filename);
%         flag_exist = 1;
%     else
%         fprintf('*** The queried data are NOT found. *** \n *** Now trying to derive it... ***\n');
%         
%     end
% end

% if ismember(save_or_load, {'save'}) || flag_exist == 0
        
%/ check if data exists
if isfile(hs_data_filename) && recompute == 0
    fprintf('*** The queried data are found. *** \n*** Loading: %s *** \n', hs_data_filename)
    load(hs_data_filename);   
else
    lon = WSV.lon; 
    lat = WSV.lat; 
    nlon = length(lon);
    nlat = length(lat); 
    lon(lon > 180) = lon(lon > 180) - 360;                             %/ for inpolygon later.

    [lon_2D, lat_2D] = meshgrid(lon, lat);
    lon_2D = lon_2D';  lat_2D = lat_2D';            %/ lon x lat
    lon_2Dto1D       = reshape(lon_2D, [], 1);
    lat_2Dto1D       = reshape(lat_2D, [], 1);

    ind_land = which_on_land('pos_lon_array', lon_2Dto1D, 'pos_lat_array', lat_2Dto1D, 'pos_hgt_array',[], 'hgt_range', []);
    ind_ocean = setdiff([1:length(lon_2Dto1D)]', ind_land);

    reg_2D_ocean           = zeros(nlon*nlat, 1) - 999;                    %/ -999 --> oceans
    reg_2D_ocean(ind_land) = nan;  
    reg_2D_ocean           = reshape(reg_2D_ocean, nlon, nlat);

    reg_2D_land            = nan(nlon*nlat, 1); 
    reg_2D_land(ind_land)  = 999;                                          %/ +999 --> land
    reg_2D_land            = reshape(reg_2D_land, nlon, nlat);

    %---- Prepare 95-p filtered data -----%
    X_daily = WSV.(strcat(slct_var, '_daily'));                        %/ for deriving the daily contribution of the hotspots
    X_tavg = mean(X_daily, 3);    
    if ~isempty(prct)                                                  %/ filtered hotspot grids for later identifying hotspots!!!
        prctile_val = prctile(reshape(X_tavg, 1,[]), prct); 
        X_tavg_filter = X_tavg;
        X_tavg_filter(X_tavg_filter <= prctile_val) = 0;                           
    end

    %----------------- Hotspots by ocean basins --------------------------
    if include_oceans 
%         slct_oceans = {'AS', 'BoB', 'NTAO', 'SWTAO', 'GoG', 'GoM'};   %/ these are hotspots that I observed.
%         slct_oceans_longname = {'Arabian Sea', 'Bay of Bengal', 'N. Tropical Atlantic', 'SW. Tropical Atlantic', 'Gulf of Guinea', 'Gulf of Mexico'};   %/ these are hotspots that I observed.

        slct_oceans = {'GoG', 'NTAO', 'STAO', 'WSAO', 'ESAO', 'SO', 'GoM', 'CaribeanSea', 'WNAO', 'ENAO', 'HudsonBay', 'AO',...
                       'MeditSea', 'BlackSea', 'CaspianSea', 'AS', 'BoB', 'RedSea', 'PersianOmanGulf', 'WIO', 'EIO', 'SCS', 'GoT',...
                       'WTPO', 'WSPO', 'YSECS', 'SoJ', 'WNPO', 'ENPO', 'ETPO', 'ESPO'};
        slct_oceans_longname = slct_oceans;

        id = 1:length(slct_oceans);
        name = slct_oceans_longname;
        bndry_data = {};
        for o = 1:length(slct_oceans)

            reg_2D_ocean_bc = reg_2D_ocean;

%             %/ separating land from oceans
%             hotspot_var_2Dto1D           = reshape(hotspot_var', [], 1);
%             hotspot_var_2Dto1D(ind_land) = nan;
%             hotspot_var_ocean            = reshape(hotspot_var_2Dto1D, size(hotspot_var', 1), size(hotspot_var', 2));
%             hotspot_var_land             = hotspot_var' - hotspot_var_ocean;

            box_reg = box_region(slct_oceans(o));

            mask_2D = inpolygon(lon_2D, lat_2D, box_reg(:,1), box_reg(:,2));
            reg_2D_ocean_bc(~mask_2D) = nan;                               %/ Since the boundary is not perfect, use land/ocean mask data to obtain the correct grids

            %/ point data for plotting (if needed)
    %                 [reg_lat_ind, reg_lon_ind] = find(~isnan(hotspot_var_ocean));
    %                 point_data = [point_data; [lon_2D(1,reg_lon_ind)', lat_2D(reg_lat_ind, 1)]];

            %/ final logical matrix for outlining boudaries
            logical_2D = reg_2D_ocean_bc;
            logical_2D(~isnan(logical_2D)) = 1;   %/ first set nonnan to be 1 (as values can be zeros even in the target region)
            logical_2D(isnan(logical_2D))  = 0;   %/ then set nan to be 0

            %/ outline ocean boundaries (%/ glb_data_mode == 1 to draw boundaries across meridian)
            bndry_cell = get_bndry_from_logi('logical_2D', logical_2D, 'bndry_lon', lon, 'bndry_lat', lat,...
                                             'glb_data_mode', 1, 'outputcell', 0, 'draw_rings', 0, 'hole_mode', 'noholes');  
            bndry_data(o,1) = {bndry_cell};
        end

        %---- Ocean and land masks on the 95-p filtered data -----%
        %/ lon x lat (95p-filtered, over ocean)
        X_tavg_filter_ocean           = reshape(X_tavg_filter, [], 1);
        X_tavg_filter_ocean(ind_land) = 0;                                 %/ set it to zeros 
        X_tavg_filter_ocean           = reshape(X_tavg_filter_ocean, nlon, nlat);

        %/ lon x lat x time (over ocean)
        X_daily_ocean           = reshape(X_daily, nlon*nlat, []);
        X_daily_ocean(ind_land) = 0;                                       %/ set it to zeros 
        X_daily_ocean           = reshape(X_daily_ocean, nlon, nlat, []);

        ocean_hotspot = create_hs_catalog('bndry_data', bndry_data, 'slct_var', slct_var,...
                                          'X_tavg_filter', X_tavg_filter_ocean, 'X_daily', X_daily_ocean,...
                                          'lon',      lon, 'lat',       lat, 'date',  WSV.slct_dates, 'contin', [], 'id', id, 'name', name,...
                                          'criteria',  criteria);
    end

    %----------------- Hotspots by river basins --------------------------
    hydrosheds_folder = '/disk/r059/tfchengac/hydrosheds/';
    contin_id = {'Africa', 'Europe', 'Siberia', 'Asia', 'Australia', 'South America', 'North America', 'Arctic', 'Greenland'}; %/ do not change this ordering.
    file_code = {'af', 'ar', 'as', 'au' ,'eu', 'gr', 'na', 'sa', 'si', 'na_mex'}'; %/ arbitrary ordering, 'na_mex' is made by me.
    
    bndry_data = {}; point_data = []; cnt = 0; text_data = {};
    for k = 1:length(file_code)
        if isequal(file_code{k}, 'na_mex') 
            level = 4;
            basin_file = sprintf('hybas_%s_lev%02d_v1c.shp', 'na', level);
        else
            level = 3;
            basin_file = sprintf('hybas_%s_lev%02d_v1c.shp', file_code{k}, level);
        end
        fprintf('*** Reading watershed shapefile: %s *** \n', basin_file)

        [S,A] = shaperead(strcat(hydrosheds_folder, basin_file),'UseGeoCoords',true);

        basin_name = hydroshed_name('noOfbasin', length(S), 'contin', file_code{k}, 'level', level);   %/ match the river basin name (can modify)

        for i = 1:length(S)
            %/ skip the unwanted basins
            if (i == 1 || i == 23)     && isequal(file_code{k}, 'na')      continue; end   %/ skip the level3 mexico.
            if ~ismember(i, [1, 2, 5]) && isequal(file_code{k}, 'na_mex')  continue; end   %/ skip those not from mexico.
            if ismember(i, [11:13])    && isequal(file_code{k}, 'sa')      continue; end   %/ skip those not from mexico.

            %/ separating land from oceans
%                 hotspot_var_2Dto1D           = reshape(X_tavg', [], 1);
%                 hotspot_var_2Dto1D(ind_land) = nan;
%                 hotspot_var_ocean            = reshape(hotspot_var_2Dto1D, size(X_tavg', 1), size(X_tavg', 2));
% 
%                 hotspot_var_land = X_tavg';
%                 hotspot_var_land(~isnan(hotspot_var_ocean))  = nan;

            if i == 10 && isequal(file_code{k}, 'sa')              
                box_reg = box_region({'NESA'});           flag = 1;            %/ Replace those river basins with the customized one
                basin_name{i} = 'NE. South America';                           %/ update basin name

            elseif i == 1 && isequal(file_code{k}, 'na_mex')
                box_reg = box_region({'SMexico'}); 	      flag = 1;
                basin_name{i} = 'S. Mexico';                                   %/ update basin name

            elseif i == 2 && isequal(file_code{k}, 'na_mex')
                box_reg = box_region({'CMexico'});        flag = 1;
                basin_name{i} = 'C. Mexico';                                   %/ update basin name

            elseif i == 5 && isequal(file_code{k}, 'na_mex')
                box_reg = box_region({'WNMexico'}); 	  flag = 1;
                basin_name{i} = 'WN. Mexico';                                  %/ update basin name

            elseif i == 17 && isequal(file_code{k}, 'au')
                %/ original one is too heavy)
                box_reg = box_region({'NewGuinea'}); 	  flag = 0;            %/ flag = 0: do not auto-outline.

            elseif i == 10 && isequal(file_code{k}, 'as')
                %/ original one is too heavy)
                box_reg = box_region({'Pearl'}); 	      flag = 0;            %/ flag = 0: do not auto-outline.
            else
                box_reg = [S(i).Lon; S(i).Lat]';          flag = 0;
            end

            reg_2D_land_bc = reg_2D_land;

            %/ obtain the points from bndry in land-only hotspot var.
            mask_2D = inpolygon(lon_2D, lat_2D, box_reg(:,1), box_reg(:,2));
            reg_2D_land_bc(~mask_2D) = nan;                                  %/ set grids outside the baisns to be nan

%                 [reg_lat_ind, reg_lon_ind] = find(~isnan(reg_2D_land_bc));
            % point_data = [point_data; [lon_2D(1,reg_lon_ind)', lat_2D(reg_lat_ind, 1)]];

            %/ create the final logical matrix for outlining boudaries
            logical_2D = reg_2D_land_bc;
            logical_2D(~isnan(logical_2D)) = 1;                                %/ first set nonnan to be 1 (as values can be zeros even in the target region)
            logical_2D(isnan(logical_2D))  = 0;                                %/ then set nan to be 0

            %/ outline boundaries
            cnt = cnt + 1;
            if flag == 1                                                       %/ use the outlined boundary 
                bndry_data{cnt,1} = get_bndry_from_logi('logical_2D', logical_2D, 'bndry_lon', lon, 'bndry_lat', lat,...
                                                        'glb_data_mode', 1, 'outputcell', 0, 'draw_rings', 0, 'hole_mode', 'noholes');
            else
                bndry_data(cnt,1) = {box_reg};                                 %/ use the river basin's boundary
            end
            id(cnt,1)         = A(i).HYBAS_ID;
            contin{cnt,1}     = contin_id{floor(A(i).HYBAS_ID/1e9)};
            text_data{cnt,1} = mean([S(i).BoundingBox(:,1)], 'omitnan');
            text_data{cnt,2} = mean([S(i).BoundingBox(:,2)], 'omitnan');
    %         text_data{cnt,3} = num2str(i);                %/ for checking: to locate the index with basins first.
            text_data(cnt,3) = basin_name(i);

        end
    end

    %/ lon x lat (95p-filtered, over land)
    X_tavg_filter_land            = reshape(X_tavg_filter, [], 1);
    X_tavg_filter_land(ind_ocean) = 0;                                 %/ set it to zeros 
    X_tavg_filter_land            = reshape(X_tavg_filter_land, nlon, nlat);

    %/ lon x lat x time (over land)
    X_daily_land            = reshape(X_daily, nlon*nlat, []);
    X_daily_land(ind_ocean) = 0;                                       %/ set it to zeros 
    X_daily_land            = reshape(X_daily_land, nlon, nlat, []);

    basin_hotspot = create_basin_catalog('bndry_data',    bndry_data, 'slct_var',   slct_var, ...
                                      'X_tavg_filter', X_tavg_filter_land, 'X_daily', X_daily_land, ...
                                      'lon',  lon,   'lat',   lat,   'date',  WSV.slct_dates,...
                                      'contin',   contin,       'id',         id, 'name',  text_data(:,3), 'criteria',  criteria);       

    if include_oceans
        %/ combine oceanic and basin hotspots
        hotspot = [ocean_hotspot; basin_hotspot];
    else
        hotspot = basin_hotspot;
    end

    %/ rank by strength.
    T       = struct2table(hotspot);                                       %/ convert to a table matrix first
    sortedT = sortrows(T, strcat(slct_var, '_', criteria), 'descend');     %/ sort by the field value
    hotspot = table2struct(sortedT);                                       %/ change it back to struct array

    %/ Simplify the boundary data!!
    for i = 1:length(hotspot)
        a = hotspot(i).bndry;
        ind_nan = find(isnan(a(:,1))); %/ since nans are mean to cut the boundary, we need to keep those nans when simplifying the boundary!
        ind_nan = [0; ind_nan; length(a)];
        b = [];
        for j = 1:length(ind_nan)-1
            a_subset = a(ind_nan(j)+1:ind_nan(j+1),:);
            if isempty(a_subset)  continue;   end

            intvl = fix(length(a_subset)*0.01);  %/ simplied to a size of just ~100
            if intvl ~= 0
                b = [b; a_subset(1,:); a_subset(2:intvl:end-1,:); a_subset(end,:)];
            else
                b = [b; a_subset];  %/ otherwise just copy it.
            end
        end
        hotspot(i).bndry_simp = b;
    end

    %/ Save the land hotspots
    fprintf('*** Saving hotspots data: %s *** \n', hs_data_filename)
    save(hs_data_filename, 'hotspot', '-v7.3');
end
% end

%/ output the simplified boundary data
if isempty(top)
    bndry_data = {hotspot(1:topx).bndry_simp}';
else
    bndry_data = {hotspot(top).bndry_simp}';
end
    
end