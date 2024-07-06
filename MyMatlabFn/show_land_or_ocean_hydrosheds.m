%%
function [output, cond_land, cond_ocean] = show_land_or_ocean_hydrosheds(varargin)
    
    pnames = {'matrix', 'pos_lon_array', 'pos_lat_array', 'lon_grids', 'lat_grids',...
              'land_or_ocean', 'plateau_hgt', 'output_as_logical', 'cond_land', 'cond_ocean'};
    dflts  = cell(1, length(pnames));
    
    [          matrix,   pos_lon_array,   pos_lat_array,   lon_grids,   lat_grids,...
               land_or_ocean,   plateau_hgt,   output_as_logical,   cond_land,   cond_ocean] = internal.stats.parseArgs(pnames, dflts, varargin{:});
    
    %/========================================================================================
    %/ Author: Fandy CHENG
    %/ Date of last update: 26 Jan 2023
    %/
    %/ This function is considered *more robust and intuitive* than 'which_on_land' function.
    %/========================================================================================
    %/ INPUT:       
    %/              'matrix': a 2D or 3D matrix
    %/       'pos_lon_array': a 1D array of longitudinal positions
    %/       'pos_lat_array': a 1D array of latitudnial  positions
    %/           'cond_land': a preprocessed cond_land -> save time and Map_Toolbox license! (Strongly recommended)
    %/          'cond_ocean': a preprocessed cond_land -> save time and Map_Toolbox license! (Strongly recommended)
    %/
    %/ OUTPUT: 
    %/              'output': data with land or oceanic values shown.
    %/           'cond_land': a logical matrix for land  (based on avg_topo_map from 5-min topography)
    %/          'cond_ocean': a logical matrix for ocean (based on avg_topo_map from 5-min topography)
    %/
    %/========================================================================================
%%
    %=============== Create random data for testing ===============
%     warning('''show_land_or_ocean'' in testing mode!!');
%     rng(0,'twister');
%     test_lon_range = [-180 180];
%     a = test_lon_range(1); b = test_lon_range(2);
%     pos_lon_array = (b-a).*rand(50000,1) + a;
%     min(pos_lon_array)
%     max(pos_lon_array)
    
%     test_lat_range = [-90 90];
%     a = test_lat_range(1); b = test_lat_range(2);
%     pos_lat_array = (b-a).*rand(50000,1) + a;
%     min(pos_lat_array)
%     max(pos_lat_array)
    %=================================================================
%     tic;
    if ~isempty(matrix) && (~isempty(pos_lon_array) || ~isempty(pos_lat_array))
        error('Either input a ''matrix'' or a pair of ''pos_lon_array'' and ''pos_lat_array''!');
    end
    
    if isempty(matrix) && isempty(pos_lon_array) && isempty(pos_lat_array)
        %/ output will be empty.
        %/ But will output cond_land & cond_ocean based on lon_grids & lat_grids.
        output = [];
    end
    
    if isempty(plateau_hgt)  
        plateau_hgt = 3000;  
    end

    %/========== Obtain Topo Map (now only for capturing  Antarctica that is absent from HydroSHEDS database) ===========%
    %/ Based on the input lon/lat grids averaged from 5-min topography.
    if isvector(lon_grids) && isvector(lat_grids)
        flag_2D_lonlat = 0;
    else
        flag_2D_lonlat = 1;
    end
    lon_grids(lon_grids < 0) = lon_grids(lon_grids < 0) + 360;  %/ always make lon_grids in [0, 360).

%     tic
    %/========== Obtain cond_land and cond_ocean (now based on HydroSHEDS) ===========%
    if isempty(cond_land) || isempty(cond_ocean)
        %/ Get the averaged topo map (for capturing Antarctica)
        avg_topo_map = interp_topo('lon_grids', lon_grids, 'lat_grids', lat_grids);
        topo_gt0 = (avg_topo_map > 0);
        
        %/ Get a list of basin boundaries by contenation
        labels_hydrosheds = {'Africa', 'Europe', 'Siberia', 'Asia', 'Australia', 'South America', 'North America', 'Arctic', 'Greenland'};
        slct_reg_bc = [labels_hydrosheds, 'Antarctica'];
        reg_bndry_list = cell(length(slct_reg_bc), 1); %/ pre-allocation
        reg_name_list  = cell(length(slct_reg_bc), 1); %/ pre-allocation
        for k = 1:length(slct_reg_bc)
            if any(ismember(slct_reg_bc{k}, labels_hydrosheds))
                %/ Load Level 3 river basin from hydrosheds
                [bndry_data, ~, basin_name_list, ~] = retrieve_hydrosheds('slct_conti', slct_reg_bc{k}, 'simplify_bndry', 0);
                reg_name_list{k} = basin_name_list;
                if length(bndry_data) ~= length(basin_name_list)
                    error('inconsistent length!')
                end
            elseif isequal(slct_reg_bc{k}, 'Antarctica')
                %/ Load my prescribed ploygon regions
                bndry_data     = box_region(slct_reg_bc(k));
                reg_name_list{k} = slct_reg_bc{k};
            else
                error('check you code!')
            end
            reg_bndry_list{k} = bndry_data;
        end
        reg_bndry_list = cat(1, reg_bndry_list{:}); %/ concatenate the nested cells
        reg_name_list  = cat(1, reg_name_list{:});  %/ concatenate the nested cells
%         length(reg_bndry_list)
%         length(reg_name_list)
        
        if flag_2D_lonlat
            lon_bc_2D = lon_grids;
            lat_bc_2D = lat_grids;
            lon_bc_2D(lon_bc_2D > 180) = lon_bc_2D(lon_bc_2D > 180) - 360; %/ as HydroSHEDS boundaries are from (-180 to 180].
            [nlon, nlat] = size(lon_bc_2D);
        else
            lon_bc = lon_grids;
            lat_bc = lat_grids;
            lon_bc(lon_bc > 180) = lon_bc(lon_bc > 180) - 360; %/ as HydroSHEDS boundaries are from (-180 to 180].
            [lon_bc_2D, lat_bc_2D] = meshgrid(lon_bc, lat_bc);
            lon_bc_2D     = lon_bc_2D';  
            lat_bc_2D     = lat_bc_2D';
            nlon = length(lon_grids);
            nlat = length(lat_grids);
        end
        lon_bc_2Dto1D = reshape(lon_bc_2D, [], 1);
        lat_bc_2Dto1D = reshape(lat_bc_2D, [], 1);
        cond_land     = nan(nlon, nlat);
        for i = 1:length(reg_bndry_list)
%             fprintf('Region name: %s\n', reg_name_list{i})
            bndry_data = reg_bndry_list{i}; 
            try
                [mask_2D, ~] = inpoly2([lon_bc_2Dto1D, lat_bc_2Dto1D], bndry_data);                          %/ inpoly2 is 600xx faster than inpolygon!! 
                mask_2D = reshape(mask_2D, nlon, nlat);
            catch
%             catch E
%                 warning(E.identifier, 'Error caught!: %s. Proceed with inpolygon().', E.message)
                [mask_2D, ~] = inpolygon(lon_bc_2D, lat_bc_2D, bndry_data(:,1), bndry_data(:,2));    %/ sometimes inpoly2 could run into a bug, then we use inpolygon.
            end
                                
            if isequal(reg_name_list{i}, 'Antarctica')
                cond_land(mask_2D & topo_gt0) = 1;  %/ since Antarctica boundaries are only sketchy outlines. It needs to combine with topo > 0.
            else
                cond_land(mask_2D) = 1;
            end
        end
        cond_land(isnan(cond_land)) = 0; %/ finally, set the reamining nan as 0.
        cond_land = logical(cond_land);
        cond_ocean = ~cond_land;
    end
    
    if land_or_ocean == 3
        avg_topo_map = interp_topo('lon_grids', lon_grids, 'lat_grids', lat_grids);
        
        if ~isempty(find(isnan(avg_topo_map), 1))  
            error('avg_topo contains NaN!');   
        end
        cond_plateau = (avg_topo_map >= plateau_hgt);
    end
%     fprintf('*** Time cost by computing cond_land & cond_ocean: %.2f s ***\n', toc)    
    
    %/========== Output land/ocean/plateau masked data ===========%
    %/ Matrix (only logical)
    if ~isempty(matrix)
        cond_land_3D = repmat(cond_land, 1, 1, size(matrix,3));  %/ repmat to 3D for easy coding.
        output = matrix;
        if land_or_ocean == 1
            output(~cond_land_3D) = nan; %/ set oceanic vals to nan

        elseif land_or_ocean == 2
            output(~cond_ocean) = nan; %/ set land vals to nan
            
        elseif land_or_ocean == 3 
            output(~cond_plateau) = nan; %/ set non-plateau vals to nan
        else
            error('invalid input of ''land_or_ocean''!');
        end
    end
    
    %/ Points (indices or logical) (NOT for uneven lonlat grid cell!).
    if ~isempty(pos_lon_array) && flag_2D_lonlat == 0
        %/ Retrieve the grid resolution (it can be different for lon and lat)
        res_lon = unique(abs(diff(lon_grids)));
        res_lat = unique(abs(diff(lat_grids)));
    
        [lon_grids_2D, lat_grids_2D] = meshgrid(lon_grids, lat_grids);
        lon_grids_2D        = lon_grids_2D'; lat_grids_2D = lat_grids_2D';
        lon_grids_2Dto1D    = reshape(lon_grids_2D, [], 1);
        lat_grids_2Dto1D    = reshape(lat_grids_2D, [], 1);

        %/ Obtain all pairs of (lon,lat) over land
        if land_or_ocean == 1
            cond_land_2Dto1D    = reshape(cond_land,    [], 1);   %/ NoOfGrids x 1
            select_grids        = [lon_grids_2Dto1D(cond_land_2Dto1D), lat_grids_2Dto1D(cond_land_2Dto1D)];
            
        elseif land_or_ocean == 2
            cond_ocean_2Dto1D   = reshape(cond_ocean,   [], 1);   %/ NoOfGrids x 1
            select_grids        = [lon_grids_2Dto1D(cond_ocean_2Dto1D), lat_grids_2Dto1D(cond_ocean_2Dto1D)];
            
        elseif land_or_ocean == 3
            cond_plateau_2Dto1D = reshape(cond_plateau,   [], 1);   %/ NoOfGrids x 1
            select_grids        = [lon_grids_2Dto1D(cond_plateau_2Dto1D), lat_grids_2Dto1D(cond_plateau_2Dto1D)];
        else
            error('invalid input of ''land_or_ocean''!');
        end
        
        %/ Transpose if not a column vector
        if size(pos_lon_array, 1) == 1   
            pos_lon_array = pos_lon_array'; 
        end
        if size(pos_lat_array, 1) == 1   
            pos_lat_array = pos_lat_array'; 
        end
        pos_lon_array(pos_lon_array < 0) = pos_lon_array(pos_lon_array < 0) + 360; %/ Always make lon_grids in [0, 360).
        
        %/ Round up to the nearest grid (lon_grids, lat_grids)
        nearest_lon = round(pos_lon_array/res_lon, 0)*res_lon; %/ Feel free to validate this trick!
        nearest_lat = round(pos_lat_array/res_lat, 0)*res_lat;
        
%         a = [pos_lon_array, nearest_lon];  %/ checking
%         b = [pos_lat_array, nearest_lat];  %/ checking
        
        nearest_lon(nearest_lon == 360) = 0; %/ Restore 360E (if any) to 0E after round-up. 
        part_pos_nearest_grids = [nearest_lon, nearest_lat]; 
        
        if output_as_logical
            output = ismember(part_pos_nearest_grids, select_grids, 'rows');
        else
            output = find(ismember(part_pos_nearest_grids, select_grids, 'rows'));                    %/ find the matching rows (i.e. matching lon/lat of land grids!)
        end
    end
%     fprintf('*** Time cost by Output land/ocean masked data: %.2f s ***\n', toc)
%     fprintf('*** Time cost by ''show_land_or_ocean_hydrosheds'': %.5f s ***\n', toc)
end

%/ Alternatively, you can use landmask function from Matlab FileExchange. 

% lon_grids = [0:180];    res_lon = diff(lon_grids(1:2));
% lat_grids = [-45:45];   res_lat = diff(lat_grids(1:2));
% 
% [lon_2D,lat_2D] = meshgrid(lon_grids, lat_grids);
%                 
% z = zeros(size(lon_2D)) + 100;           
% 
% ocean = ~landmask(lat_2D,lon_2D); 
% land  = ~ocean;
% 
% z(ocean) = NaN; 
% 
% close all
% figure
% 
% map_proj = 'Miller Cylindrical'; 
% m_proj(map_proj,'longitudes',[min(lon_grids)  max(lon_grids)], 'latitudes', [min(lat_grids)  max(lat_grids)]);
% 
% m_pcolor(lon_2D-res_lon/2, lat_2D-res_lon/2, z)
% m_coast;
% m_grid;
    