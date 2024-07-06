%%
function [output, cond_land, cond_ocean] = show_land_or_ocean(varargin)
    
    pnames = {'matrix', 'pos_lon_array', 'pos_lat_array', 'lon_grids', 'lat_grids',...
              'land_or_ocean', 'plateau_hgt', 'output_as_logical', 'cond_land', 'cond_ocean'};
    dflts  = cell(1, length(pnames));
    
    [          matrix,    pos_lon_array,   pos_lat_array,   lon_grids,   lat_grids,...
               land_or_ocean,   plateau_hgt,   output_as_logical,   cond_land,  cond_ocean] = internal.stats.parseArgs(pnames, dflts, varargin{:});
    
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
    error('For consistency and accuracy, use ''show_land_or_ocean_hydrosheds()'' instead!')


    if ~isempty(matrix) && (~isempty(pos_lon_array) || ~isempty(pos_lat_array))
        error('Either input a ''matrix'' or a pair of ''pos_lon_array'' and ''pos_lat_array''!');
    end
    
    if isempty(matrix) && isempty(pos_lon_array) && isempty(pos_lat_array)
        %/ output will be empty.
        %/ But will output cond_land & cond_ocean based on lon_grids & lat_grids.
        output = [];
    end
    
    if isempty(plateau_hgt)  plateau_hgt = 3000;  end

    %/========== Obtain Topo Map (now only for cond_plateau) ===========%
    %/ Based on the input lon/lat grids averaged from 5-min topography.
    lon_grids(lon_grids < 0) = lon_grids(lon_grids < 0) + 360;  %/ always make lon_grids in [0, 360).

    %/ Retrieve the grid resolution (it can be different for lon and lat)
    res_lon = unique(abs(diff(lon_grids)));
    res_lat = unique(abs(diff(lat_grids)));
    
    if isempty(cond_land) || isempty(cond_ocean)
        avg_topo_map = interp_topo('lon_grids', lon_grids, 'lat_grids', lat_grids);
        
        if ~isempty(find(isnan(avg_topo_map)))  error('avg_topo contains NaN!');   end
        cond_land    = (avg_topo_map > 0);  %/ set >0 instead of >= 0 can avoid mispresenting sea-ice as land.
        cond_ocean   = (avg_topo_map <= 0);
    end
    
    if land_or_ocean == 3
        avg_topo_map = interp_topo('lon_grids', lon_grids, 'lat_grids', lat_grids);
        
        if ~isempty(find(isnan(avg_topo_map)))  error('avg_topo contains NaN!');   end
        cond_plateau   = (avg_topo_map >= plateau_hgt);
    end

%     fprintf('*** Time cost by computing cond_land & cond_ocean: %.2f s ***\n', toc)
    
    
    %/========== Output land/ocean masked data ===========%
%     tic;
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
    
    %/ Points (indices or logical) 
    if ~isempty(pos_lon_array)
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
        if size(pos_lon_array, 1) == 1   pos_lon_array = pos_lon_array'; end
        if size(pos_lat_array, 1) == 1   pos_lat_array = pos_lat_array'; end
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
    
%     fprintf('*** Time cost by ''show_land_or_ocean'': %.5f s ***\n', toc)
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
    