%%
function ind = which_on_land(varargin)

    pnames = {'pos_lon_array', 'pos_lat_array', 'pos_hgt_array', 'hgt_range', 'lon_topo', 'lat_topo', 'full_logic_array', 'land_grids'};
    dflts  = cell(1, length(pnames));
    [pos_lon_array, pos_lat_array, pos_hgt_array, hgt_range, lon_topo, lat_topo, full_logic_array, land_grids] = internal.stats.parseArgs(pnames, dflts, varargin{:});
    
    %/ WARNING: This function is designed for 1x1 degree mapping, despite that
    %/          here it used very high res. (5-min) topo data to verify.
    %/          Thus, the default lon_topo and lat_topo would be 0:359 and -89:89
    
    %/          Tips: if you use this function in a big for-loop, suggest feeding the function with land_grids preloaded. 
    %/                It avoids reading data again and again. 
    
    error('The function ''which_on_land'' was deprecated. Use ''show_land_or_ocean'' function instead! ')
    
    if (isempty(lon_topo) && isempty(lat_topo)) || (isequal(lon_topo, 0:359) && isequal(lat_topo, -89:89))
        lon_topo = 0:359;    
        lat_topo = -89:89;
        
        %/ NOTE:
        %/ load the 5-minutes topo and subset it to 1-degree!
        %/ load the processed data --> MUCH faster than calling m_tbase !!!
        %/ input directly land_grids from outside should be even faster !!!
        if isempty(land_grids)
            filename = '/disk/r059/tfchengac/FLEXPART/domfill_CERA_MPI/prcssd_data_4plotting/land_grids.mat';
            if isfile(filename)
                fprintf('*** Loading %s ***\n', filename)
                load(filename, 'land_grids');
            else
                m_coord('geographic');   %/ call this to avoid "No coordinate system initialized" when calling m_tbase for the first time.
                [topo, lon_topo_2D, lat_topo_2D] = m_tbase([lon_topo(1) lon_topo(end) lat_topo(1) lat_topo(end)]); %/ Its lat dim is from 89 to -89 !!
                ind_lon = findismember_loop(lon_topo_2D(1,:),  lon_topo);                                          %/ Use findismember_loop to avoid auto sorting !!
                ind_lat = findismember_loop(lat_topo_2D(:,1)', lat_topo);
                topo = topo(ind_lat, ind_lon)';
                [ind_lon_land, ind_lat_land] = find(topo >= 0);
                land_grids = [lon_topo(ind_lon_land); lat_topo(ind_lat_land)]';
                save(filename, 'land_grids');
            end
        end
    else
        error('This function is not set for the input lon_topo and lat_topo!!')
    end

    
    %/ Transpose it if it is not a column vector
    if size(pos_lon_array, 1) == 1   pos_lon_array = pos_lon_array'; end
    if size(pos_lat_array, 1) == 1   pos_lat_array = pos_lat_array'; end
    if size(pos_hgt_array, 1) == 1   pos_hgt_array = pos_hgt_array'; end
    
    %/ IMPORTANT: Convert to [0 359] (otherwise I have to consider an opposite way to round the lon location in the western hemisphere!)
    if ~isempty(find(pos_lon_array < 0))
        pos_lon_array(pos_lon_array < 0) = pos_lon_array(pos_lon_array < 0) + 360;
    end

    %/ round up to the nearest grid (assume 1x1 deg!!)
    nearest_lon = round(pos_lon_array, 0);
    nearest_lat = round(pos_lat_array, 0);
    nearest_lon(nearest_lon == 360) = 0;                 %/ set it back to 0E if being rounded up to 360E !! 
    part_pos_nearest_grids = [nearest_lon, nearest_lat]; 

%     tic;
    if ~isempty(pos_hgt_array) && ~isempty(hgt_range)    %/ then we also select those within the height range.
        if full_logic_array
            ind = ismember(part_pos_nearest_grids, land_grids, 'rows') & ...                         %/ find the matching rows (i.e. matching lon/lat of land grids!)
                                    pos_hgt_array >= hgt_range(1) & pos_hgt_array <= hgt_range(2);   %/ and below ?m      
        else
            
            ind = find(ismember(part_pos_nearest_grids, land_grids, 'rows') & ...                    %/ find the matching rows (i.e. matching lon/lat of land grids!)
                                    pos_hgt_array >= hgt_range(1) & pos_hgt_array <= hgt_range(2));  %/ and below ?m         
        end
    else
        if full_logic_array
            ind = ismember(part_pos_nearest_grids, land_grids, 'rows');
        else
            ind = find(ismember(part_pos_nearest_grids, land_grids, 'rows'));                    %/ find the matching rows (i.e. matching lon/lat of land grids!)
        end
    end
%     toc;

end

%     lon_topo = 0.5:359.5;
%     lat_topo = -89.5:89.5;
%     [topo, ~, ~] = m_elev([lon_topo(1) lon_topo(end) lat_topo(1) lat_topo(end)]);
%     topo = topo';
%     [ind_lon_land, ind_lat_land] = find(topo >= 0);
%     land_grids = [lon_topo(ind_lon_land); lat_topo(ind_lat_land)]';

    
    %/ As topo_lon & topo_lat start from 0.5 and -89.5, we first offset
    %/ particle's position by 0.5 and offset it back after rounding.
%     offset = 0.5;   
%     nearest_lon = round((pos_lon_array + offset), 0) - offset; %/ round up to the nearest int
%     nearest_lat = round((pos_lat_array + offset), 0) - offset; %/ round up to the nearest int
%     
%     nearest_lon(nearest_lon == 360.5) = 0.5;         %/ set it back to 0.5 deg if being rounded up to 180.5 deg!! 
