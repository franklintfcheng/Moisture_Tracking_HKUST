%%
function [avg_topo_map] = interp_topo(varargin)
    
    pnames = {'lon_grids', 'lat_grids'};
    dflts  = cell(1, length(pnames));
    
    [         lon_grids,   lat_grids] = internal.stats.parseArgs(pnames, dflts, varargin{:});
%%
    %/========================================================================================
    %/ Author: Fandy CHENG
    %/ Date of last update: 26 Jan 2023
    %/
    %/ This function is considered *more robust and intuitive* than 'interp('nearest')' function.
    %/========================================================================================
    
    %/ Based on the input lon/lat grids averaged from 5-min topography.
    lon_grids_bc = lon_grids;
    lat_grids_bc = lat_grids;
    lon_grids_bc(lon_grids_bc < 0) = lon_grids_bc(lon_grids_bc < 0) + 360;  %/ always make lon_grids in [0, 360).
    
    res_lon = 'no_res';
    res_lat = 'no_res';
    flag_2D_lonlat = 0;
    if isvector(lon_grids_bc) && isvector(lat_grids_bc)
        %/ If lon_grids, lat_grids are 1D vectors
        nlon = length(lon_grids_bc);
        nlat = length(lat_grids_bc);
        
        %/ Turn into a column vector
        if size(lon_grids_bc, 1) == 1  lon_grids_bc = lon_grids_bc';  end
        if size(lat_grids_bc, 1) == 1  lat_grids_bc = lat_grids_bc';  end   

        %/ Retrieve the grid resolution (it can be different for lon and lat)
        if numel(lon_grids_bc) ~= 1  res_lon = unique(abs(diff(lon_grids_bc(1:2))));  end 
        if numel(lat_grids_bc) ~= 1  res_lat = unique(abs(diff(lat_grids_bc(1:2))));  end
    else
        %/ If lon_grids, lat_grids are 2D matrices,
        flag_2D_lonlat = 1;
        [nlon, nlat]   = size(lon_grids_bc); %/ Assume in (lon, lat)
        
        res_lon = min(abs(diff(lon_grids_bc,[],1)),[],[1,2]); %/ get the highest resolution from the uneven lon grid
        res_lat = min(abs(diff(lat_grids_bc,[],2)),[],[1,2]); %/ get the highest resolution from the uneven lat grid
        
    end
    
    if isequal(res_lon, 'no_res') && isequal(res_lat, 'no_res')
        error('No resolution can be determined from the input ''lon_grids_bc'' and ''lat_grids_bc''!');
    else
        if isequal(res_lon, 'no_res')   res_lon = res_lat;   end %/ *borrow* the res from another dim.
        if isequal(res_lat, 'no_res')   res_lat = res_lon;   end %/ *borrow* the res from another dim.
    end
    
%     %/ ==== Special treatment of lon_grids ====
%     %/ Sometimes we may input lon_grids in [-169:190],
%     %/ res_lon returns [1, 359], but since the array is not strictly increasing --> this is fine.
%     if length(res_lon) ~= 1     
%         if issorted(lon_grids_bc) == 0
%             res_lon = min(res_lon); %/ take the smallest one.
%         else
%             disp(res_lon)
%             error('res_lon is not constant! Check your input lon_grids!'); 
%         end
%     end
%     if length(res_lat) ~= 1     error('res_lat is not constant! Check your input lat_grids!'); end
    
    m_coord('geographic');                                      %/ required before calling m_tbase.
    [topo, lon_topo_2D, lat_topo_2D] = m_tbase([0 360 -90 90]); %/ whatever input lon_grids_bc lat_grids_bc are, get the global topo (5-min) first.
    topo     = topo'; 
    lon_topo = lon_topo_2D(1,:); 
    lat_topo = lat_topo_2D(:,1);
    ind_360           = find(lon_topo == 360);  %/ remove lon_grids_bc = 360. (since 360E == 0E).
    topo(ind_360, :)  = [];
    lon_topo(ind_360) = [];
    lon_topo_2D(:,ind_360) = [];
    lat_topo_2D(:,ind_360) = [];
%     unique(diff(lon_topo))
%     unique(diff(lat_topo))
    
    if res_lon < abs(diff(lon_topo(1:2))) || res_lat < abs(diff(lat_topo(1:2))) || flag_2D_lonlat 
        if flag_2D_lonlat
            fprintf('Detected lon lat are 2D matrices; automatically perform interpolation of topo to this gridding\n')
        else
            warning('The queried lon_grids_bc res is even higher than the topo lon_grids_bc res (%.2f)!', abs(diff(lon_topo(1:2))));
            warning('This leaves us no choices but simply interpolating topo onto the target gridding..')
        end
        
        if flag_2D_lonlat
            target_lon_2D = lon_grids_bc; %/(lon, lat)
            target_lat_2D = lat_grids_bc;
            
            %/ Use 'griddata' instead of 'interp2' for uneven gridding!!
            avg_topo_map = griddata(lon_topo_2D, lat_topo_2D, topo', target_lon_2D', target_lat_2D', 'linear')';
        else
            [target_lon_2D, target_lat_2D] = meshgrid(lon_grids_bc, lat_grids_bc);
            %/ Use 'griddata' instead of 'interp2' for uneven gridding!!
            avg_topo_map = griddata(lon_topo_2D, lat_topo_2D, topo', target_lon_2D, target_lat_2D, 'linear')';
        end
        
         
    else
        %/ *Average* all topo grids within the queried grid cel (Dumpest way)
        avg_topo_map = nan(nlon, nlat);
        for i = 1:nlon
            %/ [ lon1 )[ lon2 )[ lon3 )[ lon4 )...
            ind_lon = find(lon_grids_bc(i)-res_lon/2 <= lon_topo & lon_topo < lon_grids_bc(i)+res_lon/2);
            if isempty(ind_lon)
                error('empty ind_lon!');
            end

            if lon_grids_bc(i)-res_lon/2 < 0
    %             fprintf('< 0: i = %d \n', i)
                ind_lon = [ind_lon, find(lon_grids_bc(i)-res_lon/2+360 <= lon_topo)];  %/ append the missing ones at lon_grids_bc lower boundary.

            elseif lon_grids_bc(i)+res_lon/2 > 360  
    %             fprintf('> 360: i = %d \n', i)
                ind_lon = [ind_lon, find(lon_topo < lon_grids_bc(i)+res_lon/2-360)];   %/ append the missing ones at lon_grids_bc upper boundary.
            end

            for j = 1:nlat
                %/ [ lat1 )[ lat2 )[ lat3 )[ lat4 )...
                ind_lat = find(lat_grids_bc(j)-res_lat/2 <= lat_topo & lat_topo < lat_grids_bc(j)+res_lat/2); %/ NOTE: at poles, the # of topo grids will halve.
                avg_topo_map(i,j) = mean(topo(ind_lon, ind_lat), 'all');
            end
        end
    end
    %/ Convert the lon/lat ordering format of avg_topo_map into that of the
    %/ enquired lon_grids, lat_grids.
%     avg_topo_map = avg_topo_map(I_sorted_lon,:);
%     avg_topo_map = avg_topo_map(:,I_sorted_lat);

end