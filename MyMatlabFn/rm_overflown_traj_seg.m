function traj_data_new = rm_overflown_traj_seg(varargin)

    pnames = {'traj_data', 'dim_lon', 'dim_lat', 'map_lon_lower', 'map_lon_upper', 'map_lat_lower', 'map_lat_upper'};
    dflts  = cell(length(pnames), 1);
    [          traj_data,   dim_lon,   dim_lat,   map_lon_lower,   map_lon_upper,   map_lat_lower,   map_lat_upper ] ...
                 = internal.stats.parseArgs(pnames, dflts, varargin{:});

    %%
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: 2 Feb 2024
    %/
    %/ Description:
    %/   - The overall trick is to add [NaN  NaN] to break the traj while 
    %/     deleting the "overflown" segment of the traj outside the map.
    %/
    %/   - Assuming 'traj_data' a n x var matrix, with 'dim_lon', 'dim_lat' 
    %/     indicating the dimensions of lon and lat.
    %/=====================================================================
    
    %/ testing
    % map_lon_lower = 100; map_lon_upper = 125; map_lat_lower = 30; map_lat_upper = 40;
    % traj_data = [100, 30;
    %               70, 40;
    %               80, 50;
    %              120, 30;
    %              120, 25;
    %              134, 30;
    %              120, 35;];
    sz = size(traj_data);
    traj_data_new = traj_data;
    
    ind_overflown_seg = find(traj_data(:,dim_lon) < map_lon_lower | traj_data(:,dim_lon) > map_lon_upper | ...
                             traj_data(:,dim_lat) < map_lat_lower | traj_data(:,dim_lat) > map_lat_upper);
    
    traj_data_new(ind_overflown_seg, :) = nan(length(ind_overflown_seg), sz(2));

%     %/ 1. Check the points outside the *rightmost* edge.
%     ind = find(traj_data(:,1) > map_lon_upper);
%     if ~isempty(ind)
% %         warning('There are %d timesteps with lon greater than 180! Correct it to [-179, 180] and add NaNs to indicate the discontinuity', length(ind));
% 
%         traj_data(ind,1) = traj_data(ind,1) - lon_extent; %/ To avoid bug, do NOT write 360, since our plot is from -179 to 180 (due to a better viz of contf data).
%         nan_padding    = nan(length(ind), size(traj_data,2));
%         traj_data = insertrows(traj_data, nan_padding, ind); %/ [IMPORTANT] Insert a row of NaNs after each discontinuity point.
%     end
% 
%     %/ 2. Check the points outside the *leftmost* edge.
%     ind = find(traj_data(:,1) < map_lon_lower);
%     if ~isempty(ind)
% %         warning('There are %d particles whose lon is smaller than -179! Correct it to [-179, 180] and add NaNs to indicate the discontinuity', length(ind));
% 
%         traj_data(ind,1) = traj_data(ind,1) + lon_extent; %/ To avoid bug, do NOT write 360, since our plot is from -179 to 180 (due to a better viz of contf data).
%         nan_padding    = nan(length(ind), size(traj_data,2));
%         traj_data = insertrows(traj_data, nan_padding, ind); %/ [IMPORTANT] Insert a row of NaNs after each discontinuity point.
%     end
% 
%     %/ 3. At last, check if there are pairs of points that travel across the dateline, we have to add NaNs to separate them. 
%     %/    CAREFUL: the index of insertion is different between diff > 100 and diff < -100.
%     ind_across_R2L = find(diff(traj_data(:,1)) < -200);
%     nan_padding    = nan(length(ind_across_R2L), size(traj_data, 2));
%     traj_data = insertrows(traj_data, nan_padding, ind_across_R2L); %/ [IMPORTANT] Insert a row of NaNs between these pairs.
% 
%     ind_across_L2R = find(diff(traj_data(:,1)) > 200);
%     nan_padding    = nan(length(ind_across_L2R), size(traj_data, 2));
%     traj_data = insertrows(traj_data, nan_padding, ind_across_L2R); %/ [IMPORTANT] Insert a row of NaNs between these pairs.

end