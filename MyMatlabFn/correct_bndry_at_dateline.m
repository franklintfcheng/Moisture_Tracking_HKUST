%%
function bndry_data = correct_bndry_at_dateline(bndry_data, map_lon_upper, map_lon_lower)

    %/ The overall trick is to add [NaN  NaN] to break the verge without
    %/ changing the coordinates.
    
    %/ CAVEAT:
    %/ This function is not yet able to correct the boundary data that is already set
    %/ by [-179, 180] in a different global lon range (say, [-169 190]) 

    lon_extent = map_lon_upper - map_lon_lower; %/ ~= 359 if global
    
    %/ 1. Check the points outside the *rightmost* edge.
    ind = find(bndry_data(:,1) > map_lon_upper);
    if ~isempty(ind)
%         warning('There are %d timesteps with lon greater than 180! Correct it to [-179, 180] and add NaNs to indicate the discontinuity', length(ind));

        bndry_data(ind,1) = bndry_data(ind,1) - lon_extent; %/ To avoid bug, do NOT write 360, since our plot is from -179 to 180 (due to a better viz of contf data).
        nan_padding    = nan(length(ind), size(bndry_data,2));
        bndry_data = insertrows(bndry_data, nan_padding, ind); %/ [IMPORTANT] Insert a row of NaNs after each discontinuity point.
    end

    %/ 2. Check the points outside the *leftmost* edge.
    ind = find(bndry_data(:,1) < map_lon_lower);
    if ~isempty(ind)
%         warning('There are %d particles whose lon is smaller than -179! Correct it to [-179, 180] and add NaNs to indicate the discontinuity', length(ind));

        bndry_data(ind,1) = bndry_data(ind,1) + lon_extent; %/ To avoid bug, do NOT write 360, since our plot is from -179 to 180 (due to a better viz of contf data).
        nan_padding    = nan(length(ind), size(bndry_data,2));
        bndry_data = insertrows(bndry_data, nan_padding, ind); %/ [IMPORTANT] Insert a row of NaNs after each discontinuity point.
    end

    %/ 3. At last, check if there are pairs of points that travel across the dateline, we have to add NaNs to separate them. 
    %/    CAREFUL: the index of insertion is different between diff > 100 and diff < -100.
    ind_across_R2L = find(diff(bndry_data(:,1)) < -200);
    nan_padding    = nan(length(ind_across_R2L), size(bndry_data, 2));
    bndry_data = insertrows(bndry_data, nan_padding, ind_across_R2L); %/ [IMPORTANT] Insert a row of NaNs between these pairs.

    ind_across_L2R = find(diff(bndry_data(:,1)) > 200);
    nan_padding    = nan(length(ind_across_L2R), size(bndry_data, 2));
    bndry_data = insertrows(bndry_data, nan_padding, ind_across_L2R); %/ [IMPORTANT] Insert a row of NaNs between these pairs.

end