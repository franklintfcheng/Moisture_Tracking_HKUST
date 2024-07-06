%%

function [closest_grid, closest_ind] = find_nearest_val(grid_arr, pts, vartype)

    %/       a:   the reference grid array
    %/       b:   the point array to find closet values in a.
    %/ vartype:   if vartype == 'lon', we handle the lon data in [-179,180)

    %/ handle the lon data in [-179,180) --> more convenient for a range
    %/ across the prime meridian (0 deg)
    if nargin == 3 && isequal(vartype, 'lon')
        if ~isempty(grid_arr(grid_arr < 0))   %/ then assume in [-179, 180)
            pts(pts > 180) = pts(pts > 180) - 360;
            
        else                                  %/ otherwise assume in [0, 360)
            pts(pts < 0) = pts(pts < 0) + 360;
        end
%         grid_arr       = conv_to_lon_m179_180(grid_arr);
%         pts(pts > 180) = pts(pts > 180) - 360;
    end
    
    if any(diff(grid_arr) < 0)                              error('grid_arr after the conv_to_lon_m179_180 function should be strictly increasing!');  end 
    if size(grid_arr, 1) ~= 1 && size(grid_arr, 2) ~= 1     error('The input array is not a vector!');   end
    if size(pts, 1) ~= 1      && size(pts, 2) ~= 1          error('The input array is not a vector!');   end
    
    if size(grid_arr, 2) == 1     grid_arr = grid_arr';   end   %/ always a row vector 
    if size(pts, 1) == 1          pts      = pts';        end   %/ always a col vector 
    

    %/ Consider pts within the range.
    intvl  = abs(grid_arr(2)-grid_arr(1));  
    ind_in = find((pts < max(grid_arr) + intvl) & (pts > min(grid_arr) - intvl));
    
    %/ find the closest grids
    closest_grid = nan(size(pts));
    closest_ind  = nan(size(pts));
    [~, ind]        = min(abs(grid_arr - pts(ind_in)), [], 2); %/ (1 x N) - (M x 1) = (M x N) - Much faster than for-loop
    closest_grid(ind_in) = grid_arr(ind);
    closest_ind (ind_in) = ind;
end













