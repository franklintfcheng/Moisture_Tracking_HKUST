%%
function lon_m179_180 = conv_to_lon_m179_180(lon)

    %/ Convert lon into the range of [-179, 180)
    flag = 0;
    if size(lon, 1) == 1    lon = lon';  flag = 1;  end   %/ convert into a column vector

    lon_m179_180      = lon;
    ind               = find(lon_m179_180 > 180);
    
    if ~isempty(ind)
        lon_m179_180(ind) = lon_m179_180(ind) - 360;                 %/ change to [0:180,-179:0] for a correct contf plot.
        ind_translate     = [ind', 1:ind(1)-1];       %/ put the left lon dim to the right. --> [-179:0, 0:180]
        lon_m179_180      = lon_m179_180(ind_translate);
    end

    if flag
        lon_m179_180 = lon_m179_180';  %/ transpose to the original dimensions.
    end
end