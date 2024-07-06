%%

function date_LST_2D = create_2D_LST(date_UTC, data_2D, lon, lat)

    %/ NOTE: This function is to create a 2D map of local solar time (LST)
    %/       based on the longitudes of non-NaN grids in data_2D.

    [lon_2D, ~]         = meshgrid(lon, lat);
    lon_2D              = lon_2D';
    
    if size(data_2D) ~= size(lon_2D)   error('Inconsistent Dimensions!');    end
    
    lon_array           = nan(size(data_2D));
    cond                = ~isnan(data_2D);                  %/ 2D logical
    lon_array(cond)     = lon_2D(cond);                     %/ subset only the lon (on 2D map) related to MCS as it's triggering max prcp.
    date_LST_2D         = UTCtoLST(date_UTC, lon_array);    %/ estimate the local solar time (LST)
    date_LST_2D         = datetime2int(date_LST_2D);

end