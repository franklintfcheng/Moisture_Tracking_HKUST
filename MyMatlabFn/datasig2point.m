
function point_data = datasig2point(data_sig, lon, lat)

    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: Feb 28, 2024
    %/
    %/ Description: This function is designed to convert a 2D data map with
    %/              significant values (i.e. non-NaN values) into a
    %/              point_data with lon and lat positions of the corresponding
    %/              non-NaN values.
    %/        
    %/              The output 'point_data' fits the input of the 'plot_contfmap' function.
    %/=====================================================================       

    [lon_2D, lat_2D] = meshgrid(lon, lat);
    lon_2D = lon_2D'; 
    lat_2D = lat_2D';

    data_sig_2Dto1D = reshape(data_sig, [], 1);
    lon_2Dto1D      = reshape(lon_2D, [], 1);
    lat_2Dto1D      = reshape(lat_2D, [], 1);

    cond_sig = ~isnan(data_sig_2Dto1D);
    lon_2Dto1D = lon_2Dto1D(cond_sig);
    lat_2Dto1D = lat_2Dto1D(cond_sig);

    point_data = [lon_2Dto1D, lat_2Dto1D];
end