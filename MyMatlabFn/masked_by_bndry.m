%%
function X_masked = masked_by_bndry(varargin)

    pnames = {'X', 'lon', 'lat', 'bndry_data', 'output_logical'};
    dflts =  {[],   [],     [],            [],               0 };
    [X, lon, lat, bndry_data, output_logical] = internal.stats.parseArgs(pnames, dflts, varargin{:});

    %==================================================================================================
    %/       Author: Franklin Cheng (fandycheng@ust.hk)
    %/  Last Update: October 10, 2023
    %/
    %/  Description: This function returns the matrix with values enclosed
    %/                  by the boundary data (bndry_data) only. 
    %/                  Values outside the boundary are set NaN.
    %==================================================================================================
    
    [lon_2D, lat_2D] = meshgrid(lon, lat);
    lon_2Dto1D = reshape(lon_2D', [], 1);
    lat_2Dto1D = reshape(lat_2D', [], 1);
    
    [in, ~] = inpoly2([lon_2Dto1D, lat_2Dto1D], bndry_data);  %/ 1D
    ind_out = find(in == 0);

    X_1D = reshape(X, [], 1);                             %/ reshape to 1D
    X_1D(ind_out) = nan;                   
    X_masked = reshape(X_1D, length(lon), length(lat));   %/ reshape back to 2D

    if output_logical
        X_masked(isnan(X_masked)) = 0;
        X_masked = logical(X_masked);
    end
end