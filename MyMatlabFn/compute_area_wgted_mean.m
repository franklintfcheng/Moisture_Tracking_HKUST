function [AWM, M, area] = compute_area_wgted_mean(varargin)

    pnames = {'data', 'lon', 'lat', 'lon_extent', 'lat_extent'};

    dflts = repmat([], 1, length(pnames));

    [          data,   lon,   lat,   lon_extent,   lat_extent] = ...
                        internal.stats.parseArgs(pnames, dflts, varargin{:});
    %%
    %==================================================================================================
    %/       Author: Franklin Cheng (fandycheng@ust.hk)
    %/  Last Update: Mar 21, 2024
    %/
    %/       Output: area-weighted mean
    %==================================================================================================
    sz = size(data);
    if length(sz) > 2
        error('data can only be 2D!');
    end

    %/ If lon_extent starts/ends with -ve values, convert lon to [-179, 180]
    if any(lon_extent < 0)
        lon(lon > 180) = lon(lon > 180) - 360;
    end
    % ind = length(find(lon >= min(lon_extent) & (lo(1n <= max(lon_extent))))

    if isvector(lon) && isvector(lat)
        [lon_2D, lat_2D] = meshgrid(lon, lat);
        lon_2D = lon_2D';
        lat_2D = lat_2D';
    else
        lon_2D = lon;
        lat_2D = lat;
    end

    if ~isempty(lon_extent) && ~isempty(lat_extent)
        logi_lon_2D  = (lon_2D >= min(lon_extent)) & (lon_2D <= max(lon_extent));
        logi_lat_2D  = (lat_2D >= min(lat_extent)) & (lat_2D <= max(lat_extent));
    else
        logi_lon_2D = true(size(lon_2D));  %/ Otherwise, compute the weighted average for the entire matrix
        logi_lat_2D = true(size(lat_2D));  
    end
    
    area = calc_grid_area('lon', lon, 'lat', lat);
    area(~(logi_lon_2D & logi_lat_2D)) = nan;

    data_bc = data;
    data_bc(~(logi_lon_2D & logi_lat_2D)) = nan;

    wgt = area./sum(area, 'all', 'omitnan');
    AWM = sum(data_bc.*wgt, 'all', 'omitnan'); %/ Area-weighted arithmetric mean

    M = mean(data_bc, 'all', 'omitnan'); %/ Arithmetric mean (just for checking how different it is from AWM)

    %/ Double checking
    % figure
    % imagesc(area)
    % min(area, [],'all','omitnan')
    % max(area, [],'all','omitnan')

    % 
    % figure
    % imagesc(data_bc)
    % 
    % min(wgt,[],'all','omitnan')
    % max(wgt,[],'all','omitnan')
    % sum(wgt,'all','omitnan')
    % size(wgt)

end
