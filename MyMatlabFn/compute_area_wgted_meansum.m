function [AWM, M, area, area_unit, str_extent] = compute_area_wgted_meansum(varargin)

    pnames = {'data', 'lon', 'lat', 'lon_extent', 'lat_extent', 'mean_or_sum'};

    dflts  = cell(length(pnames), 1);

    [          data,   lon,   lat,   lon_extent,   lat_extent,   mean_or_sum] ...
                      = internal.stats.parseArgs(pnames, dflts, varargin{:});
    %%
    %==================================================================================================
    %/       Author: Franklin Cheng (fandycheng@ust.hk)
    %/  Last Update: Jul 5, 2024
    %/
    %/       Output: area-weighted mean or sum
    %==================================================================================================
    sz = size(data);
    if length(sz) > 2
        error('data can only be 2D!');
    end

    if isempty(lon_extent) && isempty(lat_extent)
        lon_extent = [min(lon), max(lon)];
        lat_extent = [min(lat), max(lat)];
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
    
    % if ~isempty(lon_extent) && ~isempty(lat_extent)
    logi_lon_2D  = (lon_2D >= min(lon_extent)) & (lon_2D <= max(lon_extent));
    logi_lat_2D  = (lat_2D >= min(lat_extent)) & (lat_2D <= max(lat_extent));
    % else
    %     logi_lon_2D = true(size(lon_2D));  %/ Otherwise, compute the weighted average for the entire matrix
    %     logi_lat_2D = true(size(lat_2D));  
    % end

    area = calc_grid_area('lon', lon, 'lat', lat);
    area(~(logi_lon_2D & logi_lat_2D)) = nan;

    data_bc = data;
    data_bc(~(logi_lon_2D & logi_lat_2D)) = nan;

    if isempty(mean_or_sum)
        mean_or_sum = 'mean';
    elseif ~isequal(mean_or_sum, 'sum') && ~isequal(mean_or_sum, 'mean')
        error('invalid input of ''mean_or_sum''. Set it to ''sum'' or ''mean''.');
    end

    if isequal(mean_or_sum, 'sum')
        AWM = sum(data_bc.*area, 'all', 'omitnan')*1e-6; %/ Multiply by area (km^2)
        M   = [];
        area_unit = 'km^{2}';

    elseif isequal(mean_or_sum, 'mean')
        wgt = area./sum(area, 'all', 'omitnan');
        AWM = sum(data_bc.*wgt, 'all', 'omitnan'); %/ Area-weighted arithmetric mean
        M   = mean(data_bc, 'all', 'omitnan'); %/ Arithmetric mean (just for checking how different it is from AWM)
        area_unit = '';
    end

    if lon_extent(1) > 0 
        lon1_unit = 'E'; 
    elseif lon_extent(1) < 0 
        lon1_unit = 'W';
    else
        lon1_unit = '';
    end
    if lon_extent(2) > 0 
        lon2_unit = 'E'; 
    elseif lon_extent(2) < 0 
        lon2_unit = 'W';
    else
        lon2_unit = '';
    end
    if lat_extent(1) > 0 
        lat1_unit = 'N'; 
    elseif lat_extent(1) < 0 
        lat1_unit = 'S';
    else
        lat1_unit = '';
    end
    if lat_extent(2) > 0 
        lat2_unit = 'N'; 
    elseif lat_extent(2) < 0 
        lat2_unit = 'S';
    else
        lat2_unit = '';
    end
    str_extent = sprintf(' (%s%s-%s%s, %s%s-%s%s)', num2str(abs(lon_extent(1))), lon1_unit, num2str(abs(lon_extent(2))), lon2_unit,...
                                                    num2str(abs(lat_extent(1))), lat1_unit, num2str(abs(lat_extent(2))), lat2_unit);

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
