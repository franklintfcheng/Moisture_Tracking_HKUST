%%
function [output, hori_array] = compute_zm_mm(varargin)
    
    pnames       = {   'geo_data',       'lon',  'lat',    'lon_range',    'lat_range',      'zm_or_mm',  'mean_or_sum'};
    dflts        = cell(length(pnames), 1);
    [                   geo_data,         lon,    lat,      lon_range,      lat_range,        zm_or_mm,    mean_or_sum] = ...
            internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments
    
    %%
    %/ 'geo_data': can be 3D or 4D.
    fprintf('*** Running compute_zm_mm... ***\n')
    Ndim = length(size(geo_data));
    
    if isempty(mean_or_sum)     mean_or_sum = 1;    end %/ take zonal / merid. mean by default.
    if size(lat, 2) ~= 1        lat = lat';         end
    
    if size(geo_data, 1) ~= length(lon)   error('lon dimension inconsitent!');  end
    if size(geo_data, 2) ~= length(lat)   error('lat dimension inconsitent!');  end
    
    ind_lon = find(lon >= min(lon_range) & lon <= max(lon_range));
    ind_lat = find(lat >= min(lat_range) & lat <= max(lat_range));
    
    lon_new      = lon(ind_lon);
    lat_new      = lat(ind_lat);
    
    if Ndim == 2
        geo_data_new = geo_data(ind_lon, ind_lat);
    elseif Ndim == 3
        geo_data_new = geo_data(ind_lon, ind_lat, :);
    elseif Ndim == 4
        geo_data_new = geo_data(ind_lon, ind_lat, :, :);
    end
    
    %/ Only for meridional average should we do latitude-weighted average!
    if mean_or_sum == 1
        if zm_or_mm == 2
            wgt = cosd(lat_new);                                           %/ Create area weights first because the grid area decreases with latitude
            wgt = repmat(wgt', length(lon_new), 1);                        %/ Reshape to perform 3D element-wise mult.
            wgt = wgt./sum(wgt, 2, 'omitnan');                             %/ Normalize wgt so that its sum along each longitude is 1 

            output = squeeze(sum(geo_data_new.*wgt, zm_or_mm, 'omitnan')); %/ since wgt has been normalized, here we use sum (or nansum) to get the weighted mean.

        elseif zm_or_mm == 1
            %/ grids along the same latitude have an equal area -> arithmetic average.
            output = squeeze(mean(geo_data_new, zm_or_mm, 'omitnan'));
        end
        
    elseif mean_or_sum == 2
        %/ Do NOT apply area weighting to zonal / meridional sum.
        output = squeeze(sum(geo_data_new, zm_or_mm, 'omitnan'));
    end
    
    if     zm_or_mm == 1      hori_array = lat_new;
    elseif zm_or_mm == 2      hori_array = lon_new;           end
    
    % fprintf('*** Finished compute_zm_mm. ***\n')
%     size(output)
end