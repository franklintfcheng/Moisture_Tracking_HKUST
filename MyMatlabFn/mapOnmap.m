%%
function freq_2D_nn = mapOnmap(varargin)
    
    pnames = {'data2D', 'lon', 'lat', 'dom_lon', 'dom_lat', 'unique_index'};  
    dflts  = cell(1, length(pnames));
    [          data2D,   lon,   lat,   dom_lon,   dom_lat,   unique_index] = ...
                        internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments

    %/=====================================================================
    %/ Author: Fandy Cheng
    %/ Date:   7 Mar 2022
    %/
    %/ Description:  'mapOnmap' is to map the original gridding (lon, lat) onto
    %/               the new gridding of the domain (dom_lon, dom_lat) without interpolation.
    %/
    %/               It is considered MORE accurate & intuitive THAN the Matlab
    %/               built-in function: interp('nearest'), despite it may be slower.
    %/
    %/ ------ input -------
    %/           data2D:     2D matrices with *original* (lon,lat) gridding
    %/                       Note that it ignores NaNs, but not zeros.
    %/
    %/ dom_lon, dom_lat:     *New* lon lat gridding of the domain (vectors)
    %/
    %/ ------ Output -------
    %/ freq_2D_nn:           A 2D frequency map on the new gridding.
    %/                       
    %/=====================================================================          
                    
    res_lon = unique(diff(dom_lon));
    res_lat = unique(diff(dom_lat));    %/ res_lat may not = res_lon.
    if length(res_lon) > 1     error('Non-uniform resolution of dom_lon detected!');   end
    if length(res_lat) > 1     error('Non-uniform resolution of dom_lat detected!');   end
    
    
    lon_range = [min(dom_lon) - res_lon/2, max(dom_lon) + res_lon/2];
    lat_range = [min(dom_lat) - res_lat/2; max(dom_lat) + res_lat/2];

    %/ [IMPORTANT] Get study domain boundaries for subsetting MCSs.
    study_dom_bndry = [lon_range(1), lat_range(1);
                       lon_range(1), lat_range(2);
                       lon_range(2), lat_range(2);
                       lon_range(2), lat_range(1);
                       lon_range(1), lat_range(1);];

    if isempty(find(isnan(data2D))) warning(' ''mapOnmap'': No NaN is found in data2D. The function handles all non-nans. See if this is what you want.');  end
    
    
    %/ [IMPORTANT] Find lon lat of *non-NaNs* of data2D
    [ind_lon_nonnan, ind_lat_nonnan] = find(~isnan(data2D));
    [in, ~] = inpolygon(lon(ind_lon_nonnan), lat(ind_lat_nonnan), study_dom_bndry(:,1), study_dom_bndry(:,2));   %/ 1D

    ind_lon_nonnan_in = ind_lon_nonnan(in);
    ind_lat_nonnan_in = ind_lat_nonnan(in);

    lon_in = nan(length(ind_lon_nonnan_in), 1);
    lat_in = nan(length(ind_lat_nonnan_in), 1);
    for i = 1:length(ind_lon_nonnan_in)
        lon_in(i)  = lon(ind_lon_nonnan_in(i));
        lat_in(i)  = lat(ind_lat_nonnan_in(i));
    end

    %/ Finally, map the subsetted lon, lat onto the new gridding of the domain
    [ind_lon, ind_lat] = point2map('pt_x', lon_in, 'pt_y', lat_in, 'lon', dom_lon, 'lat', dom_lat);
    
    if unique_index   %/ e.g., XX elements in a new grid cell will reduce to 1 elements. (save computation time!)
        [C, ~, ~] = unique([ind_lon, ind_lat], 'stable', 'rows');
        ind_lon = C(:,1);
        ind_lat = C(:,2);
    end
    
    %/ nn = nearest neighbor methods
    freq_2D_nn = zeros(length(dom_lon), length(dom_lat));
    
    for i = 1:length(ind_lon)
        freq_2D_nn(ind_lon(i), ind_lat(i)) = freq_2D_nn(ind_lon(i), ind_lat(i)) + 1;
    end
end