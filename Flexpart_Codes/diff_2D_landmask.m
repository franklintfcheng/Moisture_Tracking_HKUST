function [diff_data, SSIM, r, pval, RMSE, diff_data_lon, diff_data_lat] = diff_2D_landmask(varargin)

    pnames = {'data1', 'data1_lon', 'data1_lat', 'data2', 'data2_lon', 'data2_lat',...
              'show_prct_diff', 'domain_lon_range', 'domain_lat_range', 'interp_method', 'land_or_ocean',...
              'show_TP_only', 'data_folder'};
    dflts  = cell(1, length(pnames));
    [         data1, data1_lon, data1_lat, data2, data2_lon, data2_lat,...
              show_prct_diff, domain_lon_range, domain_lat_range, interp_method, land_or_ocean,...
              show_TP_only,    data_folder] = internal.stats.parseArgs(pnames, dflts, varargin{:});

    %%
    %/ Description
    %/ 'show_prct_diff':   1]: (data1 - data2)/data2*100%;  0]: data1 - data2
    data1_bc     = data1;
    data1_lon_bc = data1_lon;
    data1_lat_bc = data1_lat;
    data2_bc     = data2;
    data2_lon_bc = data2_lon;
    data2_lat_bc = data2_lat;
    
    %/ Subset the data based on domain_lon_range, domain_lat_range
    if isempty(domain_lon_range)   
        warning('''domain_lon_range'' not specified! Assume it to be the range of data1_lon.');   
        domain_lon_range = [min(data1_lon), max(data1_lon)];
    end
    if isempty(domain_lat_range)   
        warning('''domain_lat_range'' not specified! Assume it to be the range of data1_lat.');   
        domain_lat_range = [min(data1_lat), max(data1_lat)];
    end

    if isvector(data1_lon_bc) && isvector(data1_lat_bc)
        %/ transpose to column vetor
        if size(data1_lon_bc, 1) == 1      data1_lon_bc = data1_lon_bc';   end
        if size(data1_lat_bc, 1) == 1      data1_lat_bc = data1_lat_bc';   end
        if any(diff(data1_lat_bc) < 0)
            warning('Correct the lat of data 1 to be from -89 to +89.')
            data1_bc     = flip(data1_bc, 2);
            data1_lat_bc = flip(data1_lat_bc, 1);
        end
        [data1_lon_bc_2D, data1_lat_bc_2D] = meshgrid(data1_lon_bc, data1_lat_bc);
    else
        data1_lon_bc_2D = data1_lon_bc;
        data1_lat_bc_2D = data1_lat_bc;
    end
    
    if isvector(data2_lon_bc) && isvector(data2_lat_bc)
        if size(data2_lon_bc, 1) == 1      data2_lon_bc = data2_lon_bc';   end
        if size(data2_lat_bc, 1) == 1      data2_lat_bc = data2_lat_bc';   end
        if any(diff(data2_lat_bc) < 0)
            warning('Correct the lat of data 2 to be from -89 to +89.')
            data2_bc     = flip(data2_bc, 2);
            data2_lat_bc = flip(data2_lat_bc, 1);
        end
        [data2_lon_bc_2D, data2_lat_bc_2D] = meshgrid(data2_lon_bc, data2_lat_bc); %/ lat x lon
    else
        data2_lon_bc_2D = data2_lon_bc';  %/ lat x lon
        data2_lat_bc_2D = data2_lat_bc';  %/ lat x Lon 
    end
    lon_dim = 2;
    
%     ind_lon      = find(data1_lon_bc >= domain_lon_range(1) & data1_lon_bc <= domain_lon_range(2));
%     ind_lat      = find(data1_lat_bc >= domain_lat_range(1) & data1_lat_bc <= domain_lat_range(2));
%     data1_bc     = data1_bc(ind_lon,ind_lat);
%     data1_lon_bc = data1_lon_bc(ind_lon);
%     data1_lat_bc = data1_lat_bc(ind_lat);
%     
%     ind_lon      = find(data2_lon_bc >= domain_lon_range(1) & data2_lon_bc <= domain_lon_range(2));
%     ind_lat      = find(data2_lat_bc >= domain_lat_range(1) & data2_lat_bc <= domain_lat_range(2));
%     data2_bc     = data2_bc(ind_lon,ind_lat);
%     data2_lon_bc = data2_lon_bc(ind_lon);
%     data2_lat_bc = data2_lat_bc(ind_lat);
    
    %/ If data size is not consistent, do interpolation to data1_bc's gridding if not
    if ~isequal(size(data1_bc), size(data2_bc))    
        warning('Sizes of data1_bc and data2_bc are not consistent! Interpolating data2_bc to data1_bc''s gridding...');  
        if isempty(interp_method)  
            warning('''interp_method'' not specified! Assume it to be ''linear''. CAVEAT:''nearest'' will make up non-NaN values at grid cells where it shouldn''t');
            interp_method = 'linear'; 
        end
        
        %/ griddata() seems can handle unregular gridding, though I have
        %/ written my_interp() myself to handle it...
%         data2_bc     = griddata(data2_lon_bc_2D, data2_lat_bc_2D, data2_bc', data1_lon_bc_2D, data1_lat_bc_2D, interp_method)';
        

        disp(size(data2_lon_bc_2D))
        disp(size(data2_lat_bc_2D))
        disp(size(data2_bc))
        disp(size(data1_lon_bc_2D))
        disp(size(data1_lat_bc_2D))

        data2_bc     = my_interp('lon_old', data2_lon_bc_2D, 'lat_old', data2_lat_bc_2D, 'data', data2_bc,...
                                 'lon_new', data1_lon_bc_2D, 'lat_new', data1_lat_bc_2D, 'lon_dim', lon_dim, 'interp_method', interp_method);

        data2_lon_bc = data1_lon_bc;
        data2_lat_bc = data1_lat_bc;
    end
    
    %/ After interpolation, we subset the domain as queried
    ind_lon      = find(data1_lon_bc >= domain_lon_range(1) & data1_lon_bc <= domain_lon_range(2));
    ind_lat      = find(data1_lat_bc >= domain_lat_range(1) & data1_lat_bc <= domain_lat_range(2));
    data1_bc     = data1_bc(ind_lon,ind_lat);
    data1_lon_bc = data1_lon_bc(ind_lon);
    data1_lat_bc = data1_lat_bc(ind_lat);
    
    ind_lon      = find(data2_lon_bc >= domain_lon_range(1) & data2_lon_bc <= domain_lon_range(2));
    ind_lat      = find(data2_lat_bc >= domain_lat_range(1) & data2_lat_bc <= domain_lat_range(2));
    data2_bc     = data2_bc(ind_lon,ind_lat);
    data2_lon_bc = data2_lon_bc(ind_lon);
    data2_lat_bc = data2_lat_bc(ind_lat);
    
    if ~isempty(land_or_ocean)
        data1_bc = show_land_or_ocean_hydrosheds('matrix', data1_bc, 'lon_grids', data1_lon_bc, 'lat_grids', data1_lat_bc, 'land_or_ocean', land_or_ocean);
        data2_bc = show_land_or_ocean_hydrosheds('matrix', data2_bc, 'lon_grids', data2_lon_bc, 'lat_grids', data2_lat_bc, 'land_or_ocean', land_or_ocean);
    end

    %/ If show_TP_only, we compute SSIM, r and RMSE based on TP grid cells only.
    if show_TP_only
        %/ [NOTE] We *need* to reload TP, since length(contf_lat) is 181 for
        %/        reanalysis, but 179 for WSV data
        [reg_2D, ~, ~, ~, ~] = reg_extractor('lon', data1_lon_bc, 'lat', data1_lat_bc, 'slct_reg', 'TP',...
                                             'data_folder', data_folder, 'saveload_cond_landocean', 0, 'savemat', 1, 'recompute', 0);
        reg_2D(isnan(reg_2D)) = 0;
        cond = logical(reg_2D);
        data1_bc(~cond) = nan;
        data2_bc(~cond) = nan;
    end
              
    %/ Compute similarity, pattern correlation, RMSE
    %/ (based on the grid cells that both data show non-NaN values)
    cond_1 = ~isnan(data1_bc);
    cond_2 = ~isnan(data2_bc);
    cond   = (cond_1 & cond_2);
    
    data1_noNaN = data1_bc(cond);  
    data2_noNaN = data2_bc(cond);
%     length(data1_noNaN)
%     length(data2_noNaN)
    if length(data1_noNaN) ~= length(data2_noNaN)
        error('Check data availbility!');
    end          
                        
    SSIM        = ssim(data1_noNaN, data2_noNaN);
    [r, pval]   = corr(data1_noNaN, data2_noNaN);   %/ direct correlation (uncentered)
    RMSE        = sqrt(immse(data1_noNaN, data2_noNaN));
    
    if show_prct_diff
        diff_data = (data1_bc - data2_bc)./data2_bc*100;   %/ [%]
        diff_data(isinf(diff_data)) = nan;        %/ Convert inf (if any) into nan (since data2_bc can be 0).
    else
        diff_data = data1_bc - data2_bc;
    end
    
    diff_data_lon = data1_lon_bc;
    diff_data_lat = data1_lat_bc;
    
    fprintf('*** [diff_2D_landmask]: SSIM = %.2f, PCC = %.2f (p = %2.f), RMSE = %.2f ***\n', SSIM, r, pval, RMSE);

end