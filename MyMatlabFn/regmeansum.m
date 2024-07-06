function [regional_data, unit, reg_name_list] = regmeansum(varargin)

    %/ create a set of valid parameters and their default value
    pnames = {'slct_data', 'data', 'area_mean_or_sum', 'lon', 'lat', 'dates', 'slct_reg', 'DivByAreaOf', 'mth', 'regime_ver', 'data_folder', }; 
    dflts  = cell(length(pnames), 1);
    %/ parse function arguments
    [           slct_data,   data,   area_mean_or_sum,   lon,   lat,   dates,   slct_reg,   DivByAreaOf,   mth,   regime_ver,   data_folder,  ] = internal.stats.parseArgs(pnames, dflts, varargin{:});

%%
    %======================================================================
    %/ Author: Franklin Cheng
    %/ Last update: 23 May 2024
    %/
    %/ Description:
    %/             'slct_reg': The region in which the water budget is summed or averaged
    %/
    %/       'DivByAreaOf': Normalize the total water budget (only when 'area_mean_or_sum' = 'sum')
    %/                         by dividing it by the target area. 
    %/                         Useful when you are dealing with "source contribution."
    %/
    %======================================================================
    if ~(isequal(area_mean_or_sum, 'mean') || isequal(area_mean_or_sum, 'sum'))
        error('Ivalid ''area_mean_or_sum''!')
    end
    Area_global = calc_grid_area('lon', lon, 'lat', lat);
    nlon = length(lon);
    nlat = length(lat);
    
    %/ [IMPORTANT] Add the new names to the list of precipitation and evaporation products
    ListOfPE = {'Pm', 'Pm_BL', 'Pm_frac_real', 'Cf_map', 'Cf_map_frac_real', 'P_LA',...
                'pr', 'evspsbl', 'CERA_P', 'ERA5_P', 'ERA5_P_05', 'CMORPH_P', 'IMERG_P', 'GPCC_P', 'GPCP_P', 'CRU_P', 'HARv2_P', 'TPR_P', 'EM_P',...
                'CERA_E', 'ERA5_E', 'GLEAM_E', 'HARv2_E', 'TPR_E', 'EM_E', 'cmip6_pr', 'cmip6_evspsbl'};

    %/ Get the area of the reference region ('DivByAreaOf')
    if ~isempty(DivByAreaOf)      
        [reg_2D_ref, bndry_data_ref, ~, ~, ~] = reg_extractor('slct_reg', DivByAreaOf, 'lon', lon, 'lat', lat,...
                                                   'mth', mth, 'regime_ver', regime_ver, 'data_folder', data_folder, 'savemat', 0, 'recompute', 0,...
                                                   'verbose', 0);
        logi_2D_ref                   = reg_2D_ref;
        % logi_2D_ref(logi_2D_ref ~= 1) = 0;
        logi_2D_ref(isnan(logi_2D_ref)) = 0;
        logi_2D_ref                   = logical(logi_2D_ref);
        Area_ref                      = Area_global;
        Area_ref(~logi_2D_ref)        = 0;  
        Area_ref_1D                   = reshape(Area_ref, [], 1);
        
        if isequal(area_mean_or_sum, 'sum')
            fprintf('*** [regmeansum]: Detected non-empty ''DivByAreaOf''. Normalizing the summed water budget by the target area... ***\n')
        elseif isequal(area_mean_or_sum, 'mean')
            warning('DivByAreaOf is requested, but ''area_mean_or_sum'' = ''mean'', not ''sum''. Hence no normalization will be performed.')
        end
    end
    
    %/ Obtain reg_2D (based on the data lon, lat)
    [reg_2D, ~, reg_name_list, ~, ~] = reg_extractor('slct_reg', slct_reg, 'lon', lon, 'lat', lat,...
                                                     'mth', mth, 'regime_ver', regime_ver, 'data_folder', data_folder, 'savemat', 0, 'recompute', 0,...
                                                     'verbose', 0);

    Ndate           = length(dates);
    regional_data   = nan(length(reg_name_list), Ndate);

    %/ Loop the reg_name_list to compute the quantitfy for the source area!
    for i = 1:length(reg_name_list)
        logi_src_2D                   = reg_2D;
        logi_src_2D(logi_src_2D ~= i) = 0;
        logi_src_2D                   = logical(logi_src_2D);

        %/ Reshape to 1D is much *faster* than making up a 3D logical matrix!
        logi_src_1D            = reshape(logi_src_2D, [], 1);    %/ Grids x Dates
        data_bc                = data;
        data_bc                = reshape(data_bc, nlon*nlat, Ndate);   %/ Grids x Dates
        Area_src               = Area_global;
        Area_src(~logi_src_2D) = 0;  
        Area_src_1D            = reshape(Area_src, [], 1); %/ Grids x 1

        if isequal(area_mean_or_sum, 'sum')
            if ismember(slct_data, {'Pm_frac', 'Pm_frac_adj'}) || isequal(slct_data, 'CMF') 
            % if contains(slct_data, {'Pm_frac'}) || isequal(slct_data, 'CMF')  %/ Old line! Previous calculation of Pm_frac_real did not consider area-weighting when amassing moisture contributions

                regional_data(i, :) = sum(data_bc(logi_src_1D,:), 1, 'omitnan'); %/ They are ratio contributions, sum them up.
                
                if ismember(slct_data, {'CMF'})
                    unit = ' ';
                else
                    unit = '%';
                end
            else
                %/ [IMPORTANT] Compute the volume/mass before amassing the moisture contributions
                data_bc = data_bc.*Area_src_1D;  %/ Grids x Dates, [unit] -> [unit] m^2

                if ismember(slct_data, ListOfPE)
                    regional_data(i, :) = sum(data_bc(logi_src_1D,:), 1, 'omitnan')*1e-12;   %/ 1 x Dates, km^3/day
                    unit = 'km^{3} day^{-1)';
                    
                    %/ If to divide the volume by the area of the reference sink region  
                    if ~isempty(DivByAreaOf)  
                        regional_data(i, :) = regional_data(i, :)./(sum(Area_ref_1D)*1e-12); %/ 1 x Dates, mm/day
                        unit = 'mm day^{-1}';  
                    end
                else
                    error('It seems wrong to set area_mean_or_sum == ''%s'' for slct_data == ''%s''!', area_mean_or_sum, slct_data)
                end
            end

        elseif isequal(area_mean_or_sum, 'mean')
            %/ Do area-weighted mean based on grid area
            %/ [CAVEAT]: Let say if you want to compute the "areal mean"
            %/           contribution from westerly regime, it will weigh all the
            %/           westerly grids globally -> this is why the value can be very small!

            % size(data_bc)
            % size(Area_regional_1D)

            data_bc             = (data_bc.*Area_src_1D)./sum(Area_src_1D, 'omitnan');   %/ multiply the weights
            regional_data(i, :) = sum(data_bc(logi_src_1D,:), 1, 'omitnan');                        %/ summing up the weighted terms
            
            % slct_data
            if contains(slct_data, {'T2m'}) || contains(slct_data, {'tas'}) 
                unit = 'K';  
            elseif contains(slct_data, {'SM_layer'}) || ismember(slct_data, {'GLEAM_SMsurf'})
                unit = 'm^{3} m^{-3} day^{-1}';
            elseif contains(slct_data, {'IVT'})
                unit = 'kg/m/s';   
            elseif ismember(slct_data, ListOfPE)
                unit = 'mm/day';
            elseif contains(slct_data, {'ua', 'va'})
                unit = 'm/s';
            else
                unit = '(?)';  %/ by default (may be wrong)
            end
        end
    end
end