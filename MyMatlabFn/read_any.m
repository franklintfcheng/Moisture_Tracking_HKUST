function dataset = read_any(varargin)

    %/ Create a set of valid parameters and their default value
    pnames = {'dataset',  'param', 'select_field',  'slct_data',  'slct_year',...
              'time_dim',  'noleap', 'stlevel', 'edlevel', 'slct_level', 'lon_range', 'lat_range',...
              'data_folder', 'savemat', 'recompute', 'set_conden_to_zero', 'NumWorkers',...
              'load_or_not', 'plot_to_check', 'savefig'};

    dflts  = cell(length(pnames), 1);

    %/ Parse function arguments
             [ dataset,      param,  select_field,  slct_data, slct_year,...
               time_dim,   noleap, stlevel,   edlevel,   slct_level, lon_range,   lat_range,...
               data_folder, savemat,   recompute, set_conden_to_zero,   NumWorkers,...
               load_or_not,  plot_to_check,  savefig] ...
                   = internal.stats.parseArgs(pnames, dflts, varargin{:});
%%
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: May 22, 2024
    %/
    %/ This function is designed to read any kind of dataset, 
    %/   including
    %/            - Reanalysis (ERA5, ERA5-Land, CERA-20C)
    %/            - Satellite products
    %/            - Multi-product ensemble mean
    %/=====================================================================

    %/ Load the parameter module
    run(param);

    if ischar(select_field)
        select_field = {select_field};
    end

    %/ By default, load global data (except when lon_range or lat_range is given)
    stlon =   0;  edlon = 360;
    stlat = -90;  edlat = 90;
    if ~isempty(lon_range) 
        stlon = lon_range(1); edlon = lon_range(end);
    end
    if ~isempty(lat_range)
        stlat = lat_range(1); edlat = lat_range(end);
    end
    if isequal(stlon, 0) && isequal(edlon, 360) && isequal(stlat, -90) && isequal(edlat, 90)
        str_domain = '_global';
    else
        str_domain = sprintf('_%g-%gE_%g-%gN', stlon, edlon, stlat, edlat);  %/ '%g' removes trailing zeros
    end

    if ischar(slct_data)
        slct_data = {slct_data};
    end
    
    select_data = findismember_loop(dataname, slct_data); 
    if isempty(select_data)    
        error('The input data is not found in dataname. Update %s!', param);     
    end

    if set_conden_to_zero      
        str_conden = [];    
    else     
        str_conden = '_with_conden';  
    end
    
    %==================== set years and case_date ========================%
    seq = sequence('data', slct_year);
    str_years_bc = '';
    for i = 1:length(seq)
        a = strjoin({num2str(seq{i}(1)), num2str(seq{i}(end))}, '-');
        str_years_bc = strjoin({str_years_bc, a}, '_');
    end
    case_date = date_array_gen('year_list', slct_year,  'output_date_format', 'yyyymmdd');  
    %=====================================================================%

    for k = select_data
        disp(dataname{k})
        for i = 1:length(select_field)
            tic
            if contains(select_field{i}, 'monthly')
                flds_date          = 'date_yyyymm';

            elseif contains(select_field{i}, 'daily')
                flds_date          = 'date_yyyymmdd_AllYr';

            elseif contains(select_field{i}, 'subdaily')
                flds_date          = 'date_UTC_AllYr';
            else
                error('code not set!'); 
            end

            mainfield_list = {select_field{i}, 'lon', 'lat', 'date_yyyymmdd_AllYr', flds_date, 'units'};
            mainfield_list = unique(mainfield_list, 'stable');

            if ismember(data_path_parts{whichfolder(k),2}, {'ERA5-Land_'}) && ismember(datatype{k}, {'acc'}) 
                timeshift = -1; %/ [day], time = 00:00 means daily accumulation of the previous day in ERA5-Land.
            else
                timeshift = 0;
            end

            if noleap 
                if contains(mainfield_list{1}, 'daily')
                    str_noleap = '_noleap';
                elseif contains(mainfield_list{1}, 'pentad')
                    noleap     = 0;
                    str_noleap = '';
                    warning('noleap is automatically turned off for %s data', mainfield_list{1})
                elseif contains(mainfield_list{1}, 'monthly')
                    noleap     = 0;
                    str_noleap = '';
                    warning('noleap is automatically turned off for %s data', mainfield_list{1})
                else
                    error('noleap is not set for %s data!', mainfield_list{1});
                end
            else
                str_noleap = '';
            end

            str_level_range = '';
            if ~isempty(slct_level)
                dataname_bc = sprintf('%s%d',dataname{k}, slct_level); %/ UPDATE; e.g., q -> q1000 (slct_level == 1000)
            else
                dataname_bc = dataname{k};
                if ~isempty(stlevel) && ~isempty(edlevel)
                    if stlevel > edlevel
                        str_level_range = sprintf('_lv%d-%d', edlevel, stlevel);
                    else
                        str_level_range = sprintf('_lv%d-%d', stlevel, edlevel);
                    end
                    mainfield_list{end+1} = 'level';  %/ Query the level field
                end
            end

            %/ Check if file exists
            if ismember(dataname_bc, {'HiVeg','LoVeg'})   %/ invariants
                filename = strcat(data_folder, dataname_bc, '_', mainfield_list{1}, str_domain, str_conden, '.mat');
            else
                filename = strcat(data_folder, dataname_bc,str_level_range, '_', mainfield_list{1}, str_domain, str_years_bc, str_noleap, str_conden, '.mat');
            end
            
            %/ Data Processing
            if isfile(filename) && recompute == 0
                if load_or_not
                    fprintf('*** Loading %s ***\n', filename);
                    load(filename, 'S');
                    fld = fieldnames(S);
                    for f = 1:length(fld)
                         dataset.(dataname_bc).(fld{f}) = S.(fld{f});
                    end
                else
                    fprintf('!!! File exists. Skip loading as per requested. !!!\n');
                end
            else
                p_Wang1988   = [300, 500, 700];                      %/ The two-layer model of Wang 1988 
                p2           = p_Wang1988(2)*100;                    %/ hPa to Pa
                dp           = (p_Wang1988(3) - p_Wang1988(1))*100;  %/ hPa to Pa
                parts        = strsplit(dataname_bc, '_');
                whichdataset = strcat(parts{1}, '_');                %/ e.g., 'ERA5_'

                dataset = []; dataset.placeholder = [];
                if contains(dataname_bc, {'GLEAM', 'CMORPH', 'IMERG', 'HARv2', 'TPR', 'GPCC', 'HadCRUT5', 'CRU', 'GPCP', 'GLASS', 'NOAA_CDR', 'NOAA_Interp'}) % Read satellite data
                    if isequal(dataname_bc, 'TPR_RO')
                        filename_SRO = strcat(data_folder, 'TPR_SRO', '_', mainfield_list{1}, str_domain, str_years_bc, str_noleap, str_conden, '.mat');
                        fprintf('*** Loading %s... ***\n', filename_SRO);
                        load(filename_SRO, 'S');
                        dataset.('TPR_SRO').((mainfield_list{1})) = S.(mainfield_list{1});
                        clear S;
                        
                        filename_SSRO = strcat(data_folder, 'TPR_SSRO', '_', mainfield_list{1}, str_domain, str_years_bc, str_noleap, str_conden, '.mat');
                        fprintf('*** Loading %s... ***\n', filename_SSRO);
                        load(filename_SSRO, 'S');
                        dataset.('TPR_SSRO').((mainfield_list{1})) = S.(mainfield_list{1});
    
                        dataset.(dataname_bc).(mainfield_list{1}) = dataset.('TPR_SRO').((mainfield_list{1})) ...
                                                      + dataset.('TPR_SSRO').((mainfield_list{1})); %/ Add up to total runoff
                        dataset.(dataname_bc).lon = S.lon;
                        dataset.(dataname_bc).lat = S.lat;
                        dataset.(dataname_bc).date_yyyymmdd_AllYr = S.date_yyyymmdd_AllYr;
                        clear S;
                    else
                        sate_product = dataname_bc;
                        dataset = read_satedata('sate_product', sate_product, 'slct_year', slct_year, 'select_field', mainfield_list{1},...
                                                'lon_range', lon_range, 'lat_range', lat_range, 'NumWorkers', NumWorkers); %/ Load CMORPH prcp
                    end
                    
                elseif ismember(dataname_bc, {'EM_P', 'EM_E'}) %/ Multi-Product Ensemble Mean
                    if ismember(select_field, {'daily'})
                        error('Only monthly ensemble mean is permitted now!');
                        
                    elseif ismember(select_field, {'monthly'})
                        output_dates = date_array_gen('year_list', slct_year,  'dt_slct_mo', 1, 'output_date_format', 'yyyyMM');  
                        if isequal(dataname_bc, 'EM_P')
                            products = {'ERA5_P', 'GPCC_P', 'GPCP_P', 'CRU_P', 'HARv2_P', 'TPR_P'};
                        elseif isequal(dataname_bc, 'EM_E')
                            products = {'ERA5_E', 'GLEAM_E', 'HARv2_E', 'TPR_E',};
                        else
                            error('code not set!');
                        end
                    else
                        error('code not set!');
                    end
                    lon_new = (0:359)'; lat_new = (90:-1:-90)';  %/ 1x1
                    counts  = zeros(length(lon_new), length(lat_new), length(output_dates));
                    EM_data = zeros(length(lon_new), length(lat_new), length(output_dates));
                    for ii = 1:length(products)
                        filename_P = strcat(data_folder, products{ii}, '_', mainfield_list{1}, str_domain, str_years_bc, str_noleap, str_conden, '.mat');
                        fprintf('*** Loading %s... ***\n', filename_P);
                        load(filename_P, 'S');
                        lon_old = S.lon;
                        lat_old = S.lat;
                        units   = S.units;

                        %/ Standardize lon to degree east only
                        lon_old(lon_old < 0) = lon_old(lon_old < 0) + 360;

    %                     if contains(products{ii}, {'TPR', 'HARv2'})
    %                         interp_method = 'grid2grid';
    %                     else
                        interp_method = 'nearest';
    %                     end
                        %/ Interpolate into 1x1 gridding
                        NumWorkers_interp = 40;
                        [S.(mainfield_list{1}), ~, ~] = my_interp('lon_old', lon_old, 'lat_old', lat_old, 'data', S.(mainfield_list{1}),...
                                                                  'lon_new', lon_new, 'lat_new', lat_new, 'interp_method', interp_method,...
                                                                  'NumWorkers', NumWorkers_interp);
                        
                        %/ [IMPORTANT] Subset the TP grids only for those regionally downscaled products (HARv2, TPR)
                        if contains(products{ii}, {'TPR', 'HARv2'})
                            [cond, ~, ~, ~, ~] = reg_extractor('slct_reg', 'TP', 'lon', lon_new, 'lat', lat_new,...
                                                               'data_folder', data_folder, 'savemat', 0, 'recompute', recompute);
                            cond(isnan(cond)) = 0;
                            cond = logical(cond);
                            cond = repmat(cond, 1, 1, size(S.(mainfield_list{1}),3));
                            S.(mainfield_list{1})(~cond) = NaN;
                        end
                        
                        %/ Plot to check
                        if plot_to_check
                            contf_data = mean(S.(mainfield_list{1}), 3, 'omitnan');
                            contf_lon = lon_new;
                            contf_lat = lat_new;
                            if ismember(select_field, {'monthly'})
                                contf_levels = (0:11)*30;
                            else
                                contf_levels = 0:11; 
                            end
                            cbar_interval = 1;
                            colmap        = my_colormap(length(contf_levels)-1, 'precip3_11lev');
                            colmap(1,:)   = [.7 .7 .7];
                            color         = [ 0 0 0];
                            fontsize      = 22;
                            grid_mode     = 3;
                            coast_col     = [.1 .1 .1];
                            coast_wi      = 2.5;
                            % bndry_data    = [basin_catalog.bndry, TP_bndry];
                            bndry_data = [];
                            % map_lat_lower = -6;
                            % map_lat_upper = 46; 
                            % map_lon_lower = 26;
                            % map_lon_upper = map_lon_lower + (107-64)/(44-25)*(map_lat_upper-map_lat_lower);
                            map_lat_lower = -90;
                            map_lat_upper = 90; 
                            map_lon_lower = -179;
                            map_lon_upper = 180;
                            
                            titlename     = strcat({'mean '}, strrep(products{ii}, '_', ' '), {' interp 1x1'});
                            if savefig  savepath = strcat(plotting_folder, titlename);  else  savepath = [];  end
                            plot_contfmap('contf_data', contf_data, 'contf_lon', contf_lon, 'contf_lat', contf_lat,...
                                          'contf_levels', contf_levels, 'cbar_interval', cbar_interval, 'colmap', colmap, 'pcolor_mode', 1,...
                                          'linewi', 3, 'color', color, 'bndry_data', bndry_data, 'titlename', titlename, 'savepath', savepath, 'fig_fmt', 'png',...
                                          'map_lon_lower', map_lon_lower, 'map_lon_upper', map_lon_upper, 'map_lat_lower', map_lat_lower, 'map_lat_upper', map_lat_upper,...
                                          'coast_col', coast_col, 'coast_wi', coast_wi, 'fontsize', fontsize, 'grid_mode', grid_mode, ...
                                          'cbar_mode', 1, 'cbar_position', 'eastoutside');
                        end
    
                        %/ Get the correct indices of the dates 
                        ind_date = findismember_loop(output_dates, S.date_yyyymm);
                        if isempty(ind_date)  
                            error('Empty ind_date! Check your code!'); 
                        end
                        
                        %/ Replace nan with zeros (for reduction)
                        cond_zero = (S.(mainfield_list{1}) == 0);
                        cond_nan  = isnan(S.(mainfield_list{1}));
                        S.(mainfield_list{1})(cond_nan) = 0;
                        counts_each = logical(S.(mainfield_list{1}));
                        counts_each(cond_zero) = 1;  %/ count the *real* zeros
                        
                        %/ Reduction
                        counts(:,:,ind_date)  = counts(:,:,ind_date)  + double(counts_each);  
                        EM_data(:,:,ind_date) = EM_data(:,:,ind_date) + S.(mainfield_list{1});                  
                        S = [];
                    end
                    EM_data = EM_data./counts;  %/ Get the ensemble mean (will contain NaN if no product in that index)
                    dataset.(dataname_bc).(mainfield_list{1}) = EM_data;
                    dataset.(dataname_bc).lon                 = lon_new;
                    dataset.(dataname_bc).lat                 = lat_new;
                    dataset.(dataname_bc).(flds_date)         = output_dates;
                    dataset.(dataname_bc).units               = units;
                    
                    EM_data = []; %/ spare memory
                    
                    %/ Plot to check
                    if plot_to_check
                        contf_data = mean(dataset.(dataname_bc).(mainfield_list{1}), 3, 'omitnan');
                        contf_lon = lon_new;
                        contf_lat = lat_new;
                        if ismember(select_field, {'monthly'})
                            contf_levels = [0:11]*30;
                        else
                            contf_levels = 0:11; 
                        end
                        colmap = my_colormap(length(contf_levels)-1, 'precip3_11lev');
                        colmap(1,:) = [.7 .7 .7];
                        color = [ 0 0 0];
                        fontsize      = 22;
                        grid_mode     = 3;
                        coast_col     = [.1 .1 .1];
                        coast_wi      = 2.5;
                        % bndry_data = [basin_catalog.bndry, TP_bndry];
                        bndry_data = [];
                        % map_lat_lower = -6;
                        % map_lat_upper = 46; 
                        % map_lon_lower = 26;
                        % map_lon_upper = map_lon_lower + (107-64)/(44-25)*(map_lat_upper-map_lat_lower);
                        map_lat_lower = -90;
                        map_lat_upper = 90; 
                        map_lon_lower = -179;
                        map_lon_upper = 180;

                        titlename     = strcat({'mean '}, strrep(dataname_bc, '_', ' '), {' '}, select_field, {' interp 1x1'});
                        if savefig  savepath = strcat(plotting_folder, titlename);  else  savepath = [];  end
                        plot_contfmap('contf_data', contf_data, 'contf_lon', contf_lon, 'contf_lat', contf_lat,...
                                      'contf_levels', contf_levels, 'colmap', colmap, 'pcolor_mode', 1,...
                                      'linewi', 3, 'color', color, 'bndry_data', bndry_data, 'titlename', titlename, 'savepath', savepath, 'fig_fmt', 'png',...
                                      'map_lon_lower', map_lon_lower, 'map_lon_upper', map_lon_upper, 'map_lat_lower', map_lat_lower, 'map_lat_upper', map_lat_upper,...
                                      'coast_col', coast_col, 'coast_wi', coast_wi, 'fontsize', fontsize, 'grid_mode', grid_mode,...
                                      'cbar_mode', 1, 'cbar_position', 'eastoutside');
                    end
                    
                elseif ismember(dataname_bc, {'HiVeg','LoVeg'})  %/ Read vegetation type (CERA-20C)
                    data_file = sprintf('/disk/r128/tfchengac/fandyr128SOM/cera20c_data_Annual/cera20c_%s.nc', dataname_bc);
                    ncdisp(data_file);
                    data = ncread(data_file, varname{k});
                    data = int64(round(data(:,:,1),0)); %/ retrieve the 1st ensemble. Veg Type is the same for all ensemble.
    
                    dataset.(dataname_bc).(select_field{i}) = data;
                    dataset.(dataname_bc).lon = double(ncread(data_file, 'longitude'));
                    dataset.(dataname_bc).lat = double(ncread(data_file, 'latitude'));
                    dataset.(dataname_bc).typeindex = unique(reshape(data, [], 1));
                    if isequal(dataname_bc, 'HiVeg')
                        dataset.(dataname_bc).typename = {'NoVeg', 'Evergreen needleleaf trees', 'Deciduous needleleaf trees',...
                            'Deciduous broadleaf trees', 'Evergreen broadleaf trees', 'Mixed forest/woodland', 'Interrupted forest'};
                    elseif isequal(dataname_bc, 'LoVeg')
                        dataset.(dataname_bc).typename = {'NoVeg', 'Crops', 'Grass', 'Tall grass', 'Tundra', 'Irrigated crops',...
                                                           'Semidesert', 'Bogs and marshes', 'Evergeen shrubs', 'Deciduous shrubs'};
                        
                    end
                    
                elseif contains(dataname_bc, {'E_minus_P'})
                    E_filename = strcat(data_folder, whichdataset, 'E_', select_field{i}, str_domain, str_years_bc, str_noleap, '.mat');
                    load(E_filename, 'S')
                    E = S.(select_field{i});
    
                    P_filename = strcat(data_folder, whichdataset, 'P_', select_field{i}, str_domain, str_years_bc, str_noleap, '.mat');
                    load(P_filename, 'S')
                    P = S.(select_field{i});
    
                    dataset.(dataname_bc).(select_field{i}) = E - P;
                    dataset.(dataname_bc).lon = S.lon;
                    dataset.(dataname_bc).lat = S.lat;
                    dataset.(dataname_bc).date_yyyymmdd_AllYr = S.date_yyyymmdd_AllYr;   
    
                elseif contains(dataname_bc, {'dTWS'})
                    P_filename = strcat(data_folder, whichdataset, 'P_', select_field{i}, str_domain, str_years_bc, str_noleap, '.mat');
                    load(P_filename, 'S')
                    P = S.(select_field{i});

                    E_filename = strcat(data_folder, whichdataset, 'E_', select_field{i}, str_domain, str_years_bc, str_noleap, '.mat');
                    load(E_filename, 'S')
                    E = S.(select_field{i});

                    R_filename = strcat(data_folder, whichdataset, 'SRO_', select_field{i}, str_domain, str_years_bc, str_noleap, '.mat');
                    load(R_filename, 'S')
                    R = S.(select_field{i});
    
                    dataset.(dataname_bc).(select_field{i}) = P - E - R;
                    dataset.(dataname_bc).lon = S.lon;
                    dataset.(dataname_bc).lat = S.lat;
                    dataset.(dataname_bc).date_yyyymmdd_AllYr = S.date_yyyymmdd_AllYr;   

                elseif contains(dataname_bc, {'S500', 'S2_Wang1988'})
                    filename_nest = strcat(data_folder, whichdataset, 'T300_', select_field{i}, str_domain, str_years_bc, str_noleap, '.mat');
                    load(filename_nest, 'S');
                    T300 = S.(select_field{i});

                    filename_nest = strcat(data_folder, whichdataset, 'T500_', select_field{i}, str_domain, str_years_bc, str_noleap, '.mat');
                    load(filename_nest, 'S');
                    T500 = S.(select_field{i});

                    filename_nest = strcat(data_folder, whichdataset, 'T700_', select_field{i}, str_domain, str_years_bc, str_noleap, '.mat');
                    load(filename_nest, 'S');
                    T700 = S.(select_field{i});

                    p_unit = 'hPa';
                    PT = [];
                    PT.p300 = compute_PT('T', T300, 'p', 300, 'p_unit', p_unit);
                    PT.p500 = compute_PT('T', T500, 'p', 500, 'p_unit', p_unit);
                    PT.p700 = compute_PT('T', T700, 'p', 700, 'p_unit', p_unit);

                    if contains(dataname_bc, {'S500'})
                        dataset.(dataname_bc).(select_field{i}) = (T500./PT.p500).*(PT.p300 - PT.p700)/(400*100);
                        dataset.(dataname_bc).unit = 'K Pa^{-1}';  

                    elseif contains(dataname_bc, {'S2_Wang1988'})
                        %/ Compute dry static stability at 500 hPa using Eq. (3.2) in Wang (1988), 
                        %/ useful when computing the nondimensional I.

                        %/ Assume for dry air, the air density can be
                        %/ obtained from ideal gas law (rho = P/RT)
                        R      = 287;                %/ J kg-1 K-1
                        rho500 = 500*100./(R.*T500); %/ kg m^-3
                        
                        dataset.(dataname_bc).(select_field{i}) = (1./PT.p500./rho500).*(PT.p300 - PT.p700)/dp;
                        dataset.(dataname_bc).unit = 'm^{2} s^{2} Pa^{-2}';  
                        clear rho500;
                    end

                    %/ Store the basic info and clear the data
                    dataset.(dataname_bc).lon = S.lon;
                    dataset.(dataname_bc).lat = S.lat;
                    if contains(select_field{i}, 'monthly')
                        dataset.(dataname_bc).date_yyyymm = S.date_yyyymm;   
                    else
                        dataset.(dataname_bc).date_yyyymmdd_AllYr = S.date_yyyymmdd_AllYr;   
                    end

                    %/ Spare memory
                    clear S T300 T500 T700

                elseif contains(dataname_bc, {'C0_Wang1988', 'alpha_Wang1988'}) 
                    % slct_data_nest = strcat(whichdataset, {'S2_Wang1988'});
                    % filename_nest = strcat(data_folder, slct_data_nest, '_', select_field{i}, str_domain, str_years_bc, str_noleap, '.mat');
                    % load(filename_nest, 'S');
                    % S2 = S.(select_field{i});

                    slct_data_nest = strcat(whichdataset, {'S2_Wang1988'});
                    dataset_temp = []; 
                    stlevel_nest = []; 
                    edlevel_nest = [];
                    for ii = 1:length(slct_data_nest)
                        dataset_temp = read_any('dataset',  dataset_temp, 'param', param, 'select_field', select_field{i}, 'slct_data', slct_data_nest{ii}, 'slct_year', slct_year,...
                                                'time_dim', 4, 'noleap', noleap, 'stlevel', stlevel_nest, 'edlevel', edlevel_nest, 'slct_level', [], 'lon_range', lon_range, 'lat_range', lat_range,...
                                                'data_folder', data_folder, 'savemat', savemat, 'recompute', 0, 'set_conden_to_zero', set_conden_to_zero, 'NumWorkers', NumWorkers,...
                                                'load_or_not', 1, 'plot_to_check', 0, 'savefig', 0);
                    end
                    S2    = dataset_temp.(slct_data_nest{1}).(select_field{i});
                    C0    = sqrt(1/2*S2*dp^2);             %/ m s-1

                    Cp    = 1004;                          %/ J kg-1 K-1
                    R     = 287;                           %/ J kg-1 K-1
                    Lc    = 2.26e6;                        %/ J kg-1
                    b     = 0.9;                           %/ non-dimensional
                    alpha = (2*Cp*p2*C0.^2)/(R*b*Lc*dp);   %/ non-dimensional
                    
                    if contains(dataname_bc, {'C0_Wang1988'})
                        dataset.(dataname_bc).(select_field{i}) = C0;
                        dataset.(dataname_bc).unit = 'm s^{-1}';  

                    elseif contains(dataname_bc, {'alpha_Wang1988'})
                        dataset.(dataname_bc).(select_field{i}) = alpha;
                        dataset.(dataname_bc).unit = '';  
                    end
                    
                    %/ Copy the rest of the fields
                    other_flds = {'lon', 'lat', flds_date};
                    for f = 1:length(other_flds)
                        dataset.(dataname_bc).(other_flds{f}) = dataset_temp.(slct_data_nest{1}).(other_flds{f});
                    end
                    clear dataset_temp S S2 C0 alpha

                elseif contains(dataname_bc, {'q3_bar_Wang1988', 'q1_bar_Wang1988'}) 
                    slct_data_nest = strcat(whichdataset, {'q'});

                    dataset_temp = [];
                    if contains(dataname_bc, {'q3_bar_Wang1988'}) 
                        stlevel_nest = 500;
                        edlevel_nest = 925;  %/ To be consistent with the available plevel of CMIP6 model output
                    elseif contains(dataname_bc, {'q1_bar_Wang1988'}) 
                        stlevel_nest = 100;
                        edlevel_nest = 500;
                    end
                    
                    for ii = 1:length(slct_data_nest)
                        dataset_temp = read_any('dataset',  dataset_temp, 'param', param, 'select_field', select_field{i}, 'slct_data', slct_data_nest{ii}, 'slct_year', slct_year,...
                                                'time_dim', 4, 'noleap', noleap, 'stlevel', stlevel_nest, 'edlevel', edlevel_nest, 'slct_level', [], 'lon_range', lon_range, 'lat_range', lat_range,...
                                                'data_folder', data_folder, 'savemat', savemat, 'recompute', 0, 'set_conden_to_zero', set_conden_to_zero, 'NumWorkers', NumWorkers,...
                                                'load_or_not', 1, 'plot_to_check', 0, 'savefig', 0);
                    end

                    level = dataset_temp.(slct_data_nest{1}).level;
                    disp(level)

                    take_vertmean = 1;  %/ Perform vertical mean
                    level_unit = 'Pa';  %/ The default unit of level in CMIP6 data 
                    
                    dataset.(dataname_bc).(select_field{i}) = vertinte('data', dataset_temp.(slct_data_nest{1}).(select_field{i}), 'level', level, 'level_unit', level_unit, 'take_vertmean', take_vertmean);
                    dataset.(dataname_bc).unit = 'kg/kg';  

                    other_flds = {'lon', 'lat', flds_date};
                    for f = 1:length(other_flds)
                        dataset.(dataname_bc).(other_flds{f}) = dataset_temp.(slct_data_nest{1}).(other_flds{f});
                    end
                    clear dataset_temp;

                elseif contains(dataname_bc, {'dq_bar_Wang1988'}) 
                    slct_data_nest = strcat(whichdataset, {'q3_bar_Wang1988', 'q1_bar_Wang1988'});
                    dataset_temp = []; 
                    stlevel_nest = []; 
                    edlevel_nest = [];
                    for ii = 1:length(slct_data_nest)
                        dataset_temp = read_any('dataset',  dataset_temp, 'param', param, 'select_field', select_field{i}, 'slct_data', slct_data_nest{ii}, 'slct_year', slct_year,...
                                                'time_dim', 4, 'noleap', noleap, 'stlevel', stlevel_nest, 'edlevel', edlevel_nest, 'slct_level', [], 'lon_range', lon_range, 'lat_range', lat_range,...
                                                'data_folder', data_folder, 'savemat', savemat, 'recompute', 0, 'set_conden_to_zero', set_conden_to_zero, 'NumWorkers', NumWorkers,...
                                                'load_or_not', 1, 'plot_to_check', 0, 'savefig', 0);
                    end
                    
                    dataset.(dataname_bc).(select_field{i}) = (dataset_temp.(slct_data_nest{1}).(select_field{i}) - dataset_temp.(slct_data_nest{2}).(select_field{i}));
                    dataset.(dataname_bc).unit = 'kg/kg';  

                    %/ Copy the rest of the fields
                    other_flds = {'lon', 'lat', flds_date};
                    for f = 1:length(other_flds)
                        dataset.(dataname_bc).(other_flds{f}) = dataset_temp.(slct_data_nest{1}).(other_flds{f});
                    end
                    clear dataset_temp;

                elseif contains(dataname_bc, {'I_Wang1988'}) 
                    slct_data_nest = strcat(whichdataset, {'q3_bar_Wang1988', 'q1_bar_Wang1988', 'alpha_Wang1988'});
                    dataset_temp = []; 
                    stlevel_nest = []; 
                    edlevel_nest = [];
                    for ii = 1:length(slct_data_nest)
                        dataset_temp = read_any('dataset',  dataset_temp, 'param', param, 'select_field', select_field{i}, 'slct_data', slct_data_nest{ii}, 'slct_year', slct_year,...
                                                'time_dim', 4, 'noleap', noleap, 'stlevel', stlevel_nest, 'edlevel', edlevel_nest, 'slct_level', [], 'lon_range', lon_range, 'lat_range', lat_range,...
                                                'data_folder', data_folder, 'savemat', savemat, 'recompute', 0, 'set_conden_to_zero', set_conden_to_zero, 'NumWorkers', NumWorkers,...
                                                'load_or_not', 1, 'plot_to_check', 0, 'savefig', 0);
                    end
                    
                    q3_bar = dataset_temp.(slct_data_nest{1}).(select_field{i});
                    q1_bar = dataset_temp.(slct_data_nest{2}).(select_field{i});
                    alpha  = dataset_temp.(slct_data_nest{3}).(select_field{i});

                    dataset.(dataname_bc).(select_field{i}) = (q3_bar - q1_bar)./alpha;
                    dataset.(dataname_bc).unit = '';  

                    %/ Copy the rest of the fields
                    other_flds = {'lon', 'lat', flds_date};
                    for f = 1:length(other_flds)
                        dataset.(dataname_bc).(other_flds{f}) = dataset_temp.(slct_data_nest{1}).(other_flds{f});
                    end
                    clear dataset_temp;

                elseif ismember(dataname_bc, {'ERA5_S'})  %/ 4D dry static stability (S)
                    %/ Load the 4D Temp data
                    filename_nest = strcat(data_folder, 'ERA5_T', str_level_range, '_', select_field{i}, str_domain, str_years_bc, str_noleap, '.mat');
                    fprintf('*** Loading %s... ***\n', filename_nest);
                    load(filename_nest, 'S');

                    %/ Store the basic info and clear the data
                    p = S.level;
                    T = S.(select_field{i});
                    dataset.(dataname_bc).lon   = S.lon;
                    dataset.(dataname_bc).lat   = S.lat;
                    dataset.(dataname_bc).level = S.level;
                    if contains(select_field{i}, 'monthly')
                        dataset.(dataname_bc).date_yyyymm = S.date_yyyymm;   
                    else
                        dataset.(dataname_bc).date_yyyymmdd_AllYr = S.date_yyyymmdd_AllYr;   
                    end
                    clear S;
                    
                    p_unit = 'hPa';  CONV = 100;
                    p_dim  = 3;
                    debug  = 1;
                    PT  = compute_PT('T', T, 'p_dim', p_dim, 'p', p, 'p_unit', p_unit, 'debug', debug);

                    %/ Static Stability
                    if p_dim == 3
                        [~, ~, dPT_dA, ~] = gradient(PT);                      %/ Chain rule to do derivate with *non-uniform* spacings.
                    else
                        error('code not set!');
                    end
                    dp_dA     = gradient(p * CONV);                        %/ get non-uniform gradient of pressure level first.
                    dp_dA     = reshape(dp_dA, [ones(1, p_dim-1), length(dp_dA)]);      %/ 1 x 1 x plev
                    dPT_dP    = dPT_dA./dp_dA;                         %/ J kg-1 Pa-1   elementwise division.  
                    SS        = -T./PT.*dPT_dP; 
                    
                    dataset.(dataname_bc).(select_field{i}) = SS;
                    dataset.(dataname_bc).unit = 'K Pa^{-1}';  
                    %/ Spare memory
                    clear T PT dPT_dP SS;

                elseif ismember(dataname_bc, {'ERA5_SxW'})  %/ 4D dry static stability (S) times pressure velocity (W)
                    %/ Load the 4D Temp data
                    filename_nest = strcat(data_folder, 'ERA5_S', str_level_range, '_', select_field{i}, str_domain, str_years_bc, str_noleap, '.mat');
                    fprintf('*** Loading %s... ***\n', filename_nest);
                    load(filename_nest, 'S');
                    SS = S.(select_field{i});

                    filename_nest = strcat(data_folder, 'ERA5_W', str_level_range, '_', select_field{i}, str_domain, str_years_bc, str_noleap, '.mat');
                    fprintf('*** Loading %s... ***\n', filename_nest);
                    load(filename_nest, 'S');
                    W = S.(select_field{i});

                    %/ Static stability (K/Pa) times pressure velocity (Pa/s)
                    dataset.(dataname_bc).(select_field{i}) = SS.*W;
                    dataset.(dataname_bc).unit = 'K s^{-1}';  
                    dataset.(dataname_bc).lon   = S.lon;
                    dataset.(dataname_bc).lat   = S.lat;
                    dataset.(dataname_bc).level = S.level;
                    if contains(select_field{i}, 'monthly')
                        dataset.(dataname_bc).date_yyyymm = S.date_yyyymm;   
                    else
                        dataset.(dataname_bc).date_yyyymmdd_AllYr = S.date_yyyymmdd_AllYr;   
                    end
                    
                    %/ Spare memory
                    clear S SS W

                elseif contains(dataname_bc, {'ERA5', 'ERA5Land'})
                    % %/ Debugging
                    % datafolder = data_path_parts{whichfolder(k),1};
                    % filename_prefix = data_path_parts{whichfolder(k),2};
                    % filename_part = data_path_parts{whichfolder(k),3};
                    % dataname = dataname_bc;
                    % dataname_callfile = dataname_callfile{k};
                    % varname = varname{k};
                    % dataunitconv = dataunitconv(k);
                    % datatype = datatype{k};
                    % dataunit = dataunit{k};
                    % StepsInADay = StepsInADay(k);

                    %/ Read ERA5 / ERA5-Land data
                    dataset.(dataname_bc) = read_ERA5('datafolder',   data_path_parts{whichfolder(k),1}, 'filename_prefix', data_path_parts{whichfolder(k),2},...
                                                      'filename_part',data_path_parts{whichfolder(k),3}, ...
                                                      'dataname', dataname_bc,'dataname_callfile', dataname_callfile{k}, 'varname', varname{k},...
                                                      'dataunitconv',dataunitconv(k),'datatype',datatype{k}, 'dataunit', dataunit{k},...
                                                      'stlon', stlon,'edlon', edlon, 'stlat', stlat, 'edlat', edlat,...
                                                      'stlevel', stlevel,'edlevel', edlevel, 'slct_level', slct_level,...
                                                      'case_date', case_date, 'StepsInADay', StepsInADay(k), 'timeshift', timeshift, 'select_field', mainfield_list,...
                                                      'NumWorkers', NumWorkers); 
    
                elseif ismember(dataname_bc, {'q500Omega300'})
                    raw_data_list = {'q500', 'Omega300'};
                    for ii = 1:length(raw_data_list)
                        loadfile = strcat(data_folder, raw_data_list{ii}, '_', mainfield_list{1}, '_global',str_years_bc,'.mat');
                        disp(['Loading ', loadfile, ' ...']);
                        load(loadfile, 'S');
                        fld = fieldnames(S);
                        for f = 1:length(fld)
                            dataset.(raw_data_list{ii}).(fld{f}) = S.(fld{f});
                        end
                    end
                    clear S;

                    %/ Unit conversion 
                    %    Since Omega300 is in unit of hPa/day, [1/g * q500 * Omega300] = [100 mm day-1]
                    unit_conv = 100;
    
                    %/ Since q500 is regional, but Omega300 is global
                    ind_lon = findismember_loop(dataset.(raw_data_list{2}).lon, dataset.(raw_data_list{1}).lon);
                    ind_lat = findismember_loop(dataset.(raw_data_list{2}).lat, dataset.(raw_data_list{1}).lat);
                    Omega   = dataset.(raw_data_list{2}).(mainfield_list{1})(ind_lon,ind_lat,:);
                    g       = 9.81;
    
                    dataset.(dataname_bc).(mainfield_list{1}) = unit_conv * 1/g * dataset.(raw_data_list{1}).(mainfield_list{1}).* Omega;
                    dataset.(dataname_bc).lon = dataset.(raw_data_list{1}).lon;
                    dataset.(dataname_bc).lat = dataset.(raw_data_list{1}).lat;
                    dataset.(dataname_bc).date_yyyymmdd_AllYr = dataset.(raw_data_list{1}).date_yyyymmdd_AllYr;  
    
                elseif ismember(dataname_bc, {'BLH', 'RH500', 'RH850'}) %/ they are monthly data, directly load from nc file.
                    ncfilepath = strcat(data_path_parts(whichfolder(k), 1), data_path_parts(whichfolder(k), 2),...
                                        dataname_callfile{k}, data_path_parts(whichfolder(k), 3), str_years(2:end), '.nc');
                    disp(ncfilepath)
                    ncdisp(ncfilepath{1})
                    dataset.(dataname_bc).monthly  = double(ncread(ncfilepath{1}, varname{k}))*dataunitconv(k);
                    dataset.(dataname_bc).lon      = double(ncread(ncfilepath{1}, 'longitude'));
                    dataset.(dataname_bc).lat      = double(ncread(ncfilepath{1}, 'latitude'));
    
                    dn1900 = datenum(1900,1,1);
                    time_ori = double(ncread(ncfilepath{1},'time'));
                    dataset.(dataname_bc).date = str2num(datestr(dn1900 + time_ori/24,'yyyymm'));
    
                    if ismember(dataname_bc, {'RH500'})
                        dataset.(dataname_bc).level = double(ncread(ncfilepath{1}, 'level'));
                        ind_lv = find(dataset.(dataname_bc).level == 500);
                        dataset.(dataname_bc).monthly = squeeze(dataset.(dataname_bc).monthly(:,:,ind_lv,:));
    
                    elseif ismember(dataname_bc, {'RH850'})
                        dataset.(dataname_bc).level = double(ncread(ncfilepath{1}, 'level'));
                        ind_lv = find(dataset.(dataname_bc).level == 850);
                        dataset.(dataname_bc).monthly = squeeze(dataset.(dataname_bc).monthly(:,:,ind_lv,:));
    
                    elseif ismember(dataname_bc, {'SST'})
                        [lon_2D, lat_2D] = meshgrid(dataset.(dataname_bc).lon, dataset.(dataname_bc).lat);
                        lon_2D = lon_2D'; lat_2D = lat_2D';
                        lon_2Dto1D = reshape(lon_2D, [], 1);
                        lat_2Dto1D = reshape(lat_2D, [], 1);
    
                        ind_land = which_on_land('pos_lon_array', lon_2Dto1D, 'pos_lat_array', lat_2Dto1D, 'pos_hgt_array', [], 'hgt_range', []);
    
                        SST_3Dto2D = reshape(dataset.(dataname_bc).monthly, [], size(dataset.(dataname_bc).monthly,3));
    
                        SST_3Dto2D(ind_land, :) = nan;  %/ set nan on land grids
                        dataset.(dataname_bc).monthly = reshape(SST_3Dto2D, size(dataset.('SST').monthly,1), size(dataset.('SST').monthly,2), size(dataset.('SST').monthly,3));
    
                    end
    
                elseif ismember(dataname_bc, {'T2mSST'})
                    loadfile = strcat(data_folder, 'T2m_', mainfield_list{1}, '_global',str_years_bc,'.mat');
                    disp(['Loading ', loadfile, ' ...']);
                    load(loadfile, 'S');
                    fld = fieldnames(S);
                    for f = 1:length(fld)
                        dataset.('T2m').(fld{f}) = S.(fld{f});
                    end
    
                    loadfile = strcat(data_folder, 'SST_', mainfield_list{1}, '_global',str_years_bc,'.mat');
                    disp(['Loading ', loadfile, ' ...']);
                    load(loadfile, 'S');
                    fld = fieldnames(S);
                    for f = 1:length(fld)
                        dataset.('SST').(fld{f}) = S.(fld{f});
                    end
                    clear S;
    
                    %/ obtain indices of land grids and merge SST with T2m
                    [lon_2D, lat_2D] = meshgrid(dataset.('SST').lon, dataset.('SST').lat);
                    lon_2D = lon_2D'; lat_2D = lat_2D';
                    lon_2Dto1D = reshape(lon_2D, [], 1);
                    lat_2Dto1D = reshape(lat_2D, [], 1);
    
                    ind_land = which_on_land('pos_lon_array', lon_2Dto1D, 'pos_lat_array', lat_2Dto1D, 'pos_hgt_array', [], 'hgt_range', []);
    
                    SST_3Dto2D = reshape(dataset.('SST').monthly, [], size(dataset.('SST').monthly,3));
                    T2m_3Dto2D = reshape(dataset.('T2m').monthly, [], size(dataset.('T2m').monthly,3));
    
                    SST_3Dto2D(ind_land, :) = T2m_3Dto2D(ind_land, :);  %/ replace SST values with T2m values on land grids
            %         SST_3Dto2D(ind_land, :) = -9999;  %/ testing
                    dataset.(dataname_bc).monthly = reshape(SST_3Dto2D, size(dataset.('SST').monthly,1), size(dataset.('SST').monthly,2), size(dataset.('SST').monthly,3));
                    dataset.(dataname_bc).lon     = dataset.('SST').lon;
                    dataset.(dataname_bc).lat     = dataset.('SST').lat;
                    dataset.(dataname_bc).date    = dataset.('SST').date;
    
                    %/ plot to check the assimilated results
                    contf_data = mean(dataset.(dataname_bc).monthly, 3);
                    contf_levels = 240:5:320;
                    contf_lon = dataset.('SST').lon;
                    contf_lat = dataset.('SST').lat;
                    contf_unit = '';
                    colmap = brewermap(length(contf_levels)-1, '*Spectral');
                    cbar_interval = 2;
                    pcolor_mode = 1;
                    fontsize = 18;
                    create_fig = 1;
                    grid_mode = 0;
                    cbar_mode = 1;
                    coast_col = [.4 .4 .4];
                    coast_wi = 1.5;
                    fig_fmt = 'png';
                    titlename = dataname_bc;
                    if savefig  savepath = strcat(plotting_folder, strrep(titlename, ' ', '_'));
                    else        savepath = []; end
    
                    close all;
                    %---- use m_map ----%
                    plot_contfmap('contf_data', contf_data, 'contf_lon', contf_lon, 'contf_lat', contf_lat, 'contf_levels', contf_levels,...
                                  'contf_unit', contf_unit, 'colmap', colmap, 'cbar_interval', cbar_interval, 'pcolor_mode', pcolor_mode,...
                                  'titlename', titlename, 'savepath', savepath, 'fig_fmt', fig_fmt,...
                                  'glb_data_mode', 1, 'glb_plateau_mode', 0, 'plateau_hgt', [], 'plateau_col', [],...
                                  'map_lon_lower', -179, 'map_lon_upper', 180, 'map_lat_lower', -90, 'map_lat_upper', 90, 'coast_col', coast_col, 'coast_wi', coast_wi,...
                                  'fontsize', fontsize,  'create_fig', create_fig, 'grid_mode', grid_mode, 'cbar_mode', cbar_mode)
                    fprintf('done \n')

                else
                    %/ Read daily CERA-20C data
                    dataset.(dataname_bc) = read_CERA('time_dim', time_dim, ...
                                                      'datafolder',   data_path_parts{whichfolder(k),1}, 'filename_prefix', data_path_parts{whichfolder(k),2},...
                                                      'filename_part',data_path_parts{whichfolder(k),3}, ...
                                                      'dataname', dataname_bc,'dataname_callfile', dataname_callfile{k}, 'varname', varname{k},...
                                                      'dataunitconv',dataunitconv(k),'datatype',datatype{k},...
                                                      'stlon', stlon,'edlon', edlon, 'stlat', stlat, 'edlat', edlat,...
                                                      'stlevel', stlevel,'edlevel', edlevel, 'slct_level', slct_level,...
                                                      'case_date', case_date, 'StepsInADay', StepsInADay(k), 'timeshift', timeshift, 'select_field',mainfield_list,...
                                                      'fc_steps', fc_steps{k}, 'NumWorkers', NumWorkers);
                end
    
                %/ After unitconversion, set -ve values (condensation) = 0 (do it after interpolation, otherwise -ve values are created after the interpolation). 
                if set_conden_to_zero && ...
                    ismember(dataname_bc, {'ERA5Land_Etot', 'ERA5Land_Esoil', 'ERA5Land_Etran', 'ERA5Land_Eopwt', 'ERA5Land_Ecano',  'ERA5Land_Esnow',...
                                           'CERA_P', 'CERA_E', 'ERA5_P', 'ERA5_E', 'EM_P', 'EM_E', 'CRU_P', 'GPCP_P', 'TPR_P', 'TPR_E', 'GPCC_P',...
                                           'HARv2_P', 'HARv2_E', 'GLEAM_E', 'GLASS_E', 'CMORPH_P', 'IMERG_P'}) 

                    if ~isempty(find(dataset.(dataname_bc).(mainfield_list{1}) < 0, 1))
                        warning('!!! -ve value found in %s!. Setting to zeros... !!!', dataname_bc);
                    end
                    dataset.(dataname_bc).(mainfield_list{1})(dataset.(dataname_bc).(mainfield_list{1}) < 0) = 0; 
                end
    
                %/ If omitting leap days
                if noleap
                    if isequal(mainfield_list{1}, 'daily')
                        ind_leap = find(mod(dataset.(dataname_bc).date_yyyymmdd_AllYr, 1e4) == 229);
                        if time_dim == 3
                            dataset.(dataname_bc).(mainfield_list{1})(:,:,ind_leap) = [];
                        elseif time_dim == 4
                            dataset.(dataname_bc).(mainfield_list{1})(:,:,:,ind_leap) = [];
                        else
                            error('code not ready for time_dim == %d!', time_dim);
                        end

                        %/ Check if the date-related field exists, update it if so.
                        date_related_flds = {'date_yyyymmdd_AllYr', 'casedate_UTC', 'casedate_UTC_dt'};
                        for ii = 1:length(date_related_flds)
                            if isfield(dataset.(dataname_bc), date_related_flds{ii}) 
                                dataset.(dataname_bc).(date_related_flds{ii})(ind_leap) = [];
                            end
                        end
                    else
                        error('code not ready for noleap = 1 for %s!', mainfield_list{1});
                    end
                    fprintf('*** %d leap days are omitted (noleap = 1) ***\n', length(ind_leap))
                end

                %/ Save data in mat
                if savemat
                    S = dataset.(dataname_bc);
                    disp(['Writing into ', filename, ' ...'])
                    save(filename, 'S','-v7.3');
                    clear S;
                    disp(strcat({'Time taken: '},num2str(toc),{' sec for writting.'}))
                end
                
                %/ Show min, mean, mean(abs()) and max values
                fprintf('*** %s %s: min = %.2g ***\n',      dataname_bc, select_field{:}, min(dataset.(dataname_bc).(select_field{:}), [],  'all', 'omitnan'));
                fprintf('*** %s %s: max = %.2g ***\n',      dataname_bc, select_field{:}, max(dataset.(dataname_bc).(select_field{:}), [],  'all', 'omitnan'));
                fprintf('*** %s %s: mean = %.2g ***\n',     dataname_bc, select_field{:}, mean(dataset.(dataname_bc).(select_field{:}),     'all', 'omitnan'));
                fprintf('*** %s %s: mean abs = %.2g ***\n', dataname_bc, select_field{:}, mean(abs(dataset.(dataname_bc).(select_field{:})),'all', 'omitnan'));
                
                if contains(select_field{:}, 'monthly')
                    disp(dataset.(dataname_bc).date_yyyymm([1:3,end-2:end]));
                else
                    disp(dataset.(dataname_bc).date_yyyymmdd_AllYr([1:3,end-2:end]));
                end
            end
        end
    end
    fprintf('!!! read_any: Job Completed !!!\n');

    %/ Shut down parpool
    poolobj = gcp('nocreate');
    delete(poolobj);

end