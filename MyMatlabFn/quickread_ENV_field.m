function dataset = quickread_ENV_field(varargin)

    pnames       = {        'casedates',     'year_list',         'MCS_id',    'bndry_data',    'boxreg_lon',    'boxreg_lat',    'dataset', 'target_dataname',     'model',        'exp', 'project_name',...
                        'vector_or_not',  'select_field',    'anom_period',        'noleap',      'sig_mode',         'alpha',   'slct_reg',        'zm_or_mm', 'lon_range',  'lat_range',      'stlevel',    'edlevel',...
                           'local_grad',       'FHA_ver',    'data_folder', 'interp_to_4km', 'interp_to_8km', 'interp_method', 'NumWorkers',      'str_remark',...
                          'load_or_not',     'recompute', 'search_dataset',   'matfilename',         'debug'};
    
    dflts        = cell(length(pnames), 1);
    
    [                         casedates,      year_list,            MCS_id,      bndry_data,      boxreg_lon,      boxreg_lat,       dataset,  target_dataname,       model,          exp,   project_name,...
                          vector_or_not,   select_field,       anom_period,          noleap,        sig_mode,           alpha,      slct_reg,         zm_or_mm,   lon_range,    lat_range,        stlevel,      edlevel,...
                             local_grad,        FHA_ver,       data_folder,   interp_to_4km,   interp_to_8km,   interp_method,    NumWorkers,       str_remark,...
                            load_or_not,      recompute,    search_dataset,     matfilename,           debug] = ...
            internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments
%%
    %/============================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: Feb 28, 2024
    %/
    %/ NOTE:
    %/      This function handles only the processed reanalysis, satellite, or CMIP6 model data (stored in mat)
    %/      To retrieve CMORPH prcp data, use get_MCS_prcp_field() instead.
    %/
    %/      Output dataset w/ 'select_field', lon, lat, date_UTC_AllYr, date_HKT_AllYr
    %/
    %/      'select_field': Possible entries include 
    %/                      'daily', 'daily_clim', 'monthly', 'monthly_clim',
    %/                      '30min_025deg', 'pentad', 'pentad_clim', 
    %/                      with the possible FHA fields (as suffix): '_MS', '_AC', '_CISO', '_MJO', '_IAID', '_HF'
    %/
    %/     'interp_to_4km',  'interp_to_8km', 'interp_method': ONLY apply when 'MCS_prcp' is called.
    %/===========================================================================
    
    % local_grad = 0; search_dataset = []; %/ avoid bug

    % disp(target_dataname)
    if isempty(target_dataname)
        warning('No input target_dataname. Skip loading.');
        return
    elseif iscell(target_dataname)
        if length(target_dataname) > 1
            error('quickread_ENV_field was not designed to handle multiple target_dataname!');
        else
            target_dataname = target_dataname{:};
        end
    end
    
    %/ Double check 'select_field'
    if iscell(select_field)
        if length(select_field) ~= 1
            error('quickread_ENV_field does not handle multiple select_field!');
        else
            select_field = select_field{:}; %/ cell to string
        end
    end

    %/ [IMPORTANT] Set basic information 
    if contains(select_field, 'daily')
        freq          = 'daily';
        ts_conv       = 1;    %/ from days to days (no change)
        flds_date     = 'date_yyyymmdd_AllYr';
        str_flds_date = 'str_daily_dates';
        flds_basic    = {'lon', 'lat', 'date_yyyymmdd_AllYr', 'str_daily_dates', 'remark'};
        Lr = 8;  %/ requried length of the date (8 for yyyymmdd)
        output_date_format = 'yyyymmdd';

    elseif contains(select_field, 'pentad')
        freq          = 'pentad';
        ts_conv       = 1/5;  %/ from days to pentads
        flds_date     = 'date_yyyyptd';
        str_flds_date = 'str_ptd_dates';
        flds_basic    = {'lon', 'lat', 'date_yyyyptd', 'str_ptd_dates', 'remark'};
        Lr = 6;  %/ requried length of the date (6 for yyyyptd)
        output_date_format = 'yyyyptd';

    elseif contains(select_field, 'monthly')
        freq          = 'monthly';
        ts_conv       = 1/30;  %/ from days to months (WARNING: not test yet!)
        flds_date     = 'date_yyyymm';
        str_flds_date = 'str_monthly_dates';
        flds_basic    = {'lon', 'lat', 'date_yyyymm', 'str_monthly_dates', 'remark'};
        Lr = 6;  %/ requried length of the date (6 for yyyymm)
        output_date_format = 'yyyymm';
    else
        error('code not ready!')
    end

    %/ Determine casedates and year_list (if not given)
    if isempty(casedates) && isempty(year_list)
        error('Input either casedates or year_list!');

    elseif isempty(casedates) && ~isempty(year_list)
        casedates  = date_array_gen('year_list', year_list, 'output_date_format', output_date_format, 'noleap', noleap);
        
    else
        if ~isempty(casedates) && ~isempty(year_list)
            warning('Both casedates and year_list are given. Subsetting data based on casedates instead.');
        end
        %/ Retrieve years from the date format (e.g., yyyymmdd -> yyyy)
        L = numel(num2str(casedates(1)));
        if L-4 >= 0
            conv = 1./10^(L-4);
            year_list = unique(floor(casedates.*conv));
        else
            error('Invalid format of casedates (e.g., %d)!', casedates(1));
        end
    end
    str_years  = sprintf('%d-%d', year_list(1), year_list(end));

    if ~isempty(model) && ~isempty(exp)
        flag_cmip = 1;
        if iscell(model)
            model = model{:}; %/ convert cell to char
        end
        if iscell(exp)
            exp = exp{:};     %/ convert cell to char
        end
    else
        flag_cmip = 0;
    end

    if ~isempty(stlevel) && ~isempty(edlevel)
        if flag_cmip
            if stlevel > edlevel
                str_level_range = sprintf('%dto%d', edlevel, stlevel);
            else
                str_level_range = sprintf('%dto%d', stlevel, edlevel);
            end
        else
            if stlevel > edlevel
                str_level_range = sprintf('_lv%d-%d', edlevel, stlevel);
            else
                str_level_range = sprintf('_lv%d-%d', stlevel, edlevel);
            end
        end
    else
        str_level_range = '';
    end
    
    if isempty(NumWorkers)  
        NumWorkers = 20;  
    end


    if isempty(load_or_not)
        load_or_not = 1;
    end
    % if isempty(search_dataset)
    %     search_dataset = 0;
    % end


    if noleap 
        if contains(select_field, 'daily')
            str_noleap = '_noleap'; 
        elseif contains(select_field, 'pentad')
            warning('%s data do not contain the leap days. No need to set noleap = 1.', select_field)
            str_noleap = '';
        elseif contains(select_field, 'monthly')
            warning('noleap is automatically turned off for monthly data.')
            noleap = 0;
            str_noleap = '';
        else
            error('code not ready!');
        end
    else
        str_noleap = ''; 
    end

    if isempty(boxreg_lon)
        boxreg_lon = [0 360];
    end
    if isempty(boxreg_lat)
        boxreg_lat = [-90 90];
    end

    %/ NOTE: Do NOT set str_domain = '_global', because that will overwrite the original global data file
    str_domain = sprintf('_%d-%dE_%d-%dN', boxreg_lon(1),boxreg_lon(2),boxreg_lat(1),boxreg_lat(2)); %/ save processed 4D data -> time-savor!


    %=============================== [IMPORTANT!] =====================================
    %/ Update the list of 4D vars that was saved annually (e.g., 4D MSE profile)
    var_annual_list_ori = {'MSE', 'U', 'V', 'W', 'T', 'Z', 'q', 'Q1', 'Q2', 'W_dMSE_dP', 'V_dMSE_dy', 'U_dMSE_dx', 'dMSE_dt'};
    var_annual_list     = [var_annual_list_ori, strcat('ERA5_', var_annual_list_ori)]; 

    %/ The possible suffix of FHA fields.
    FHA_pat = {'_MS', '_AC', '_CISO', '_MJO', '_IAID', '_HF'};
    if isempty(FHA_ver)
        FHA_ver = 'AllAtOnce'; %/ 'AllAtOnce': Apply FHA on all pentads in 42 year_list. 
                               %/ 'append':    Apply FHA on pentads of each year, and then append the results.
    end
    str_mth = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
    %========================================================================================
    if contains(select_field, FHA_pat) && ~contains(select_field, '_bw') 
        middle_str = strcat('_FHA_', FHA_ver);
    else
        middle_str = [];
    end

    if ~isempty(matfilename)  
        matfilename_bc = strcat('_', matfilename);   
    else                      
        matfilename_bc = '';                       
    end
    
    str_zm_or_mm = [];
    if ~isempty(zm_or_mm)
        if zm_or_mm == 1    
            str_zm_or_mm = '_zm';
        elseif zm_or_mm == 2    
            str_zm_or_mm = '_mm';
        else   
            error('Only set zm_or_mm to 1 (zonal mean) or 2 (meridional mean)!'); 
        end
    end

    %/ Does the dataname refer to a vector field? check if files exist!f
    if vector_or_not == 1   
        if str2double(target_dataname) <= 137 || ismember(target_dataname, {'WS97_124', 'WS114_124', 'WS114_137'})   %/ then assume it's model level
            var_list = strcat({'U','V'}, target_dataname);   
        else
            var_list = strcat({'u','v'}, target_dataname);     
        end
    else
        var_list = {target_dataname};                      
    end
    
    if isempty(casedates)     
        warning('no input for casedates. Will assume to retrieve all dates.');                    
    end
    dataset.placeholder = [];

    %/ For MCS precip
    if ismember(target_dataname, {'MCS_prcp'})  
        matfilename_full = strcat(data_folder, matfilename, '_', target_dataname, '.mat');
        if isfile(matfilename_full) && recompute == 0
            fprintf('*** Data has been processed. Loading %s ***\n', matfilename_full)
            load(matfilename_full, 'S');  
            flds = fieldnames(S);
            for f = 1:length(flds)
                dataset.(target_dataname).(flds{f}) = S.(flds{f});
            end
        else
            if isempty(MCS_id)    
                error('MCS_id is required to retrieve MCS_prcp!');   
            end

            %/ Retrieve MCS-related prcp (auto interp from CMORPH 8-km to 4-km)
            MCS_fullinfo      = cell(length(casedates), 2);
            MCS_fullinfo(:,1) = num2cell(casedates);
            MCS_fullinfo(:,2) = num2cell(MCS_id);     %/ NOTE: There will have to be one MCS id for each date only (BUT dates can be repeated!)

            boxreg_lon        = [90, 125];            
            boxreg_lat        = [15,  50];    

            % NOTE: set convert_to_acc == 0 -> hourly prcp rate (mm/hr)
            [MCS_prcp_field, MCS_prcp_lon, MCS_prcp_lat, MCS_prcp_date, MCS_prcp_unit] = ...
                get_MCS_prcp_field('MCS_fullinfo', MCS_fullinfo, 'bndry_data', [], 'boxreg_lon', boxreg_lon, 'boxreg_lat', boxreg_lat, 'keep_norain', 1,...
                                   'str_remark', str_remark, 'which_prcp_data', 'CMORPH', 'convert_to_acc', 0, 'NumWorkers', NumWorkers,...
                                   'lon_range', [90, 135], 'lat_range', [15, 45], 'interp_to_4km', interp_to_4km, 'interp_to_8km', interp_to_8km, 'interp_method', interp_method);

            %/ Pack into dataset.(target_dataname)          
            dataset.(target_dataname)                = [];
            dataset.(target_dataname).subdaily       = MCS_prcp_field;
            dataset.(target_dataname).lon            = MCS_prcp_lon;
            dataset.(target_dataname).lat            = MCS_prcp_lat;
            dataset.(target_dataname).date_UTC_AllYr = MCS_prcp_date;
            dataset.(target_dataname).unit           = MCS_prcp_unit;

            %/ double check the date sequence
            if ~isequal(casedates, dataset.(target_dataname).date_UTC_AllYr)
                error('The dates of the final output %s does not match input casedates!!', target_dataname);
            end

            if ~isempty(matfilename)
                S = dataset.(target_dataname);
                save(matfilename_full, 'S', '-v7.3');
                fprintf('*** Saving the processed data into %s ***\n', matfilename_full);
            end
        end
    end
    
    %/ For OTHER types of data
    if ~ismember(target_dataname, {'MCS_prcp'})  
        for k = 1:length(var_list)
            var = var_list{k};

            % %/ Update var for special conditions (to save time and space)
            % flag_rename = 0;
            % if isequal(var, 'Lcpr')  %/ For Lc*P, we simply read P is fine. Then we multipy Lc with P in the plotting function (e.g., quickplot_weather)
            %     var = 'pr';
            %     flag_rename = 1;
            % elseif isequal(var, 'ERA5_LcP')
            %     var = 'ERA5_P';
            %     flag_rename = 1;
            % end

            if flag_cmip
                ens = read_CMIP_ens('project_name',project_name,'model_list',model,'var',var,'exp',exp,...
                                    'ori_field',select_field,'remove_nest_cell',1); %/ Read the available ens 
                var_cmip6  = strcat(model,'_',exp,'_',var,str_level_range,'_',strjoin(ens,'_'));
            else
                var_cmip6  = [];
            end

            if isequal(select_field, '30min_025deg') %/ Each time step of this field should be stored separately.
                if flag_cmip
                    error('code not ready for cmip6 models with %s field!', select_field);
                end
                matfilename_full = strcat(data_folder, matfilename, '_', var, '.mat');
                
                if isfile(matfilename_full) && recompute == 0
                    if load_or_not
                        fprintf('*** Data has been processed. Loading %s ***\n', matfilename_full)
                        load(matfilename_full, 'S');  
                        flds = fieldnames(S);
                        for f = 1:length(flds)
                            dataset.(var).(flds{f}) = S.(flds{f});
                        end
                    else
                        fprintf('*** Data has been processed. Skip Loading %s ***\n', matfilename_full)
                    end
                else
                    %/ check if casedates exist, otherwise we call an error.
                    if isempty(casedates)                   
                        error('Must input casedates for 30min_025deg data!!');                   
                    end
                    if numel(num2str(casedates(1))) ~= 12   
                        error('Input casedates must in yyyymmddHHMM format (12 digits)!!');      
                    end

                    %/ should NOT use unique, since the input casedates are always unsorted. unique() will auto sort. More importantly, it avoids mistakes!!
                    casedates_bc = casedates;
                    Nt           = length(casedates_bc);

                    %/ initialization for data appending
                    subdaily       = [];
                    date_UTC_AllYr = [];
                    date_HKT_AllYr = [];

                    %/ loading processed mat data (use parfor if NumWorkers is not empty.)
                    if ~isempty(NumWorkers)
                        if isempty(gcp('nocreate')) 
                            parpool('local', NumWorkers);  
                        end

                        parfor t = 1:length(casedates_bc)
                            loadfile = strcat(data_folder, var, '_', select_field, '_', num2str(floor(casedates_bc(t)/1e4)), '.mat');
                            fprintf('[%d/%d] Loading %s ...\n', t, Nt, loadfile);
                            S = par_load(loadfile, 'S');

                            date_UTC_AllYr_int = datetime2int(S.date_UTC_AllYr);
                            ind                = findismember_loop(date_UTC_AllYr_int, casedates_bc(t));
                            if isempty(ind)   
                                error('No data can be retrieved for the given casedates. Check the input or code!');       
                            end

                            subdaily_t = S.subdaily(:,:,ind);

                            lon = S.lon;
                            lat = S.lat;

                            %/ Mask data in the basin if bndry_data is given
                            if ~isempty(bndry_data)
                                [subdaily_t, ~, ~] = which_in_boundary('bndry_data_inpoly2', bndry_data, 'map_data', subdaily_t,...
                                                                       'lon', lon, 'lat', lat, 'output_1D_or_2D', '2D');
                            end

                            %/ Subset the data if boxreg is given (do it here to reduce memory burden!)
                            if ~isempty(boxreg_lon) && ~isempty(boxreg_lat)
                                ind_lon = find(boxreg_lon(1) <= lon & lon <= boxreg_lon(2));
                                ind_lat = find(boxreg_lat(1) <= lat & lat <= boxreg_lat(2));
                                subdaily_t = subdaily_t(ind_lon, ind_lat, :);
                            end

                            %/ Appending
                            subdaily_3Dto2D = reshape(subdaily_t, [], length(ind));
                            subdaily        = [subdaily,       subdaily_3Dto2D];
                            date_UTC_AllYr  = [date_UTC_AllYr; S.date_UTC_AllYr(ind)];
                            date_HKT_AllYr  = [date_HKT_AllYr; S.date_HKT_AllYr(ind)];
                        end
                    else
                        for t = 1:length(casedates_bc) %/ Identical for-loop, but without parfor
                            loadfile = strcat(data_folder, var, '_', select_field, '_', num2str(floor(casedates_bc(t)/1e4)), '.mat');
                            fprintf('[%d/%d] Loading %s ...\n', t, Nt, loadfile);
                            S = par_load(loadfile, 'S');

                            date_UTC_AllYr_int = datetime2int(S.date_UTC_AllYr);
                            ind                = findismember_loop(date_UTC_AllYr_int, casedates_bc(t));
                            if isempty(ind)   
                                error('No data can be retrieved for the given casedates. Check the input or code!');       
                            end

                            subdaily_t = S.subdaily(:,:,ind);

                            lon = S.lon;
                            lat = S.lat;

                            %/ Mask data in the basin if bndry_data is given
                            if ~isempty(bndry_data)
                                [subdaily_t, ~, ~] = which_in_boundary('bndry_data_inpoly2', bndry_data, 'map_data', subdaily_t,...
                                                                       'lon', lon, 'lat', lat, 'output_1D_or_2D', '2D');
                            end

                            %/ Subset the data if boxreg is given (do it here to reduce memory burden!)
                            if ~isempty(boxreg_lon) && ~isempty(boxreg_lat)
                                ind_lon = find(boxreg_lon(1) <= lon & lon <= boxreg_lon(2));
                                ind_lat = find(boxreg_lat(1) <= lat & lat <= boxreg_lat(2));
                                subdaily_t = subdaily_t(ind_lon, ind_lat, :);
                            end

                            %/ Appending
                            subdaily_3Dto2D = reshape(subdaily_t, [], length(ind));
                            subdaily        = [subdaily,       subdaily_3Dto2D];
                            date_UTC_AllYr  = [date_UTC_AllYr; S.date_UTC_AllYr(ind)];
                            date_HKT_AllYr  = [date_HKT_AllYr; S.date_HKT_AllYr(ind)];
                        end
                    end

                    %/ get the lon lat
                    loadfile = strcat(data_folder, var, '_', select_field, '_', num2str(floor(casedates_bc(1)/1e4)), '.mat');
                    S = par_load(loadfile, 'S');
                    lon = S.lon;
                    lat = S.lat;

                    %/ update lon, lat if boxreg_lon and boxreg_lat are given.
                    if ~isempty(boxreg_lon) && ~isempty(boxreg_lat)
                        ind_lon = find(boxreg_lon(1) <= lon  &  lon <= boxreg_lon(2));
                        ind_lat = find(boxreg_lat(1) <= lat  &  lat <= boxreg_lat(2));
                        dataset.(var).lon = lon(ind_lon);
                        dataset.(var).lat = lat(ind_lat);
                    else
                        dataset.(var).lon = lon;
                        dataset.(var).lat = lat;
                    end

                    %/ reshape to 3D subdaily and store it
                    dataset.(var).subdaily = reshape(subdaily, length(dataset.(var).lon), length(dataset.(var).lat), []);
                    dataset.(var).date_UTC_AllYr = date_UTC_AllYr;
                    dataset.(var).date_HKT_AllYr = date_HKT_AllYr;

                    %/ double check the date sequence
                    if ~isequal(casedates, datetime2int(dataset.(var).date_UTC_AllYr))
                        error('The dates of the final output %s does not match input casedates!!', var);
                    end

                    if ~isempty(matfilename)
                        S = dataset.(var);
                        save(matfilename_full, 'S', '-v7.3');
                        fprintf('*** Saving processed data into %s ***\n', matfilename_full);
                    end
                end
            else
                %/ First of all, check if the *requested* field has been processed.
                if recompute
                    flag_exist = 0;
                else
                    [dataset, flag_exist, ~] = check_n_load_field('dataset', dataset, 'data_folder', data_folder, 'var', var, 'model', model, 'exp', exp, 'project_name', project_name,...
                                            'select_field', select_field, 'noleap', noleap, 'middle_str', middle_str, 'str_years', str_years, 'boxreg_lon', boxreg_lon, 'boxreg_lat', boxreg_lat, 'matfilename', matfilename,...
                                            'stlevel', stlevel, 'edlevel', edlevel, 'search_dataset', search_dataset, 'search_datafolder', 1, 'load_or_not', load_or_not);
                end
                
                %/ If not existed, see whether to load 'pentad' / 'daily' data
                if flag_exist == 0
                    if contains(select_field, 'pentad')
                        if recompute
                            load_raw = 1;
                            flag_exist_pentad = 0;
                        else
                            %/ ALWAYS check and load 'pentad', 'pentad_clim'.
                            %/ This will avoid calling the time-consuming daily2pentad() function.
                            flds = {'pentad', 'pentad_clim'};
                            [dataset, flag_exist_pentad, ~] = check_n_load_field('dataset', dataset, 'data_folder', data_folder, 'var', var,  'model', model, 'exp', exp, 'project_name', project_name,...
                                                                              'select_field', flds,  'noleap', noleap, 'middle_str', middle_str, 'str_years', str_years, 'boxreg_lon', boxreg_lon, 'boxreg_lat', boxreg_lat, 'matfilename', matfilename,...
                                                                              'stlevel', stlevel, 'edlevel', edlevel, 'search_dataset', search_dataset, 'search_datafolder', 1, 'load_or_not', 1);
                            if flag_exist_pentad
                                load_raw = 0;
                            else
                                load_raw = 1;
                            end
                        end
                    else
                        load_raw = 1; 
                    end
                    
                    if load_raw
                        %/ Load the raw daily/monthly data 
                        if contains(select_field, 'monthly')
                            select_field_raw = 'monthly';   
                        else
                            select_field_raw = 'daily';     %/ Otherwise, Always load from daily field
                        end
                        boxreg_lon_raw   = []; %/ default for loading global daily data.
                        boxreg_lat_raw   = []; %/ default for loading global daily data.
                        %/ Check if the *daily* field has been processed. 
                        [dataset, flag_exist_raw, matfilename_full] = check_n_load_field('dataset', dataset, 'data_folder', data_folder, 'var', var,  'model', model, 'exp', exp, 'project_name', project_name,...
                                                                         'select_field', select_field_raw,  'noleap', noleap, 'str_years', str_years, 'boxreg_lon', boxreg_lon_raw, 'boxreg_lat', boxreg_lat_raw, 'matfilename', matfilename,...
                                                                         'stlevel', stlevel, 'edlevel', edlevel, 'search_dataset', search_dataset, 'search_datafolder', 1, 'load_or_not', 1);

                        if flag_exist_raw == 0
                            warning('No %s data is found %s!', select_field_raw, matfilename_full{:});
                            return;
                        end
                    end
                end
            end

            %/ Check
            if flag_exist == 0 || load_or_not == 1
                %/ Assume the last dim = time dim
                if flag_exist == 0
                    intermediate_field = select_field_raw;
                else
                    intermediate_field = select_field;
                end
    
                time_dim = length(size(dataset.(var).(intermediate_field)));
                if time_dim ~= 3 && time_dim ~= 4    
                    error('code not set yet for non-3D geodata!');    
                end
    
                %/ Subset by boxreg_lon & boxreg_lat (before daily/pentad/monthly/FHA -> faster)
                if flag_exist == 0 && ~isempty(boxreg_lon) && ~isempty(boxreg_lat)

                    if ~isvector(dataset.(var).lon) || ~isvector(dataset.(var).lat)
                        warning('Detected that either lon or lat is likely a matrix! Now interpolating it into monotonic gridding...')
                        %/ Check if lat and lon were created by [lat,lon] = meshgrid(lat_array,lon_array) 
                        %/ or [lon,lat] = meshgrid(lon_array,lat_array)\
                        [dlon1,dlon2] = gradient(dataset.(var).lon); 
                        [dlat1,dlat2] = gradient(dataset.(var).lat); 
                        if isequal(dlat1,zeros(size(dataset.(var).lat))) 
                            res_lon = abs(mode(dlon1, 'all'));  %/ Find the most common resolution (in case of non-monotonic matrix, abs to handle decreasing or increasing lon)
                            res_lat = abs(mode(dlat2, 'all'));  %/ Find the most common resolution (in case of non-monotonic matrix, abs to handle decreasing or increasing lat)
                        else
                            res_lon = abs(mode(dlon2, 'all'));  %/ Find the most common resolution (in case of non-monotonic matrix, abs to handle decreasing or increasing lon)
                            res_lat = abs(mode(dlat1, 'all'));  %/ Find the most common resolution (in case of non-monotonic matrix, abs to handle decreasing or increasing lat)
                        end
                        lon_new = 0:res_lon:360-res_lon;
                        lat_new = -90:res_lat:90;
                        
                        %/ Update
                        NumWorkers_interp = 100;
                        dataset.(var).(intermediate_field) = my_interp('lon_old', dataset.(var).lon, 'lat_old', dataset.(var).lat, 'data', dataset.(var).(intermediate_field), 'lon_new', lon_new, 'lat_new', lat_new, 'is_global', 0, 'lon_dim', 1,...
                                                                       'NumWorkers', NumWorkers_interp);
                        dataset.(var).lon = lon_new;
                        dataset.(var).lat = lat_new;
                        dataset.(var).remark = 'Interpolated to monotonic gridding';
                    end

                    ind_lon = find(boxreg_lon(1) <= dataset.(var).lon  &  dataset.(var).lon <= boxreg_lon(2));
                    ind_lat = find(boxreg_lat(1) <= dataset.(var).lat  &  dataset.(var).lat <= boxreg_lat(2));
                    if isempty(ind_lon) || isempty(ind_lat)
                        error('Either ind_lon or ind_lat is empty! Please check!');
                    end
                    dataset.(var).lon = dataset.(var).lon(ind_lon);
                    dataset.(var).lat = dataset.(var).lat(ind_lat);
                    
                    if time_dim == 3
                        dataset.(var).(intermediate_field) = dataset.(var).(intermediate_field)(ind_lon,ind_lat,:);
                    elseif time_dim == 4
                        dataset.(var).(intermediate_field) = dataset.(var).(intermediate_field)(ind_lon,ind_lat,:,:);
                    end
                end
    
                if flag_exist
                    dates      = double(dataset.(var).(flds_date));
                    year_list  = unique(floor(dates/(10^(Lr-4))),'stable');
                end
                
                %/ Append to 'flds_basic'
                if ismember(var, var_annual_list) || time_dim == 4
                    flds_basic = [flds_basic, {'level', 'units'}];
                end
                if flag_cmip
                    flds_basic = [flds_basic, {'model', 'exp', 'ens', 'units'}];
                end
    
                % %/ Resume the var 
                % if flag_rename
                %     if isequal(var, 'pr')  %/ For Lc*P, we simply read P is fine. Then we multipy Lc with P in the plotting function (e.g., quickplot_weather)
                %         dataset = RenameField(dataset, var, 'Lcpr');  %/ Update the field name
                %         var = 'Lcpr';  %/ Update
                % 
                %     elseif contains('ERA5_LcP')
                %         dataset = RenameField(dataset, var, 'ERA5_P');  %/ Update the field name
                %         var = 'ERA5_P';  %/ Update
                %     end
                % end

                %/ Post-processing daily data
                if flag_exist == 0 && isequal(freq, 'daily')
    
                    %/ 365 days (ignore leap day) -> consistent with pentad field processing.
                    date_mmdd_OneYr = date_array_gen('year_list', 1979, 'st_month', 1, 'st_day', 1, 'ed_month', 12, 'ed_day', 31, 'output_date_format', output_date_format);
                    date_mmdd_OneYr = mod(date_mmdd_OneYr, 1e4);  
                    tot_timestamp   = length(date_mmdd_OneYr);  
    
                    %/ Create daily string arrays.
                    dataset.(var).str_daily_dates = cell(tot_timestamp, 2);
                    dataset.(var).str_daily_dates(:, 1) = cellstr(string(date_mmdd_OneYr));
                    dataset.(var).str_daily_dates(:, 2) = {str_mth{floor(date_mmdd_OneYr/1e2)}};              %/ record the month of the middle day in a given pentad. (fairest)
    
                    dates      = double(dataset.(var).date_yyyymmdd_AllYr);  %/ make sure it's double, not int64! Otherwise floor() won't work properly.
                    year_list  = unique(floor(dates/1e4),'stable');
    
                    sz = size(dataset.(var).(freq));
                    dataset.(var).daily_clim = nan([sz(1:end-1), tot_timestamp]);  %/ nlon, nlat, 365
                    
                    %/ daily clim
                    if time_dim == 3
                        for d = 1:tot_timestamp
                            ind = findismember_loop(mod(dates, 1e4), date_mmdd_OneYr(d));
                            dataset.(var).daily_clim(:,:,d) = mean(dataset.(var).daily(:,:,ind), time_dim, 'omitnan');
                        end
                    elseif time_dim == 4
                        for d = 1:tot_timestamp
                            ind = findismember_loop(mod(dates, 1e4), date_mmdd_OneYr(d));
                            dataset.(var).daily_clim(:,:,:,d) = mean(dataset.(var).daily(:,:,:,ind), time_dim, 'omitnan');
                        end
                    end
                end
                
                %/ Post-processing pentad data
                if flag_exist == 0 && isequal(freq, 'pentad')
                    %/ Create pentad string arrays.
                    date_mmdd_OneYr = date_array_gen('year_list', 1979, 'st_month', 1, 'st_day', 1, 'ed_month', 12, 'ed_day', 31, 'output_date_format', 'yyyymmdd');
                    date_mmdd_OneYr = mod(date_mmdd_OneYr, 1e4);  %/ 365 days
    
                    tot_timestamp = length(date_mmdd_OneYr)/5;                                 %/ skipped the leap day, as LinHo and Wang 2002 did.
                    dataset.(var).str_ptd_dates = cell(tot_timestamp, 2);
                    for t = 1:tot_timestamp
                        ptd_dates = date_mmdd_OneYr((1:5)+5*(t-1));
                        dataset.(var).str_ptd_dates{t, 1} = sprintf('P%d_%d-%d', t, ptd_dates(1), ptd_dates(end));
                        dataset.(var).str_ptd_dates{t, 2} = str_mth{floor(ptd_dates(3)/1e2)};              %/ record the month of the middle day in a given pentad. (fairest)
                    end
                    
                    if all(flag_exist_pentad)
                        fprintf('!!! ''pentad'' and ''pentad_clim'' have been processed already. Skip processing them. !!!\n');
                    else
                        [pentad_clim, pentad_allyr, ~, ~, date_yyyyptd] = ...
                                    daily2pentad('daily_data', dataset.(var).daily, 'dates', dataset.(var).date_yyyymmdd_AllYr);
    
                        dataset.(var).date_yyyyptd          = double(date_yyyyptd); 
                        dataset.(var).pentad                = pentad_allyr;
                        dataset.(var).pentad_clim           = pentad_clim; 
                        if flag_cmip
                            dataset.(var).model             = model;
                            dataset.(var).exp               = exp;
                            dataset.(var).ens               = ens;
                        end
                        %/ Spare memory by clearing daily field.
                        dataset.(var) = rmfield(dataset.(var),'daily');
    
                        %/ Store data
                        flds = {'pentad', 'pentad_clim'};
                        for i = 1:length(flds)
                            S = [];
                            fld_bc = flds{i};  %/ cell to str.
                            S.(fld_bc) = dataset.(var).(fld_bc); 
    
                            %/ Save together with the basic info.
                            for j = 1:length(flds_basic)
                                if isfield(dataset.(var), flds_basic{j})
                                    S.(flds_basic{j}) = dataset.(var).(flds_basic{j});
                                else
                                    warning('Missing ''%s'' in dataset.%s. Omitted.', flds_basic{j}, var);
                                end
                            end
                            
                            if flag_cmip
                                filename_bc = strcat(data_folder, var_cmip6, '_', fld_bc, str_domain, '_', str_years, str_noleap, matfilename_bc, '.mat');
                            else
                                filename_bc = strcat(data_folder, var,       '_', fld_bc, str_domain, '_', str_years, str_noleap, matfilename_bc, '.mat');
                            end
    
                            fprintf('*** Saving the processed data into %s ***\n', filename_bc);
                            save(filename_bc, 'S', '-v7.3'); 
                        end
                    end
                end
                
                %/ Post-processing monthly data (e.g., monthly_3MA_anom, monthly_3MA_anom_detrend)
                if flag_exist == 0 && isequal(freq, 'monthly')
                    if contains(select_field, 'monthly_3MA')  
                        n_MA = 3;
                    else                                        
                        n_MA = [];
                    end

                    if isequal(select_field_raw, 'daily')
                        flds = {'monthly', select_field}; %/ Then we'll save monthly and monthly_3MA separately

                        [dataset.(var).(select_field), dataset.(var).date_yyyymm_AllYr] = ...
                            daily2monthly('daily_data', dataset.(var).daily, 'dates', dataset.(var).date_yyyymmdd_AllYr, 'n_MA', n_MA);
                        
                    elseif isequal(select_field_raw, 'monthly')
                        flds = {select_field};
                        if ~isempty(n_MA)
                            dataset.(var).(select_field) = movmean(dataset.(var).(select_field_raw), n_MA, time_dim); 
                        end

                        if contains(select_field, '_anom')  
                            dataset_temp = [];
                            select_field_temp = 'monthly';
                            
                            if isempty(anom_period)
                                error('Please set ''anom_period'' in order to compute ''%s'' data!', select_field);
                            end

                            dataset_temp = quickread_ENV_field('casedates', [], 'year_list', anom_period, 'dataset', dataset_temp, 'target_dataname', target_dataname,...
                                  'model', model, 'exp', exp, 'project_name', project_name, 'vector_or_not', vector_or_not, 'select_field', select_field_temp, 'anom_period', [],...
                                  'noleap', noleap, 'sig_mode', sig_mode, 'alpha', alpha, 'zm_or_mm', zm_or_mm, 'lon_range', lon_range, 'lat_range', lat_range, 'FHA_ver', FHA_ver,...
                                  'str_remark', str_remark, 'data_folder', data_folder, 'load_or_not', 1, 'recompute', 0, 'boxreg_lon', boxreg_lon, 'boxreg_lat', boxreg_lat,...
                                  'NumWorkers', [], 'debug', 0);

                            time_dim_temp = length(size(dataset_temp.(var).(select_field_temp)));

                            %/ Monthly anomaly = monthly data - monthly clim
                            for mth = 1:12
                                ind_mthly_clim = find(mod(dataset_temp.(var).(flds_date),1e2) == mth);
                                if isempty(ind_mthly_clim)
                                    error('ind_mthly_clim is empty!');
                                end
                                % length(ind_mthly_clim)
                                if time_dim_temp == 3
                                    mthly_clim = mean(dataset_temp.(var).(select_field_temp)(:,:,ind_mthly_clim), time_dim_temp, 'omitnan');
                                elseif time_dim_temp == 4
                                    mthly_clim = mean(dataset_temp.(var).(select_field_temp)(:,:,:,ind_mthly_clim), time_dim_temp, 'omitnan');
                                end

                                % dataset.(var)
                                % dataset.(var).(flds_date)(1:5)
                                % mod(dataset.(var).(flds_date),1e2)
                                ind_all  = find(mod(dataset.(var).(flds_date),1e2) == mth);
                                if isempty(ind_all)
                                    error('ind_all is empty!');
                                end
                                if time_dim == 3
                                    if ~isempty(n_MA)
                                        dataset.(var).(select_field)(:,:,ind_all) = dataset.(var).(select_field)(:,:,ind_all) - mthly_clim;
                                    else
                                        dataset.(var).(select_field)(:,:,ind_all) = dataset.(var).(select_field_temp)(:,:,ind_all) - mthly_clim;
                                    end
                                elseif time_dim == 4
                                    if ~isempty(n_MA)
                                        dataset.(var).(select_field)(:,:,:,ind_all) = dataset.(var).(select_field)(:,:,:,ind_all) - mthly_clim;  
                                    else
                                        dataset.(var).(select_field)(:,:,:,ind_all) = dataset.(var).(select_field_temp)(:,:,:,ind_all) - mthly_clim;  
                                    end
                                end
                            end
                            if isempty(find(~isnan(dataset.(var).(select_field)), 1))
                                error('data is empty!');
                            end

                            %/ Detrend or not (e.g., monthly_3MA_anom_detrend)
                            if contains(select_field, '_detrend')  
                                %/ NOTE: detrend() only work across the columns, so reshape is needed.
                                sz = size(dataset.(var).(select_field));
                                A_2D = reshape(dataset.(var).(select_field), [], sz(end)); %/ Assuming time dim is at the last dim
                                A_2D = detrend(A_2D, 'omitnan');
                                dataset.(var).(select_field) = reshape(A_2D, sz);
                            end
                        end
                    else
                        error('code not set for select_field_raw = %s!', select_field_raw);
                    end

                    %/ Store data
                    for i = 1:length(flds)
                        S = [];
                        fld_bc = flds{i};  %/ cell to str.
                        S.(fld_bc) = dataset.(var).(fld_bc); 
    
                        %/ Save together with basic info
                        for j = 1:length(flds_basic)
                            if isfield(dataset.(var), flds_basic{j})
                                S.(flds_basic{j}) = dataset.(var).(flds_basic{j});
                            else
                                warning('Missing ''%s'' in dataset.%s. Omitted.', flds_basic{j}, var);
                            end
                        end
                        
                        if flag_cmip
                            filename_bc = strcat(data_folder, var_cmip6, str_level_range, '_', fld_bc, str_domain, '_', str_years, str_noleap, matfilename_bc, '.mat');
                        else
                            filename_bc = strcat(data_folder, var,       str_level_range, '_', fld_bc, str_domain, '_', str_years, str_noleap, matfilename_bc, '.mat');
                        end
    
                        fprintf('*** Saving the processed data into %s ***\n', filename_bc);
                        save(filename_bc, 'S', '-v7.3'); 
                    end
                end
    
                %/ Fourier Harmonics analysis (FHA)
                if flag_exist == 0
                    if contains(select_field, '_bw')  %/ using butterworth filter
                        fprintf('*** Running butterworth filtering... ***\n')

                        if contains(select_field, 'daily_MJO_')
                            ca  = 70;         %/ band start (longer period)
                            cb  = 20;         %/ band end
                        else
                            error('code not set for %s!', select_field);
                        end
                        fca     = 1/ca;       %/ the lowest frequency 
                        fcb     = 1/cb;       %/ the highest frequency 
                        fc      = [fca, fcb]; %/ cutoff frequency 
                        fs      = 1;          %/ sampling frequency; Default is 1
                        n_order = 6;          %/ filter order (4 <= m <= 6 should be adequate for most applications); Default is 6
                        
                        %/ Permute the time dimension to the first for filtfilt()
                        if time_dim == 3
                            dimorder     = [3 1 2];
                            dimorder_rev = [2 3 1];
                        elseif time_dim == 4
                            dimorder     = [4 1 2 3];
                            dimorder_rev = [2 3 4 1];
                        else
                            error('code not set!');
                        end

                        x = permute(dataset.(var).(freq), dimorder);

                        %/ Set all infinite values (e.g., NaN, Inf) to zeros
                        %/ to avoid bug from filtfilt
                        ind = length(isfinite(x));
                        if ~isempty(ind)
                            warning('%d infinite values (e.g., NaNs or Infs) found in the data. Replacing them with zeros before filtering...', length(ind))
                            x(~isfinite(x)) = 0;
                        end

                        [b,a] = butter(n_order, fc/(fs/2));           % IIR filter design
                        dataset.(var).(select_field) = permute(filtfilt(b,a,x), dimorder_rev);   % zero-phase filtering; The function operates along the first array dimension of x unless x is a row vector.
                        
                        %/ A quick check for debugging
                        if debug
                            figure
                            ind_lon = 1;
                            ind_lat = 1;
                            ind_time = 1:500;
                            y  = squeeze(dataset.(target_dataname).('daily')(ind_lon, ind_lat, ind_time));
                            y2 = squeeze(dataset.(target_dataname).(select_field)(ind_lon, ind_lat, ind_time));
                            plot(y,'k-.'); grid on ; hold on
                            plot(y2,'LineWidth',1.5);
                            legend('daily', select_field);
                        end

                        %/ Store only the filtered fields into mat field 
                        flds = fields(dataset.(var));
                        ind  = find(ismember(flds, select_field)); %/ select the fields that contain any elements in FHA_pat.
                        for i = 1:length(ind)
                            S = [];
                            fld_bc = flds{ind(i)};  %/ cell to str.
                            S.(fld_bc) = dataset.(var).(fld_bc); 
        
                            %/ Save together with the basic info.
                            for j = 1:length(flds_basic)
                                if isfield(dataset.(var), flds_basic{j})
                                    S.(flds_basic{j}) = dataset.(var).(flds_basic{j});
                                else
                                    warning('Missing ''%s'' in dataset.%s. Omitted.', flds_basic{j}, var);
                                end
                            end
        
                            if flag_cmip
                                filename_bc = strcat(data_folder, var_cmip6, str_level_range, '_', fld_bc, str_domain, '_', str_years, str_noleap, matfilename_bc, '.mat');
                            else
                                filename_bc = strcat(data_folder, var,       str_level_range, '_', fld_bc, str_domain, '_', str_years, str_noleap, matfilename_bc, '.mat');
                            end
                            fprintf('*** Saving the processed data into %s ***\n', filename_bc);
                            save(filename_bc, 'S', '-v7.3'); 
                        end

                    elseif contains(select_field, FHA_pat)
                        if ~isequal(freq, 'daily') && ~isequal(freq, 'pentad')
                            error('Fourier Harmonics analysis is not ready for select_field = %s!', select_field);
                        end
        
                        %/ Identify which FHA component is quried. Process that only.
                        for i = 1:length(FHA_pat)
                            if contains(select_field, FHA_pat{i})
                                fha = FHA_pat{i}(2:end);  %/ e.g., '_MS' -> 'MS'
                            end
                        end
        
                        if ~isempty(find(isnan(dataset.(var).(freq)), 1))         
                            warning('%s %s contains nan!!', var, freq);         
                        end
        
                        nyr = length(year_list);
                        sz = size(dataset.(var).(freq));
                        if isequal(fha, '_MS')
                            dataset.(var).(strcat(freq,'_',fha)) = nan(sz(1:end-1));
                        else
                            dataset.(var).(strcat(freq,'_',fha)) = nan(sz);
                        end
                        time_dim = length(sz);
        
                        if isequal(FHA_ver, 'AllAtOnce')
                            %/=================================================
                            %/ For a time series with a period T, 
                            %/     n/T = the frequency of the n-th harmonics, 
                            %/     T/n = the time scale.
                            %/=================================================
        
                            %/ The queried time series
                            T = sz(end);                         %/ Since n_list changes with the length of days/pentads
                            d = 0.00001;                         %/ Dummy value to avoid overlap of time scale.
                            ts = [];
                            ts.IAID  = [T    365+d]*ts_conv;     %/ Time scale of IAID (365-all days)
                            ts.AC    = [365   90+d]*ts_conv;     %/ Time scale of AC   (90-365 days)
                            ts.CISO  = [90    20+d]*ts_conv;     %/ Time scale of CISO (20-90 days) 
                            ts.MJO   = [70    20+d]*ts_conv;     %/ Time scale of MJO  (20-90 days) 
                            ts.HF    = [20     1+d]*ts_conv;     %/ Time scale of HF   (< 20 days) 
        
                            %/ The correpsonding nth harmonic
                            n_LIST = [];
                            n_LIST.MS    = 0;
                            n_LIST.IAID  = ceil(T./ts.IAID(1)):floor(T./ts.IAID(2));
                            n_LIST.AC    = ceil(T./ts.AC(1))  :floor(T./ts.AC(2));
                            n_LIST.CISO  = ceil(T./ts.CISO(1)):floor(T./ts.CISO(2));
                            n_LIST.MJO   = ceil(T./ts.MJO(1)) :floor(T./ts.MJO(2));
                            n_LIST.HF    = ceil(T./ts.HF(1))  :floor(T./ts.HF(2));
                            % size(dataset.(var).(freq))
        
                            %/ New way to compute AC, CISO signals from the entire time series 
                            tic;
                            if isequal(fha, '_MS')
                                dataset.(var).(strcat(freq,'_',fha))  = my_fourier('f', dataset.(var).(freq), 'n_list', n_LIST.(fha));    
                                dataset.(var).(strcat(freq,'_',fha))  = mean(dataset.(var).(strcat(freq,'_',fha)), time_dim); %/ 3D to 2D. Same const for all time.
        
                            elseif isequal(fha, '_HF')  %/ [WARNING] It may blow up the memory!
                                dataset.(var).(strcat(freq,'_IAID'))  = my_fourier('f', dataset.(var).(freq), 'n_list', n_LIST.IAID);  
                                dataset.(var).(strcat(freq,'_AC'))    = my_fourier('f', dataset.(var).(freq), 'n_list', n_LIST.AC);    
                                dataset.(var).(strcat(freq,'_CISO'))  = my_fourier('f', dataset.(var).(freq), 'n_list', n_LIST.CISO);  
                                dataset.(var).(strcat(freq,'_MJO'))   = my_fourier('f', dataset.(var).(freq), 'n_list', n_LIST.MJO);   
                                dataset.(var).(strcat(freq,'_HF'))    = dataset.(var).(freq) - dataset.(var).(strcat(freq,'_IAID')) ...
                                                                        - dataset.(var).(strcat(freq,'_AC')) - dataset.(var).(strcat(freq,'_CISO')) - dataset.(var).(strcat(freq,'_MS'));
                            else
                                dataset.(var).(strcat(freq,'_',fha))  = my_fourier('f', dataset.(var).(freq), 'n_list', n_LIST.(fha));    
                            end
                            fprintf('Finished my_fourier (%s) on %s (%d-%d %s) in %.2f s\n', FHA_ver, fha, round(min(ts.(fha))), round(max(ts.(fha))), freq, toc);
                            dataset.(var) = rmfield(dataset.(var), freq);  %/ Spare memory
                            
                            %/ Derive the climatological signals (For now, using all available dates)
                            dataset.(var).(strcat(freq,'_clim_',fha))        = nan([sz(1:end-1), tot_timestamp]); %/ Initialization
                            dataset.(var).(strcat(freq,'_clim_',fha,'_sig')) = nan([sz(1:end-1), tot_timestamp]); %/ Initialization
                            for tt = 1:tot_timestamp
                                if isequal(freq, 'daily')
                                    ind_tt = mod(dataset.(var).date_yyyymmdd_AllYr, 1e4) == date_mmdd_OneYr(tt);
                                    
                                elseif isequal(freq, 'pentad')
                                    ind_tt = mod(dataset.(var).date_yyyyptd, 1e4) == tt;
                                end
        
                                if time_dim == 3
                                    dataset.(var).(strcat(freq,'_clim_',fha))(:,:,tt)        = mean(dataset.(var).(strcat(freq,'_',fha))(:,:,ind_tt), time_dim);
                                    dataset.(var).(strcat(freq,'_clim_',fha,'_sig'))(:,:,tt) = ttest_sig_fn(dataset.(var).(strcat(freq,'_',fha))(:,:,ind_tt), alpha, time_dim, 0);
                                elseif time_dim == 4
                                    dataset.(var).(strcat(freq,'_clim_',fha))(:,:,:,tt)        = mean(dataset.(var).(strcat(freq,'_',fha))(:,:,:,ind_tt), time_dim);
                                    dataset.(var).(strcat(freq,'_clim_',fha,'_sig'))(:,:,:,tt) = ttest_sig_fn(dataset.(var).(strcat(freq,'_',fha))(:,:,:,ind_tt), alpha, time_dim, 0);
                                end
                            end
        
                            % reshp_sz = [sz(1:end-1), 365*ts_conv, nyr];  %/ The desired size after reshape().
                            % dataset.(var).(strcat(freq,'_clim_',fha))    = squeeze(mean(reshape(dataset.(var).(strcat(freq,'_',fha)),  reshp_sz), length(reshp_sz)));   %/ lon x lat x annual dates
                            
                            %/ Below are one-line code with reshape() & ttest_sig_fn().
                            % dataset.(var).(strcat(freq,'_clim_',fha,'_sig'))    = ttest_sig_fn(reshape(dataset.(var).(strcat(freq,'_',fha)),   reshp_sz), alpha, length(reshp_sz), 0);
        
                            if ~isempty(find(isnan(dataset.(var).(strcat(freq,'_',fha))), 1))    
                                warning('%s %s_%s contains nan!!', var, freq, fha);    
                            end
        
                            %/ Double check
                            d = dataset.(var).(strcat(freq,'_clim_',fha)) - dataset.(var).(strcat(freq,'_clim_',fha,'_sig'));
                            d(isnan(d)) = 0;
                            if max(abs(d),[],'all') > 1e-9
                                error('Non-nan values in %s_clim_%s do not match that of %s_clim_%s_sig!! Check others also!', freq, fha, freq, fha);
                            end
        
                            % dataset.(var).(strcat(freq,'_MS_AC'))  = dataset.(var).(strcat(freq,'_AC')) + dataset.(var).(strcat(freq,'_MS'));
        
                            % %/ Double check
                            % d = dataset.(var).(strcat(freq,'_clim',fha)) - dataset.(var).(strcat(freq,'_clim_',fha,'_sig'));
                            % d(isnan(d)) = 0;
                            % if max(abs(d),[],'all') > 1e-9
                            %     error('Non-nan values in %s_clim_%s do not match that of %s_clim_%s_sig!! Check others also!', freq, fha, freq, fha);
                            % end
        
                            % CISO = abs(dataset.(var).(strcat(freq,'_clim_CISO_sig'))); 
                            % AC   = abs(dataset.(var).(strcat(freq,'_clim_AC_sig')));   
                            % r    = CISO./AC;
                            % r(isnan(CISO) & ~isnan(AC)) = 0;
                            % r(~isnan(CISO) & isnan(AC)) = 1e9; %/ positive inf
                            % dataset.(var).strcat(freq,'_clim_CISO_AC_abs_ratio') = r;
        
                        elseif isequal(FHA_ver, 'Append')
                            %/ Old way to compute AC, CISO signals from the clim time series
                            tic;
                            dataset.(var).strcat(freq, '_MS') = nan([sz(1:end-1), nyr]);
        
                            %/ The queried time series
                            T = 365*ts_conv;
                            d = 0.00001;                        %/ Dummy value to avoid overlap of time scale.
                            ts.AC    = [365  90+d]*ts_conv;     %/ Time scale of AC   (90-365 days)
                            ts.CISO  = [90   20+d]*ts_conv;     %/ Time scale of CISO (20-90 days) 
                            ts.MJO   = [70   20+d]*ts_conv;     %/ Time scale of MJO  (20-90 days) 
        
                            %/ The correpsonding nth harmonic
                            n_LIST_MS    = 0;
                            n_LIST_AC    = ceil(T./ts.AC(1))  :floor(T./ts.AC(2));
                            n_LIST_CISO  = ceil(T./ts.CISO(1)):floor(T./ts.CISO(2));
                            n_LIST_MJO   = ceil(T./ts.MJO(1)) :floor(T./ts.MJO(2));
        
                            for t = 1:nyr
                                if isequal(freq, 'daily')
                                    ind = find(floor(dates/1e4) == year_list(t));
                                elseif isequal(freq, 'pentad')
                                    ind = find(floor(dates/1e2) == year_list(t));
                                end
                                if time_dim == 3
                                    MS                                             = my_fourier('f', dataset.(var).(freq)(:,:,ind), 'n_list', n_LIST_MS);  
                                    dataset.(var).(strcat(freq,'_MS'))(:,:,t)      = MS(:,:,1);
                                    dataset.(var).(strcat(freq,'_AC'))(:,:,ind)    = my_fourier('f', dataset.(var).(freq)(:,:,ind), 'n_list', n_LIST_AC);    %/ 90-365 days
                                    dataset.(var).(strcat(freq,'_CISO'))(:,:,ind)  = my_fourier('f', dataset.(var).(freq)(:,:,ind), 'n_list', n_LIST_CISO);  %/ 20-73 days 
                                    dataset.(var).(strcat(freq,'_MJO'))(:,:,ind)   = my_fourier('f', dataset.(var).(freq)(:,:,ind), 'n_list', n_LIST_MJO);   %/ 20-61 days 
                                    dataset.(var).(strcat(freq,'_HF'))(:,:,ind)    = dataset.(var).(freq)(:,:,ind) - dataset.(var).(strcat(freq,'_AC'))(:,:,ind) ...
                                                                                     - dataset.(var).(strcat(freq,'_CISO'))(:,:,ind) - dataset.(var).(strcat(freq,'_MS'))(:,:,t); %/ < 20 days
        
                                elseif time_dim == 4
                                    MS                                             = my_fourier('f', dataset.(var).(freq)(:,:,:,ind), 'n_list', n_LIST_MS);  
                                    dataset.(var).(strcat(freq,'_MS'))(:,:,:,t)    = MS(:,:,:,1);
                                    dataset.(var).(strcat(freq,'_AC'))(:,:,:,ind)  = my_fourier('f', dataset.(var).(freq)(:,:,:,ind), 'n_list', n_LIST_AC); 
                                    dataset.(var).(strcat(freq,'_CISO'))(:,:,:,ind)= my_fourier('f', dataset.(var).(freq)(:,:,:,ind), 'n_list', n_LIST_CISO);  
                                    dataset.(var).(strcat(freq,'_MJO'))(:,:,:,ind) = my_fourier('f', dataset.(var).(freq)(:,:,:,ind), 'n_list', n_LIST_MJO);  %/ 20-61 days 
                                    dataset.(var).(strcat(freq,'_HF'))(:,:,:,ind)  = dataset.(var).(freq)(:,:,:,ind) - dataset.(var).(strcat(freq,'_AC'))(:,:,:,ind) ...
                                                                                     - dataset.(var).(strcat(freq,'_CISO'))(:,:,:,ind) - dataset.(var).(strcat(freq,'_MS'))(:,:,t);
                                end
                                dataset.(var) = rmfield(dataset.(var), freq);  %/ Spare memory
        
                                if ~isempty(find(isnan(dataset.(var).(strcat(freq,'_AC'))), 1))      
                                    warning('%s %s_AC contains nan!!', var, freq);   
                                end
                                if ~isempty(find(isnan(dataset.(var).(strcat(freq,'_CISO'))), 1))    
                                    warning('%s %s_CISO contains nan!!', var, freq);    
                                end
                                if ~isempty(find(isnan(dataset.(var).(strcat(freq,'_MJO'))), 1))    
                                    warning('%s %s_MJO contains nan!!', var, freq);
                                end
                                %/ Apply FHA on climatological data.
                                clim_MS  = my_fourier('f', dataset.(var).(strcat(freq,'_clim')), 'n_list', n_LIST_MS);
                                dataset.(var).(strcat(freq,'_clim_MS'))    = clim_MS(:,:,1);  %/ The mean state is time independent.
                                dataset.(var).(strcat(freq,'_clim_AC'))    = my_fourier('f', dataset.(var).(strcat(freq,'_clim')), 'n_list', n_LIST_AC);  
                                dataset.(var).(strcat(freq,'_clim_CISO'))  = my_fourier('f', dataset.(var).(strcat(freq,'_clim')), 'n_list', n_LIST_CISO);  
                                dataset.(var).(strcat(freq,'_clim_MJO'))   = my_fourier('f', dataset.(var).(strcat(freq,'_clim')), 'n_list', n_LIST_MJO);  
                                dataset.(var).(strcat(freq,'_clim_MS_AC')) = dataset.(var).(strcat(freq,'_clim_MS')) + dataset.(var).(strcat(freq,'_clim_AC'));     %/ Constant term (i.e. time averging)
                                dataset.(var).(strcat(freq,'_clim_HF'))    = dataset.(var).(strcat(freq,'_clim')) - dataset.(var).(strcat(freq,'_clim_MS')) ...
                                                                             - dataset.(var).(strcat(freq,'_AC')) - dataset.(var).(strcat(freq,'_CISO'));      %/ HF
                            end
                            fprintf('Finished my_fourier (%s) in %.2f s\n', FHA_ver, toc);
                        else
                            error('code not set!')
                        end
        
                        %/ Store only the FHA fields into mat field 
                        flds = fields(dataset.(var));
                        ind = find(contains(flds, fha)); %/ select the fields that contain any elements in FHA_pat.
                        for i = 1:length(ind)
                            S = [];
                            fld_bc = flds{ind(i)};  %/ cell to str.
                            S.(fld_bc) = dataset.(var).(fld_bc); 
        
                            %/ Save together with the basic info.
                            for j = 1:length(flds_basic)
                                if isfield(dataset.(var), flds_basic{j})
                                    S.(flds_basic{j}) = dataset.(var).(flds_basic{j});
                                else
                                    warning('Missing ''%s'' in dataset.%s. Omitted.', flds_basic{j}, var);
                                end
                            end
        
                            if flag_cmip
                                FHA_filename_bc = strcat(data_folder, var_cmip6, str_level_range, '_', fld_bc, '_FHA_', FHA_ver, str_domain, '_', str_years, str_noleap, matfilename_bc, '.mat');
                            else
                                FHA_filename_bc = strcat(data_folder, var, str_level_range, '_', fld_bc, '_FHA_', FHA_ver, str_domain, '_', str_years, str_noleap, matfilename_bc, '.mat');
                            end
                            fprintf('*** Saving processed fourier harmonics (%s) into %s ***\n', FHA_ver, FHA_filename_bc);
                            save(FHA_filename_bc, 'S', '-v7.3'); 
                        end
                    end
                end
    
                %/ Identify the raw data (e.g., 'pentad_MJO' for 'pentad_clim_MJO_sig')
                select_field_raw = strrep(strrep(select_field, '_sig', ''), '_clim', '');
    
                %/ Subsetting by casedates (do it after daily/pentad/monthly/FHA!)
                if ~isempty(casedates)
                    L = numel(num2str(casedates(1)));
                    if L ~= Lr
                        error('casedates for daily data must be in a %d-digit format!', Lr);
                    end
    
                    %/ 1. Update climatological dates
                    queried_dates  = unique(mod(casedates, 10^(Lr-4)),'stable');    %/ Convert yyyymmdd to mmdd, unique without auto sorting ('stable')
                    existing_dates = unique(mod(dataset.(var).(flds_date), 10^(Lr-4)), 'stable');
                    ind_clim_date  = findismember_loop(existing_dates, queried_dates);
                    dataset.(var).(strcat(flds_date, '_clim')) = existing_dates(ind_clim_date);                    %/ Additional field of the updated clim date 
                    
                    %/ NOTE: My other functions (e.g., read_CMIP) may not 
                    %/       output the field specified by str_flds_date 
                    if isfield(dataset.(var), str_flds_date)
                        dataset.(var).(str_flds_date)              = dataset.(var).(str_flds_date)(ind_clim_date,:);   %/ Update the string of dates, this is a matrix; (for plotting labels)
                    end

                    %/ 2. Update calendar dates
                    queried_dates  = casedates;
                    existing_dates = dataset.(var).(flds_date);
                    ind_date       = findismember_loop(existing_dates, queried_dates);
                    if isempty(ind_date)
                        error('empty ind_date! Check ''casedates''!');
                    end

                    %/ Update the target field data 
                    if contains(select_field, '_clim')
                        ind = ind_clim_date;
                    else
                        ind = ind_date;
                    end
                    if time_dim == 3
                        dataset.(var).(select_field) = dataset.(var).(select_field)(:,:,ind);
                    elseif time_dim == 4
                        dataset.(var).(select_field) = dataset.(var).(select_field)(:,:,:,ind);
                    end
                    dataset.(var).(flds_date) = dataset.(var).(flds_date)(ind_date);
                end
    
                %/ Get the length of clim dates (after dates being subset)
                nt = length(dataset.(var).(strcat(flds_date, '_clim')));
    
                %/ If slct_reg is given, subset the region (Do it before zm_or_mm!)
                if ~isempty(slct_reg)
                    fprintf('NOTE: Values outside of the %s region will be set to NaN.', slct_reg);
                    reg_2D = reg_extractor('lon', dataset.(var).lon, 'lat', dataset.(var).lat, 'slct_reg', slct_reg, 'data_folder', data_folder, 'saveload_cond_landocean', 0, 'savemat', 0);
                    % reg_2D
                    reg_2D(~isnan(reg_2D)) = 1;
                    reg_2D(isnan(reg_2D))  = 0;
                    
                    if isempty(find(reg_2D == 1, 1))
                        error('reg_2D is all NaNs!');
                    end

                    reg_2Dto1D = reshape(reg_2D, [], 1); %/ (lon, lat) logical -> (grid) logical
                    sz = size(dataset.(var).(select_field));
                    dataset.(var).(select_field) = reshape(dataset.(var).(select_field), [sz(1)*sz(2), sz(3:end)]); %/ (lon,lat,time) -> (grid, time)
                    
                    %/ Set values NaN if outside of it 
                    if time_dim == 3
                        dataset.(var).(select_field)(~reg_2Dto1D,:) = nan;
                    elseif time_dim == 4
                        dataset.(var).(select_field)(~reg_2Dto1D,:,:) = nan;
                    end
                    dataset.(var).(select_field) = reshape(dataset.(var).(select_field), sz); %/ reshape back to original dims
                end

                %/ Take zonal/meridional area-weighted mean
                %/ NOTE: The computational procedure is different for sig field
                %/       Why? Because zm of sig field ~= zm of raw field.
                if ~isempty(zm_or_mm)
                    if contains(select_field, FHA_pat) && contains(select_field, '_sig')
                        %/ load the raw data (by recursion)
                        fprintf('*** loading select_field_raw = %s directly from ''data_folder'' (avoid bug)... ***\n', select_field_raw)
                        dataset = quickread_ENV_field('casedates', casedates, 'dataset', dataset, 'target_dataname', target_dataname,  'vector_or_not', vector_or_not,...
                                                      'model', model, 'exp', exp, 'project_name', project_name, 'select_field', select_field_raw,...
                                                      'boxreg_lon', boxreg_lon, 'boxreg_lat', boxreg_lat, 'noleap', noleap, 'sig_mode', [], 'alpha', alpha, 'FHA_ver', FHA_ver, 'str_remark', str_remark,...
                                                      'data_folder', data_folder, 'search_dataset', 0, 'matfilename', matfilename, 'recompute', 0); 
                        
                        geo_data  = dataset.(var).(select_field_raw);  
                        dataset.(var) = rmfield(dataset.(var), select_field_raw);  %/ remove field to save memory and avoid bugs.
                        % dataset.(var) = rmfield(dataset.(var), select_field);   %/ remove field to save memory and avoid bugs.
    
                        %/ First, do zm / mm
                        [geo_data_zm_mm, dataset.(var).hori_array] = ...
                            compute_zm_mm('geo_data', geo_data, 'lon',  dataset.(var).lon, 'lat', dataset.(var).lat,...
                                          'lon_range', lon_range, 'lat_range', lat_range, 'zm_or_mm', zm_or_mm);
    
                        %/ Second, do ttest on the zm / mm fields
                        sz           = size(geo_data_zm_mm);
                        nyr          = length(year_list);
                        reshp_sz     = [sz(1:end-1), nt, nyr];  %/ The desired size after reshape().
                        dataset.(var).(strcat(select_field, str_zm_or_mm)) = ttest_sig_fn(reshape(geo_data_zm_mm, reshp_sz), alpha, length(reshp_sz), 0);
    
                        %/ Double check
                        geo_data_zm_mm_clim = mean(reshape(geo_data_zm_mm, reshp_sz), length(reshp_sz));
                        d = geo_data_zm_mm_clim - dataset.(var).(strcat(select_field, str_zm_or_mm));
                        d(isnan(d)) = 0;
                        if max(abs(d),[],'all') > 1e-9
                            error('Inconsitent results between two ways of calculations of geo_data_zm_mm_clim!');
                        end
                    else
                        geo_data = dataset.(var).(select_field);
                        dataset.(var) = rmfield(dataset.(var), select_field);  %/ remove field to save memory and avoid bugs.
    
                        [dataset.(var).(strcat(select_field, str_zm_or_mm)), dataset.(var).hori_array] = ...
                            compute_zm_mm('geo_data', geo_data, 'lon',  dataset.(var).lon, 'lat', dataset.(var).lat,...
                                          'lon_range', lon_range, 'lat_range', lat_range, 'zm_or_mm', zm_or_mm);
                    end
                end
    
                %/ Compute local gradient (after zm / mm)
                if local_grad
                    str_local_grad = {'zonal', 'meridional', 'vertical'};
                    fprintf('*** Computing local %s gradient... ***\n', str_local_grad{local_grad});
                    
                    VAR = dataset.(var).(strcat(select_field, str_zm_or_mm));
                    r   = 6371e3;
                    if zm_or_mm == 1
                        lat = dataset.(var).hori_array;
                        hy  = diff(lat(1:2))*pi/180;   %/ rmb to convert degree into radian!!!! NOTE: hy can be -ve. <- correct.
                    elseif zm_or_mm == 2
                        lon = dataset.(var).hori_array;
                        hx  = diff(lon(1:2))*pi/180;   %/ rmb to convert degree into radian!!!! NOTE: hy can be -ve. <- correct.   
                    else
                        lon = dataset.(var).lon;
                        lat = dataset.(var).lat;
                        hx  = diff(lon(1:2))*pi/180;   %/ rmb to convert degree into radian!!!!
                        hy  = diff(lat(1:2))*pi/180;   %/ rmb to convert degree into radian!!!! NOTE: hy can be -ve. <- correct.
                    end
    
                    if local_grad == 1        %/ zonal 
                        if length(lon) ~= size(VAR, 1)   
                            error('Make sure the lon dim is the 1st dim of the variable!');      
                        end
                        if zm_or_mm == 1                 
                            error('As zm_or_mm == 1, impossible to take zonal local gradient.'); 
                        end
                        %/ Since after mm, the lon dim is still the 1st dim,
                        %/ but we use the midpoint of lat range to compute unit_conv.
                        [~, dVAR_dx, ~, ~] = gradient(VAR, hx);    
                        if zm_or_mm == 2 
                            unit_conv      = (1./(r*cosd(mean(lat_range))))';       %/ a const
                        else
                            unit_conv      = (1./(r*cosd(lat)))';       %/ 1 x lat
                        end
                        dVAR_dx            = unit_conv.*dVAR_dx;        %/ J kg-1 m-1.   div_A = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D))
                        str_local_grad     = '_xgrad';
                        dataset.(var).(strcat(select_field, str_zm_or_mm, str_local_grad)) = dVAR_dx;
    
                    elseif local_grad == 2    %/ meridional 
                        if zm_or_mm == 2                   
                            error('Impossible to take zonal local gradient when zm_or_mm == 2.'); 
                        end
                        if zm_or_mm == 1
                            if length(lat) ~= size(VAR, 1) 
                                error('Make sure the lat dim is the 1st dim of the variable when zm_or_mm == 1!'); 
                            end
                        else
                            if length(lat) ~= size(VAR, 2) 
                                error('Make sure the lat dim is the 1st dim of the variable!');       
                            end
                        end
                        %/ Since after zm, the lat dim becomes the 1st dim.
                        %/ Note also that [ygrad, xgrad, ~, ~] = gradient()
                        if zm_or_mm == 1
                            [~, dVAR_dy, ~]    = gradient(VAR, hy);         
                            unit_conv          = 1./r;                     
                            dVAR_dy            = unit_conv.*dVAR_dy;        %/ J kg-1 m-1.   div_A = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D))
                        else
                            [dVAR_dy, ~, ~, ~] = gradient(VAR, hy);         
                            unit_conv          = 1./r;                     
                            dVAR_dy            = unit_conv.*dVAR_dy;        %/ J kg-1 m-1.   div_A = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D))
                        end
                        str_local_grad     = '_ygrad';
                        dataset.(var).(strcat(select_field, str_zm_or_mm, str_local_grad)) = dVAR_dy;
                        
                    elseif local_grad == 3    %/ vertical 
                        P                  = dataset.(var).level;
                        if length(P) ~= size(VAR, 3)     
                            error('Make sure the p-level dim is the 3rd dim of the variable!'); 
                        end
                        
                        [~, ~, dVAR_dA, ~] = gradient(VAR);                          %/ Chain rule to do derivate with *non-uniform* spacings.
                        dP_dA              = gradient(P * 100);                      %/ get non-uniform gradient of pressure level first.
                        dP_dA              = reshape(dP_dA, [1 1 length(dP_dA)]);    %/ 1 x 1 x plev
                        dVAR_dP            = dVAR_dA./dP_dA;                         %/ J kg-1 Pa-1   elementwise division.  
                        str_local_grad     = '_zgrad';
                        dataset.(var).(strcat(select_field, str_zm_or_mm, str_local_grad)) = dVAR_dP;
                        
                    else
                        error('Incorrect input of compute_gradient_xyz. Set ''1'' for zonal gradient, ''2'' for meridional gradient, ''3'' for vertical gradient');
                    end
                end
                
                %/ Compute AS t-test (v2) -> IT IS CORRECT.
                if isequal(sig_mode, '_ASsigv2')  
                    if contains(select_field, {'_clim_AC'}) 
                        calendar_Harmon_field = strcat(freq,'_AC');
    
                    elseif contains(select_field, {'_clim_CISO'}) 
                        calendar_Harmon_field = strcat(freq,'_CISO');
    
                    elseif contains(select_field, {'_clim_MJO'}) 
                        calendar_Harmon_field = strcat(freq,'_MJO');
                    else
                        error('code not yet set!');
                    end
                    
                    %/ take zm or mm if queried
                    if ~isempty(zm_or_mm)
                        [X, ~] = compute_zm_mm('geo_data', dataset.(var).(calendar_Harmon_field), 'lon',  dataset.(var).lon, 'lat', dataset.(var).lat,...
                                               'lon_range', lon_range, 'lat_range', lat_range, 'zm_or_mm', zm_or_mm);
                    else
                        X = dataset.(var).(calendar_Harmon_field);
                    end
    
                    dataset.(var).(strcat(select_field, str_zm_or_mm, sig_mode)) = AS_ttest_wrapper('X', X, 'dates', dates, 'alpha', alpha);
                end        
            end
        end
    end
    fclose('all'); %/ avoid the bug of 'opening too many files.'
end