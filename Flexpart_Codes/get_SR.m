function [SR, SR_filename_suffix, dataset] = get_SR(varargin)

    pnames = {'SR', 'project_name', 'param', 'dataset', 'dataname', 'data_folder', 'valid_WSV_dataname',...
              'optimal_rr', 'expmnt', 'ldirect', 'output_res', 'dt', 'RHc_dqc_scheme', 'from_basin', 'masterfolder',....
              'slct_data', 'model', 'exp', 'ens', 'area_mean_or_sum', 'ins_or_acc',...
              'grouping_mode', 'grouping_mode_labels', 'regime_ver', 'basin_name', 'basin_catalog', 'thres_agglomerate',...
              'select_field', 'DivByAreaOf', 'mth', 'str_mth', 'st_month', 'st_day', 'ed_month', 'ed_day', 'year_list', 'shared_year_list',...
              'anom_base_period', 'derive_SR', 'recompute_SR', 'savemat', 'savemat_prefix', 'savenc', 'nc_remark', 'NumWorkers'};  

    dflts  = cell(1, length(pnames));

    [          SR,   project_name,   param,   dataset,   dataname,   data_folder,   valid_WSV_dataname,...
               optimal_rr,  expmnt,   ldirect,    output_res,    dt,  RHc_dqc_scheme,   from_basin,   masterfolder,....
               slct_data,   model,   exp,   ens,   area_mean_or_sum,   ins_or_acc,...
               grouping_mode,   grouping_mode_labels,   regime_ver,   basin_name,    basin_catalog,   thres_agglomerate,...
               select_field,   DivByAreaOf,   mth,    str_mth,   st_month,  st_day,    ed_month,  ed_day,   year_list,    shared_year_list,...
               anom_base_period,    derive_SR,   recompute_SR,   savemat,   savemat_prefix,   savenc,  nc_remark,  NumWorkers] ...
                   = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments
%% 
    %==========================================================================
    %/      Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 16 Feb 2024
    %/
    %/ Description: This function is currently designd to retrieve 
    %/              moisture source-receptor (SR) matrix.
    %/              
    %/              It works for the following variables:
    %/              1. WaterSip Variables (WSV): Pm, Pm_frac, etc.
    %/              2. Any other reanalysis data saved in a correct matname format
    %/              3. Any other CMIP data saved in a correct matname format
    %==========================================================================
    if isempty(param)
        error('Input ''param''!');
    else
        run(param);  %/ load masterfolder, data parameter
    end
    if isempty(project_name)
        error('Input ''project_name''!');
    end
    if isempty(NumWorkers)  NumWorkers = 20;                            end
    if isempty(derive_SR)   derive_SR = 1;                              end
    if isempty(mth)         error('Please determine ''mth'' [mth=1-12 (monthly), mth=13-16 (seasonally), mth=0 (annual), mth=17 (Apr-Sep), mth=18 (Oct-Mar)] !');  end
    if length(mth) > 1      error('''mth'' should only be a value, not an array!');    end
    if isempty(str_mth)     str_mth = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'MAM', 'JJA', 'SON', 'DJF', 'AMJJAS', 'ONDJFM', 'JFD'}; end
    basin_name_strrep = strrep(strrep(strrep(basin_name, '-', '_'), ' ', '_'), '.', '_');
    str_basin = strcat('_', basin_name);
    if ~isempty(grouping_mode)
        str_grouping_mode = strcat('_', strrep(grouping_mode, '-', '_'));
    else
        str_grouping_mode = '';
    end

    if mth == 0 
        str_timescale_bc = '_annual';
    else
        str_timescale_bc = strcat('_', str_mth{mth});  
    end

    if ~isempty(ens)
        if ischar(ens)  ens = {ens};  end %/ always make it a cell
        ens(cellfun(@isempty,ens)) = [];  %/ remove empty cell from the query [IMPORTANT!]
    end

    %/ [for bwd Pm ONLY] Normalize the total water budget (area_mean_or_sum = 'sum') by dividing it by the basin area
    str_DivByAreaOf = [];
    if ~isempty(DivByAreaOf) 
        if isequal(DivByAreaOf, 'basin_itself') %/ Otherwise, will just load the region as given by 'DivByAreaOf'
            DivByAreaOf = basin_name;  
        end

        if isequal(area_mean_or_sum, 'sum')
            str_DivByAreaOf = sprintf('_DivByAreaOf%s', DivByAreaOf);

        elseif isequal(area_mean_or_sum, 'mean')
            str_DivByAreaOf = '';
            DivByAreaOf = [];
            warning('''DivByAreaOf'' is for regional sum. Will ignore it.');
        end
    end

    if isequal(grouping_mode, 'basin-wise')  
        str_thres_agllo = sprintf('_agglo%.1fperc', thres_agglomerate);  
    else
        str_thres_agllo = [];
    end
    str_years_meteo = sprintf('_%d-%d', year_list(1), year_list(end));

    %/ Anomaly
    if ~isempty(anom_base_period)
        if isequal(anom_base_period, 'self')  %/ No fixed period; simply take anomaly based on the data period
            anom_base_period = year_list;
        end
        str_anom_base = sprintf('_anom%d_%d', anom_base_period(1), anom_base_period(end));
    else
        str_anom_base = [];
    end
        
    if ischar(slct_data)
        slct_data = {slct_data};  %/ Convert char into cell
    end

    %/ Set slct_data_raw
    if contains(slct_data, '_global')        
        slct_data_raw = slct_data(1:end-7); %/ e.g., get 'T2m' from 'T2m_global'

    elseif contains(slct_data, {'Pm', 'Cf_map'})
        if contains(expmnt, '_CERA_')     %/ when FLEXPART forcing data is CERA-20C
            forcing_P = 'CERA_P';
            % forcing_P = 'P';
        elseif contains(expmnt, '_EA_')   %/ when FLEXPART forcing data is ERA5
            if contains(expmnt, '0.5deg')
                forcing_P = 'ERA5_P_05';
            else
                forcing_P = 'ERA5_P';
            end
        else
            error('What is the forcing data for the tracking experiment %s?', expmnt)
        end
        if contains(slct_data, {'Pm'})
            slct_data_raw = {'Pm', forcing_P};     
        elseif contains(slct_data, {'Cf_map'})
            slct_data_raw = {'Cf_map', forcing_P}; 
        end

    elseif isequal(slct_data, 'Pm_frac')
        slct_data_raw = {'Pm'};
    elseif isequal(slct_data, 'CMF')
        slct_data_raw = {'Pm'};             %/ We will derive continental moisture feedback (CMF) from Pm
    else
        slct_data_raw = slct_data;
    end

    %/ Set slct_data_callfile, slct_data_new to handle CMIP6 data    
    if ~isempty(model) && ~isempty(exp) && ~isempty(ens)
        str_ens = strcat('_', strjoin(ens, '_'));
        slct_data_callfile = strcat(model,'_',exp,'_',slct_data,str_ens);
        slct_data_new      = strcat(slct_data,'_',strrep(model,'-','_'),'_',strrep(exp,'-','_'));           %/ Note that field name does not allow '-', replacing it with '_'
    else
        slct_data_callfile = slct_data;
        slct_data_new      = slct_data;
    end

    %/ Load all seasonal values for 'domain-wise-exact' b4 getting the annual value!
    mth_list_bc = mth;      %/ by default
    if any(ismember(slct_data_raw, {'Pm', 'Cf_map'}))
        if isequal(grouping_mode, 'domain-wise-exact')
            if mth == 0
                %==============================================================
                %/ NOTE: mth = 19 -> season = 7 -> JFD(0) - such that we can
                %/       take into account all months of 40 years in order to 
                %/       compute the annual value correctly.
                %/
                %/       BUT you can still compute D(-1)JF(0) with mth == 16.
                %==============================================================
                mth_list_bc = [13:15,  19];  %/ annual
            elseif mth == 20
                mth_list_bc = [13, 15, 19];  %/ nonJJA
            end
        end
    end

    %/ str_regime_ver
    if isequal(grouping_mode, 'domain-wise-exact')
        str_regime_ver = strcat('_regime',regime_ver); %/ just to make sure SR files are different for different regime versions.
    else
        str_regime_ver = []; 
    end

    %/ Inner monthly/seasonal loop
    for mth_bc = mth_list_bc
        %/ Dates
        if ~isempty(st_month) && ~isempty(st_day) && ~isempty(ed_month) && ~isempty(ed_day)
            str_timescale = '_interannual'; %/ To imply the mean/acc of the same specific period for different year.
            str_dates = sprintf('_%d%02d%02d-%d%02d%02d', year_list(1), st_month, st_day, year_list(end), ed_month, ed_day);
        else
            if mth_bc == 0        
                str_timescale = '_annual';    
                str_dates     = str_years_meteo; %/ [CAVEAT]: Please do NOT change this line, otherwise it will affect the reading of the processed WSV data.
            else
                str_timescale = strcat('_', str_mth{mth_bc});
                str_dates     = strcat(str_years_meteo, str_timescale); 
            end
        end

        if contains(slct_data_new, {'Pm', 'Cf_map'})
            SR_filename_suffix  = strcat(savemat_prefix, '_', slct_data_callfile, '_', area_mean_or_sum, '_', ins_or_acc, str_DivByAreaOf, '_', grouping_mode, str_regime_ver, str_thres_agllo, str_dates, str_anom_base);
        else
            SR_filename_suffix  = strcat(savemat_prefix, '_', slct_data_callfile, '_', area_mean_or_sum, '_', ins_or_acc, str_DivByAreaOf, str_thres_agllo, str_dates, str_anom_base);
        end
        SR_filename = char(strcat(data_folder, 'SR_', SR_filename_suffix, '_', basin_name, '.mat'));

        %/ Check if SR matrix exists (if recompute_SR == 0)
        flag_computing = 0;
        if recompute_SR == 0
            if isfile(SR_filename) 
                fprintf('!!! SR data is found! Loading %s !!! \n', SR_filename)
                load(SR_filename, 'S');  
            else
                if derive_SR
                    flag_computing = 1;
                else
                    warning('!!! No SR data is found at %s! And not to compute as queried. \n', SR_filename)
                    return;
                end
            end
        else
            if derive_SR
                flag_computing = 1;
            else
                warning('Set ''derive_SR = 1'' if you want to recompute SR !!');
                return;
            end
        end

        if contains(select_field,'monthly') 
            date_fld = 'date_yyyymm';
        else
            date_fld = 'date_yyyymmdd_AllYr';
        end

        %/ Load tracking param here (in case it was skipped by flag_computing == 0)
        if any(ismember(slct_data_raw, valid_WSV_dataname))
            [str_RHc_dqc, str_remark, str_src_z, str_years, WSV_lon, WSV_lat,...
             maxtraj_day, str_BLH_factor, str_optimal, str_traj_rm_jump, dt_slct, data_folder,...
             plotting_folder, WSV_dir, str_expmntinfo, basin_catalog, str_domain, str_domain_trajfile, str_sharpcut, forcing] = ...
                load_tracking_param('expmnt', expmnt, 'ldirect', ldirect, 'output_res', output_res, 'dt', dt, 'year_list', year_list,...
                                    'optimal_rr', optimal_rr, 'RHc_dqc_scheme', RHc_dqc_scheme, 'from_basin', from_basin, 'masterfolder', masterfolder);

        end

        %/ Compute SR matrix
        if flag_computing
            %/ 1. Data reading
            %============= Read reanalysis fields ====================%   
            select_data = findismember_loop(dataname, setdiff(slct_data_raw, valid_WSV_dataname)); %/ unique() to avoid bug.
            for k = select_data
                str_years_meteo_bc = str_years_meteo;

                if ~isempty(model) && ~isempty(exp) && ~isempty(ens) %/ Always reload the field to avoid messing up with different CMIP models, exp and ens
                    flag_load = 1;  
                elseif isfield(dataset, dataname{k})    %/ For ordinary field, skip it if it has been loaded.
                    flag_load = 0;
                else
                    flag_load = 1;
                end

                if flag_load == 0
                    disp(['!!! ', dataname{k}, ' has already been loaded. Skip data retrieval. !!!'])
                else
                    domain_str = '_global';
                    if ismember(dataname{k}, cmip_dataname) %/ if it is a CMIP6 variable
                        str_ens = strcat('_', strjoin(ens, '_'));
                        dataname_bc = strcat(model,'_',exp,'_',dataname{k}, str_ens);
                        loadfile = strcat(data_folder, dataname_bc, '_', select_field, domain_str, str_years_meteo_bc, '.mat');
                        disp(['Loading ', loadfile, ' ...']);
                        load(loadfile, 'S');
                        fld = fieldnames(S);
                        main_fld = strcat(dataname{k},'_',strrep(exp, '-', '_'));
                        dataset.(dataname{k}).(select_field) = S.(main_fld);
                        dataset.(dataname{k}).model = model;
                        dataset.(dataname{k}).exp   = exp;
                        dataset.(dataname{k}).ens   = ens;
                        remaining_fld = setdiff(fld, main_fld);
                        for f = 1:length(remaining_fld)
                            dataset.(dataname{k}).(remaining_fld{f}) = S.(remaining_fld{f});
                        end
                    else
                        dataname_bc = dataname{k};
                        loadfile = strcat(data_folder, dataname_bc, '_', select_field, domain_str, str_years_meteo_bc, '.mat');
                        disp(['Loading ', loadfile, ' ...']);
                        load(loadfile, 'S');
                        fld = fieldnames(S);
                        for f = 1:length(fld)
                            dataset.(dataname{k}).(fld{f}) = S.(fld{f});
                        end
                    end
                    clear S;
                end

                if contains(select_field,'monthly') && isfield(dataset.(dataname{select_data}), 'date_yyyymm')
                    date_fld = 'date_yyyymm';
                else
                    date_fld = 'date_yyyymmdd_AllYr';
                end

                %/ NOTE: Only output slct_dates is enough, 
                %        the queried data will be subset according to slct_dates
                %        in L316.
                if isempty(year_list)
                    error('year_list is empty! Check your code!');
                end
                % disp(year_list([1, end]))
                disp(dataset.(dataname{k}))
                [~, ~, ~, dataset.(dataname{k}).slct_dates] = ...
                       compute_anom('dataset', dataset, 'var', dataname{k}, 'select_field', select_field,...
                       'mth', mth_bc, 'st_month', st_month, 'st_day', st_day, 'ed_month', ed_month, 'ed_day', ed_day,...
                       'year_list', year_list, 'compute_anombyAM', 0, 'alpha', []);
            end

            %============= Read WaterSip Vars (retrieve_WSV) =========%
            WSV = [];               %/ clear up WSV to spare memory first.
            select_data_WSV = findismember_loop(valid_WSV_dataname, slct_data_raw);
            for k = select_data_WSV
                WSV_name = valid_WSV_dataname{k};
                WSV_data_filename = string(strcat(data_folder, 'WSV_', WSV_name, str_basin,...
                                          strrep(str_expmntinfo, ' ', '_'), str_dates(2:end), '.mat'));
                WSV_dir_bc     = WSV_dir;
                maxtraj_day_bc = maxtraj_day;
                ldirect_bc     = ldirect;
                from_basin_bc  = from_basin;
                str_remark_bc  = str_remark;

                %/ Load the WSV data
                if isfile(WSV_data_filename)
                    if from_basin ~= 0
                        %/ ALWAYS load data since Pm field is different for different hs!
                        fprintf('!!! The queried %s for the %s region is found. Loading... !!! \n', WSV_name, basin_name);
                        L = load(WSV_data_filename, 'WSV');
                        flds = fieldnames(L.WSV);
                        for ii = 1:length(flds)
                            WSV.(flds{ii}) = L.WSV.(flds{ii});  %/ recursively load into WSV struct. Otherwsie there will always be one WSV field.
                        end
                    else
                        %/ No need to always load data if from_basin == 0;
                        if isfield(WSV, WSV_name)                                   %/ skip it if it has been loaded.
                            disp(['!!! ', WSV_name, ' has already been loaded. Skip data retrieval. !!!'])
                        else
                            fprintf('!!! The queried %s is found. Loading... !!! \n', WSV_name);
                            L = load(WSV_data_filename, 'WSV');
                            flds = fieldnames(L.WSV);
                            for ii = 1:length(flds)
                                WSV.(flds{ii}) = L.WSV.(flds{ii});  %/ recursively load into WSV struct. Otherwsie there will always be one WSV field.
                            end
                        end
                    end
                else
                    fprintf('*** %s Not Found ***\n', WSV_name)
                    fprintf('*** Now implementing retrieve_WSV... ***\n')
                    [WSV_data_daily, WSV_dates, ~] = retrieve_WSV('project_name', project_name, 'WSV_name', WSV_name,  'WSV_dir',  WSV_dir_bc, 'year_list', year_list, 'mth', mth_bc,...
                                                                'st_month', st_month, 'st_day', st_day, 'ed_month', ed_month, 'ed_day', ed_day, 'str_RHc_dqc', str_RHc_dqc, 'ldirect', ldirect_bc, 'from_basin', from_basin_bc, 'maxtraj_day', maxtraj_day_bc, 'str_optimal', str_optimal,...
                                                                'dt_slct', dt_slct, 'str_traj_rm_jump', str_traj_rm_jump,...
                                                                'str_BLH_factor', str_BLH_factor, 'str_remark', str_remark_bc, 'str_src_z', str_src_z, 'str_domain', str_domain, 'NumWorkers', NumWorkers,...
                                                                'basin_name', basin_name, 'basin_catalog', basin_catalog, 'compute_anombyAM', 0, 'forcing', forcing);
                    %/ WaterSip Var (WSV)
                    WSV.(WSV_name).(select_field)   = WSV_data_daily;   
                    WSV.(WSV_name).slct_dates       = WSV_dates;

                    if ismember(WSV_name, {'uptake', 'BL_uptake', 'optimal_trajtime',...
                                           'rr_tot_L', 'rr_tot_NLL', 'rr_tot_NLO', 'rr_tot_Land', 'rr_dominant'})
                        WSV.(WSV_name).lon  = WSV_lon;
                        WSV.(WSV_name).lat  = WSV_lat;
                    else
                        WSV.(WSV_name).lon  = WSV_lon;
                        WSV.(WSV_name).lat  = WSV_lat(2:end-1);
                    end
                    if from_basin  WSV.(WSV_name).basin_name = basin_name;  
                    else           WSV.(WSV_name).basin_name = '';              end

                    %/ Save only the whole year or seasonal WSV.
                    if savemat
                        fprintf('*** Saving WSV data: %s *** \n', WSV_data_filename{:})
                        save(WSV_data_filename, 'WSV', '-v7.3');
                    end
                end
            end

            %/ Set 'final_data_daily'
            if isempty(select_data) && isempty(select_data_WSV)
                error('Check if ''%s'' is valid a dataname from your parameter module ''%s.m''!', slct_data_raw{:}, param);

            elseif ~isempty(select_data) && isempty(select_data_WSV)   %/ if only reanalysis data is loaded, but not WSV data
                if length(select_data) > 1   
                    error('multiple reanalysis data are selected. Check your code!');  
                end
                final_data_dates = dataset.(dataname{select_data}).slct_dates;
                
                ind_date = findismember_loop(dataset.(dataname{select_data}).(date_fld), final_data_dates);
                
                % slct_data
                if isempty(final_data_dates)
                    error('final_data_dates is empty! Check your code!')
                elseif isempty(ind_date)
                    disp(dataname{select_data})
                    disp(final_data_dates([1:3,end-2:end]))
                    disp(dataset.(dataname{select_data}).(date_fld)([1:3,end-2:end]))
                    error('ind_date is empty! Check your code!')
                else
                    % fprintf('size(ind_date) = \n'); size(ind_date)
                end
                final_data_daily = dataset.(dataname{select_data}).(select_field)(:,:,ind_date);
                final_data_lon   = dataset.(dataname{select_data}).lon;
                final_data_lat   = dataset.(dataname{select_data}).lat;
            else 
                final_data_daily = WSV.(WSV_name).(select_field);  %/ no need to subset dates, retrieve_WSV has done it.
                final_data_lon   = WSV.(WSV_name).lon;
                final_data_lat   = WSV.(WSV_name).lat;
                final_data_dates = WSV.(WSV_name).slct_dates;
            end

            %/ A quick check on the data by plotting
            % close all;
            % figure
            % contf_levels = [0, 0.005, 0.02, 0.05, 0.1, 0.25, 0.5, 1, 1.5, 2, 3, 4]; %/ uneven level
            % colmap = my_colormap(length(contf_levels)-1,'precip3_11lev');
            % colmap(1,:) = [1 1 1];
            % plot_contfmap('contf_lon', final_data_lon, 'contf_lat', final_data_lat, 'contf_data', final_data_daily(:,:,9), 'contf_levels', contf_levels, 'colmap', colmap)

            if isempty(select_data) && isempty(select_data_WSV)
                error('The enquiried data is not in the list of dataname! Add it if missing!');
            elseif isempty(final_data_daily)
                error('final_data_daily is empty! Check your code!');
            end

            %/ 2. [Optional] Upscale to 100% (Use with CAUTION)
            if ismember(slct_data_new, {'Pm_frac', 'Pm_frac_adj'})
                %/ NOTE: Here we compute ratio contributions (in %), weighted by grid area.  WARNING: will output NaN when Pm = 0 everywhere on some days!
                Area_global      = calc_grid_area('lon', final_data_lon, 'lat', final_data_lat);
                final_data_daily = (final_data_daily.*Area_global)./sum(final_data_daily.*Area_global, [1,2], 'omitnan')*100;
            end

            %/ 3. Obtain reg_2D (based on the data's lon, lat)
            if contains(slct_data, '_global')
                slct_reg = {'global'};

            elseif contains(slct_data_new, {'Pm', 'Cf_map'})  
                %/ Then 'slct_reg' will contain local and nonlocal sources!
                if isequal(grouping_mode, 'global')  %/ useful for obtaining the moisture attribution rate
                    slct_reg     = {'global'};

                elseif isequal(grouping_mode, 'basin-wise')
                    slct_reg     = {'hydrosheds_TPbasins_oceans'};

                elseif contains(grouping_mode, 'domain-wise')          %/ Check also line 474 for grouping
                    if isequal(basin_name, 'TP')
                        other_TP_basins = [];
                        slct_reg        = [basin_name, other_TP_basins, grouping_mode_labels];
                        
                    elseif isequal(basin_name, 'Inner_TP')             %/ Make sure it is defined by reg_extractor.m in regmeansum.m
                        other_TP_basins = setdiff([basin_catalog.name], {'S. Inner TP', 'N. Inner TP'});
                        slct_reg        = [basin_name, other_TP_basins, grouping_mode_labels];
                        
                    elseif isequal(basin_name, 'TP-Indus_TP-Ganges')   %/ Make sure it is defined by reg_extractor.m in regmeansum.m
                        other_TP_basins = setdiff([basin_catalog.name], {'TP-Indus', 'TP-Ganges'});
                        slct_reg        = [basin_name, other_TP_basins, grouping_mode_labels];
                    else
%                             other_TP_basins = setdiff([basin_catalog.name], basin_name);
%                             slct_reg        = [basin_name, other_TP_basins, grouping_mode_labels];  %/ TP basins / grid cells + quiried regime domains
                        slct_reg =  [basin_catalog.name, grouping_mode_labels];
                    end

                elseif isequal(grouping_mode, 'land-ocean')
                    slct_reg     = {'land', 'ocean'};

                elseif isequal(grouping_mode, 'local-nonlocal')  %/ Check also line 474 for grouping
                    if isequal(basin_name, 'TP')
                        slct_reg = {'hydrosheds_TP_oceans'};
                    else
                        slct_reg = {'hydrosheds_TPbasins_oceans'}; %<- no need to change
                    end
                else
                    slct_reg = grouping_mode;
                end

            else %/ otherwise for reanalysis variables like P or E, we compute its regional sum/mean based on the given basin.
                slct_reg = basin_name; 
            end

            %/ 4.1 Get regional mean or sum (weighted by grid area)
            %/     NOTE: Here we input 'slct_data' instead of 'slct_data_new' 
            %/           to keep the original var name to get the correct unit.
            % slct_reg
            % DivByAreaOf
            [regional_data_daily, unit, reg_name_list] = regmeansum('slct_data', slct_data, 'data', final_data_daily,...
                                                    'area_mean_or_sum', area_mean_or_sum, 'lon', final_data_lon, 'lat', final_data_lat,...
                                                    'dates', final_data_dates, 'slct_reg', slct_reg, 'DivByAreaOf', DivByAreaOf, 'mth', mth_bc,...
                                                    'regime_ver', regime_ver, 'data_folder', data_folder);
            fprintf('regional_data_daily (%s) unit: %s\n', slct_data_new{:}, unit)

%             size(final_data_daily)
%             size(regional_data_daily)
%             length(find(regional_data_daily ~= 0))
            if ismember(slct_data_new, {'TPR_E'})
                %/ There are unreasonably high values on some of the days in TPR_E,
                %/ replacing them with the value nearby (for now)

                %/ final_data_dates([5563        5575        5576])
                %/ 19950325
                %/ 19950406
                %/ 19950407
                ind = find(regional_data_daily > 1000, 1);  
                if ~isempty(ind) 
                    error('there are unreasonably huge value!');
                end
            end

            %/ 4.2 Get the common dates between P and the target data 
            if ismember(slct_data_new, {'Pm', 'Pm_frac_real', 'Cf_map', 'Cf_map_frac_real'})
                P_dates             = dataset.(forcing_P).(date_fld);
                common_dates        = intersect(final_data_dates, P_dates);
                ind_Pm_frac_date    = findismember_loop(final_data_dates, common_dates);
                ind_P_date          = findismember_loop(P_dates,          common_dates);
                regional_Pm_daily   = regional_data_daily(:,ind_Pm_frac_date);  %/ Update
                final_data_dates    = final_data_dates(ind_Pm_frac_date);       %/ Update

                [regional_Pm_ts, regional_Pm_sd, date_ts, ~] = daily2any('data_daily', regional_Pm_daily, 'dates', common_dates,...
                                                                         'mth', mth_bc, 'ins_or_acc', ins_or_acc);


                P_daily    = dataset.(forcing_P).(select_field)(:,:,ind_P_date);
                P_lon      = dataset.(forcing_P).lon; 
                P_lat      = dataset.(forcing_P).lat; 

                [regional_P_daily, regional_P_daily_unit, ~] = regmeansum('slct_data', forcing_P, 'data', P_daily,...
                                                      'area_mean_or_sum', area_mean_or_sum, 'lon', P_lon, 'lat', P_lat,...
                                                      'dates', common_dates, 'slct_reg', DivByAreaOf, 'DivByAreaOf', DivByAreaOf, 'mth', mth_bc,...
                                                      'regime_ver', regime_ver, 'data_folder', data_folder);
                fprintf('regional_P_daily unit: %s\n', regional_P_daily_unit)

                %/ Since P (CERA) has no data on Jan 1 for all years, and no Dec 31 for 2010, there will be 1-2 days fewer. 
                %/ Hence, ignore the warning from 'daily2any' about the incomplete coverage.
                [regional_P_ts,   ~,  ~,  ~] = daily2any('data_daily', regional_P_daily, 'dates', common_dates,...
                                                         'mth', mth_bc, 'ins_or_acc', ins_or_acc);

                % size(regional_Pm_daily)
                % size(regional_P_daily)
                if ismember(slct_data_new, {'Pm_frac_real', 'Cf_map_frac_real'}) 
                    regional_data_daily = regional_Pm_daily./regional_P_daily*100;
                    regional_data_ts    = regional_Pm_ts./regional_P_ts*100;          %/ compute the real fraction (hence there will be Unattributed fraction due to incomplete attribution)
                    regional_data_sd    = std(regional_data_ts, 0, 2, 'omitnan'); %/ rmb to 'omitnan'!!
                    unit = '%';
                elseif ismember(slct_data_new, {'Pm', 'Cf_map'}) 
                    regional_data_daily = regional_Pm_daily;
                    regional_data_ts    = regional_Pm_ts;
                    regional_data_sd    = regional_Pm_sd;
                end

            elseif ismember(slct_data_new, {'Pm_frac_adj'})
                P_dates             = dataset.(forcing_P).(date_fld);
                common_dates        = intersect(final_data_dates, P_dates);
                ind_Pm_frac_date    = findismember_loop(final_data_dates, common_dates);
                ind_P_date          = findismember_loop(P_dates,          common_dates);
                regional_data_daily = regional_data_daily(:,ind_Pm_frac_date); %/ Update
                final_data_dates    = final_data_dates(ind_Pm_frac_date);      %/ Update
                P_daily  = dataset.(forcing_P).(select_field)(:,:,ind_P_date);
                P_lon    = dataset.(forcing_P).lon; 
                P_lat    = dataset.(forcing_P).lat; 

                [regional_P_daily, ~, ~] = regmeansum('slct_data', forcing_P, 'data', P_daily,...
                                                      'area_mean_or_sum', 'sum', 'lon', P_lon, 'lat', P_lat,...
                                                      'dates', common_dates, 'slct_reg', DivByAreaOf,  'DivByAreaOf', DivByAreaOf, 'mth', mth_bc,...
                                                      'regime_ver', regime_ver, 'data_folder', data_folder);

                %/ 4.3 Check no-rain or no-Pm days!! [IMPORTANT to avoid BUG!]
                %/      1. Due to no rain on some days, 
                %/          the total ratio contributions = 0. 
                %/          This might cause a misleadingly large 'Elsewhere'!
                %/          Set nans for those days!
                %/      2. Since P_LA does NOT perfectly match P, for some small
                %/          areas we may find non-zero Pm but zero P, or zero Pm but non-zero P!
                %/          Workaround: remove all zero-P or zero-Pm days.
                ind_no_rain = find(regional_P_daily == 0 | sum(regional_data_daily,1) == 0); %<- Important!
                if ~isempty(ind_no_rain)
                    fprintf('*** [mth_bc: %d] [basin: %s]: Detected %d days when P(basin) = 0 OR sum(Pm) = 0! Setting values on those days to be NaN. ***\n',...
                            mth_bc, basin_name, length(ind_no_rain));
                end
                regional_data_daily(:,ind_no_rain) = nan;
                regional_P_daily(:,ind_no_rain)    = nan;  %<- Don't forget to modify P data also!

                %/ 5. Special processing of Pm_frac_adj & daily2any
                PxPm_frac = regional_data_daily.*regional_P_daily;

                [regional_PxPm_frac_ts, ~, ~, ~] = daily2any('data_daily', PxPm_frac, 'dates', common_dates,...
                                                             'mth', mth_bc, 'ins_or_acc', 'acc');

                [regional_P_ts,   ~, date_ts, ~] = daily2any('data_daily', regional_P_daily, 'dates', common_dates,...
                                                             'mth', mth_bc, 'ins_or_acc', 'acc');
                regional_data_ts = regional_PxPm_frac_ts./regional_P_ts;   %/ sum(P.*Pm_frac)./sum(P)
                regional_data_sd = std(regional_PxPm_frac_ts./regional_P_ts, 0, 2, 'omitnan'); 

            else
                %/ NOTE: daily2any also works for monthly data now
                [regional_data_ts, regional_data_sd, date_ts, ~] = daily2any('data_daily', regional_data_daily, 'dates', final_data_dates,...
                                                                             'mth', mth_bc, 'ins_or_acc', ins_or_acc);
            end
            regional_data_clim = mean(regional_data_ts, 2, 'omitnan');

            %/ 5.2. Correct the unit (if needed)
            unit_ts = unit;
            if contains(unit_ts, {'day^{-1}', '/day', 'month^{-1}', '/month'}) && isequal(ins_or_acc, 'acc')
                if mth_bc == 0
                    unit_ts = strrep(unit_ts, 'day', 'yr');
                elseif ismember(mth_bc, 1:12)
                    unit_ts = strrep(unit_ts, 'day', 'month');
                elseif mth_bc >= 13
                    unit_ts = strrep(unit_ts, 'day', 'season');
                end
            end

            %/ 6. Perform groupings (for Pm, Pm_frac, Pm_frac_adj, Pm_frac_real ONLY)
            if ismember(slct_data_new, {'Pm', 'Pm_frac_real', 'Cf_map', 'Cf_map_frac_real'}) 
                % if ~isequal(area_mean_or_sum, 'sum') && ismember(slct_data_new, {'Pm', 'Pm_frac_real'})
                %     error_message = [ 'You must set "area_mean_or_sum" to "sum" for %s to get reasonable result!\n',...
                %                       'You can divide the summmed water budget by the area of your target region by setting ''DivByAreaOf''!\n' ];
                %     error(error_message, slct_data_new);
                % end
                if contains(grouping_mode, 'domain-wise')
                    %/ Make 'local' and 'nonlocal_TP' categories
                    ind_local       = findismember_loop(reg_name_list, basin_name);
                    ind_TPbasins    = findismember_loop(reg_name_list, [basin_catalog.name]);
                    ind_nonlocal_TP = setdiff(ind_TPbasins, ind_local);  
                    ind_domain      = setdiff(1:length(reg_name_list), [ind_local, ind_TPbasins']);
                    grouped_src     = ['local', 'nonlocal TP',  reg_name_list(ind_domain)', 'Unattributed']';

                    if isequal(basin_name, 'TP')
                        %/ Insert a zero row of 'nonlocal TP' for 'TP' for the convienence of coding 
                        grouped_regional_data_daily = [regional_data_daily(ind_local,:); zeros(1, size(regional_data_daily,2));...
                                                       regional_data_daily(ind_domain,:)];
                        grouped_regional_data_ts    = [regional_data_ts(ind_local,:);    zeros(1, size(regional_data_ts,2));...
                                                       regional_data_ts(ind_domain,:)];
                        grouped_regional_data_clim  = [regional_data_clim(ind_local,:);  zeros(1, size(regional_data_clim,2));...
                                                       regional_data_clim(ind_domain,:)];
                        grouped_regional_data_sd    = [regional_data_sd(ind_local,:);    zeros(1, size(regional_data_sd,2));...
                                                       regional_data_sd(ind_domain,:)];            
                    else
                        grouped_regional_data_daily = [regional_data_daily(ind_local,:); sum(regional_data_daily(ind_nonlocal_TP,:), 1);...
                                                       regional_data_daily(ind_domain,:)];
                        grouped_regional_data_ts    = [regional_data_ts(ind_local,:);    sum(regional_data_ts(ind_nonlocal_TP,:), 1);...
                                                       regional_data_ts(ind_domain,:)];
                        grouped_regional_data_clim  = [regional_data_clim(ind_local,:);  sum(regional_data_clim(ind_nonlocal_TP,:), 1);...
                                                       regional_data_clim(ind_domain,:)];
                        grouped_regional_data_sd    = [regional_data_sd(ind_local,:);    std(sum(regional_data_ts(ind_nonlocal_TP,:), 1), 'omitnan');...
                                                       regional_data_sd(ind_domain,:)];
                    end                   

                    %/ Add 'Unattributed' category
                    if ismember(slct_data_new, {'Pm', 'Cf_map'})
                        total_daily = regional_P_daily;
                        total_ts    = regional_P_ts;
                    elseif ismember(slct_data_new, {'Pm_frac_real'})
                        total_daily = 100;
                        total_ts    = 100;
                    else
                        error('code not set!')
                    end
                    Unattributed_daily = total_daily      - sum(grouped_regional_data_daily, 1);
                    Unattributed_ts    = total_ts         - sum(grouped_regional_data_ts,    1);
                    Unattributed_clim  = mean(total_ts,2) - sum(grouped_regional_data_clim,  1);
                    Unattributed_sd    = std(Unattributed_ts, 'omitnan');

                    grouped_regional_data_daily(end+1,:) = Unattributed_daily;
                    grouped_regional_data_ts(end+1,:)    = Unattributed_ts;
                    grouped_regional_data_clim(end+1,:)  = Unattributed_clim;
                    grouped_regional_data_sd(end+1,:)    = Unattributed_sd;

                elseif isequal(grouping_mode, 'basin-wise')
                    ind_TPbasins = findismember_loop(reg_name_list, [basin_catalog.name]);
                    if isequal(basin_name, 'TP')  %/ keep all TP-basins as sig. src
                        I_sig           = unique([ind_TPbasins; find(regional_data_clim >= thres_agglomerate)], 'stable');
                        I_insig         = setdiff(find(regional_data_clim < thres_agglomerate), ind_TPbasins); 

                        %/ Two sorting, one among TP basins, another among the %remaining
                        ind_remaining = setdiff(I_sig, ind_TPbasins);
                        [~, I1]   = sort(regional_data_clim(ind_TPbasins), 'descend');
                        [~, I2]   = sort(regional_data_clim(ind_remaining), 'descend');
                        I = [I1; length(ind_TPbasins)+I2];
                        fprintf('*** Local recycling of entire TP = %.2f%% ***\n', sum(regional_data_clim(ind_TPbasins)))
                    else
                        I_sig           = find(regional_data_clim >= thres_agglomerate);
                        I_insig         = find(regional_data_clim < thres_agglomerate);
                        [~, I] = sort(regional_data_clim(I_sig), 'descend');
                    end
                    %/ Format: [(Source Names), (Contributions)]
                    grouped_regional_data_daily = [regional_data_daily(I_sig(I),:); sum(regional_data_daily(I_insig,:), 1)];
                    grouped_regional_data_ts    = [regional_data_ts(I_sig(I),:);    sum(regional_data_ts(I_insig,:),    1)];
                    grouped_regional_data_clim  = [regional_data_clim(I_sig(I),:);  sum(regional_data_clim(I_insig,:),  1)];
                    grouped_src                 = [reg_name_list(I_sig(I)); 'Elsewhere'];

                elseif ismember(grouping_mode, {'land-ocean', 'global'})
                    grouped_regional_data_daily = regional_data_daily;
                    grouped_regional_data_ts    = regional_data_ts;
                    grouped_regional_data_clim  = regional_data_clim;
                    grouped_regional_data_sd    = regional_data_sd;
                    grouped_src                 = reg_name_list;

                elseif isequal(grouping_mode, 'local-nonlocal')
                    ind_local    = findismember_loop(reg_name_list, basin_name);
                    ind_nonlocal = setdiff(1:length(reg_name_list), ind_local);

                    %/ Summing up over all nonlocal contributions.
                    grouped_regional_data_daily = [regional_data_daily(ind_local,:); sum(regional_data_daily(ind_nonlocal,:),1)];
                    grouped_regional_data_ts    = [regional_data_ts(ind_local,:); sum(regional_data_ts(ind_nonlocal,:),1)];
                    grouped_regional_data_clim  = [regional_data_clim(ind_local,:); sum(regional_data_clim(ind_nonlocal,:),1)];
                    grouped_regional_data_sd    = [regional_data_sd(ind_local,:); std(sum(regional_data_ts(ind_nonlocal,:),1),'omitnan')];
                    grouped_src                 = {'local', 'nonlocal'}';  %/ make it a coloum vector
                else
                    grouped_regional_data_daily = regional_data_daily;
                    grouped_regional_data_ts    = regional_data_ts;
                    grouped_regional_data_clim  = regional_data_clim;
                    grouped_regional_data_sd    = regional_data_sd;
                    grouped_src                 = reg_name_list;  %/ make it a coloum vector
                end

            elseif ismember(slct_data_new, {'Pm_frac', 'Pm_frac_adj'}) 
                if contains(grouping_mode, 'domain-wise')
                    %/ Make 'local' and 'nonlocal_TP' categories
                    ind_local       = findismember_loop(reg_name_list, basin_name);
                    ind_TPbasins    = findismember_loop(reg_name_list, [basin_catalog.name]);
                    ind_nonlocal_TP = setdiff(ind_TPbasins, ind_local);
                    ind_domain      = setdiff(1:length(reg_name_list), [ind_local, ind_TPbasins']);
                    grouped_src = ['local', 'nonlocal TP',  reg_name_list(ind_domain)', 'Elsewhere']';

                    if isequal(basin_name, 'TP')
                        %/ Insert a zero row of 'nonlocal TP' for 'TP' for the convienence of coding 
                        grouped_regional_data_daily = [regional_data_daily(ind_local,:); zeros(1, size(regional_data_daily,2));...
                                                       regional_data_daily(ind_domain,:)];
                        grouped_regional_data_ts    = [regional_data_ts(ind_local,:);    zeros(1, size(regional_data_ts,2));...
                                                       regional_data_ts(ind_domain,:)];
                        grouped_regional_data_clim  = [regional_data_clim(ind_local,:);  zeros(1, size(regional_data_clim,2));...
                                                       regional_data_clim(ind_domain,:)];
                        grouped_regional_data_sd    = [regional_data_sd(ind_local,:);  zeros(1, size(regional_data_sd,2));...
                                                       regional_data_sd(ind_domain,:)];            

                    else
                        grouped_regional_data_daily = [regional_data_daily(ind_local,:); sum(regional_data_daily(ind_nonlocal_TP,:), 1);...
                                                       regional_data_daily(ind_domain,:)];
                        grouped_regional_data_ts    = [regional_data_ts(ind_local,:);    sum(regional_data_ts(ind_nonlocal_TP,:), 1);...
                                                       regional_data_ts(ind_domain,:)];
                        grouped_regional_data_clim  = [regional_data_clim(ind_local,:);  sum(regional_data_clim(ind_nonlocal_TP,:), 1);...
                                                       regional_data_clim(ind_domain,:)];
                        grouped_regional_data_sd    = [regional_data_sd(ind_local,:);  std(sum(regional_data_ts(ind_nonlocal_TP,:), 1), 'omitnan');...
                                                       regional_data_sd(ind_domain,:)];
                    end                   

                    %/ Add 'Elsewhere' category
                    elsewhere_daily = 100 - sum(grouped_regional_data_daily, 1);
                    elsewhere_ts    = 100 - sum(grouped_regional_data_ts,    1);
                    elsewhere_clim  = 100 - sum(grouped_regional_data_clim,  1);
                    elsewhere_sd    = std(elsewhere_ts, 'omitnan');

                    grouped_regional_data_daily(end+1,:) = elsewhere_daily;
                    grouped_regional_data_ts(end+1,:)    = elsewhere_ts;
                    grouped_regional_data_clim(end+1,:)  = elsewhere_clim;
                    grouped_regional_data_sd(end+1,:)    = elsewhere_sd;

                elseif isequal(grouping_mode, 'basin-wise')
                    ind_TPbasins = findismember_loop(reg_name_list, [basin_catalog.name]);
                    if isequal(basin_name, 'TP')  %/ keep all TP-basins as sig. src
                        I_sig           = unique([ind_TPbasins; find(regional_data_clim >= thres_agglomerate)], 'stable');
                        I_insig         = setdiff(find(regional_data_clim < thres_agglomerate), ind_TPbasins); 

                        %/ Two sorting, one among TP basins, another among the %remaining
                        ind_remaining = setdiff(I_sig, ind_TPbasins);
                        [~, I1]   = sort(regional_data_clim(ind_TPbasins), 'descend');
                        [~, I2]   = sort(regional_data_clim(ind_remaining), 'descend');
                        I = [I1; length(ind_TPbasins)+I2];
                        fprintf('*** Local recycling of entire TP = %.2f%% ***\n', sum(regional_data_clim(ind_TPbasins)))
                    else
                        I_sig           = find(regional_data_clim >= thres_agglomerate);
                        I_insig         = find(regional_data_clim < thres_agglomerate);
                        [~, I] = sort(regional_data_clim(I_sig), 'descend');
                    end
                    %/ Format: [(Source Names), (Contributions)]
                    grouped_regional_data_daily = [regional_data_daily(I_sig(I),:); sum(regional_data_daily(I_insig,:), 1)];
                    grouped_regional_data_ts    = [regional_data_ts(I_sig(I),:);    sum(regional_data_ts(I_insig,:),    1)];
                    grouped_regional_data_clim  = [regional_data_clim(I_sig(I),:);  sum(regional_data_clim(I_insig,:),  1)];
                    grouped_src                 = [reg_name_list(I_sig(I)); 'Elsewhere'];

                elseif ismember(grouping_mode, {'land-ocean', 'global'})
                    grouped_regional_data_daily = regional_data_daily;
                    grouped_regional_data_ts    = regional_data_ts;
                    grouped_regional_data_clim  = regional_data_clim;
                    grouped_regional_data_sd    = regional_data_sd;
                    grouped_src                 = reg_name_list;

                elseif isequal(grouping_mode, 'local-nonlocal')
                    if isequal(basin_name, 'Inner_TP')
                        error('code not ready');
                    elseif isequal(basin_name, 'TP-Indus_TP-Ganges')
                        error('code not ready');
                    end
                    ind_local = findismember_loop(reg_name_list, basin_name);
                    ind_nonlocal = setdiff(1:length(reg_name_list), ind_local);

                    %/ Summing up over all nonlocal contributions.
                    grouped_regional_data_daily = [regional_data_daily(ind_local,:); sum(regional_data_daily(ind_nonlocal,:),1)];
                    grouped_regional_data_ts    = [regional_data_ts(ind_local,:); sum(regional_data_ts(ind_nonlocal,:),1)];
                    grouped_regional_data_clim  = [regional_data_clim(ind_local,:); sum(regional_data_clim(ind_nonlocal,:),1)];
                    grouped_regional_data_sd    = [regional_data_sd(ind_local,:); std(sum(regional_data_ts(ind_nonlocal,:),1),'omitnan')];
                    grouped_src                 = {'local', 'nonlocal'}';  %/ make it a coloum vector
                end
            else
                grouped_regional_data_daily = regional_data_daily;
                grouped_regional_data_ts    = regional_data_ts;
                grouped_regional_data_clim  = regional_data_clim;
                grouped_regional_data_sd    = regional_data_sd;
                grouped_src                 = reg_name_list;
            end
            % size(grouped_regional_data_ts)
            
            %/ 7. Additional calculation (compute anomaly - departure from the climatology)
            if ~isempty(anom_base_period)
                if ~isempty(model) && ~isempty(exp) && ~isempty(ens)  

                    %/ For CMIP6 models, always subtract the value from the historical clim (1971–2014) 
                    exp_his       = 'historical';
                    % year_list_his = 1971:2014;
                    year_list_his = anom_base_period;
                    
                    if isequal(exp_his, exp) && ~isequal(year_list_his, year_list)
                        warning('For consistency, the period of historical data should match the default period (1971–2014)!');
                    end

                    ens_hist       = read_CMIP_ens('project_name', project_name, 'model_list', model,  'var',  slct_data,   'exp',  exp_his, 'ori_field', select_field);
                    ens_hist       = ens_hist{:};
    
                    %/ [IMPORTANT] Set anom_base_period = [] to avoid infinite looping!
                    %/ [IMPORTANT] Set dataset_temp = [] to make sure get_SR always load from file to avoid bug
                    SR_temp = []; dataset_temp = []; dataset_temp.placeholder = [];
                    [SR_temp, ~, dataset_temp] = get_SR('SR', SR_temp, 'project_name', project_name, 'param', param, 'dataset', dataset_temp, 'dataname', dataname, 'data_folder', data_folder, 'valid_WSV_dataname', valid_WSV_dataname,...
                                                'optimal_rr', optimal_rr, 'expmnt', expmnt, 'ldirect', ldirect, 'output_res', output_res, 'dt', dt, 'RHc_dqc_scheme', RHc_dqc_scheme, 'from_basin', from_basin, 'masterfolder', masterfolder,...
                                                'slct_data', slct_data, 'model', model, 'exp', exp_his, 'ens', ens_hist, 'area_mean_or_sum', area_mean_or_sum, 'ins_or_acc', ins_or_acc, 'grouping_mode', grouping_mode, 'grouping_mode_labels', grouping_mode_labels,...
                                                'regime_ver', regime_ver, 'basin_name', basin_name, 'basin_catalog', basin_catalog, 'thres_agglomerate', thres_agglomerate,...
                                                'select_field', select_field, 'DivByAreaOf', DivByAreaOf, 'mth', mth, 'str_mth', str_mth,...
                                                'st_month', st_month, 'st_day', st_day, 'ed_month', ed_month, 'ed_day', ed_day, 'year_list', year_list_his, 'shared_year_list', shared_year_list, 'anom_base_period', [],...
                                                'recompute_SR', recompute_SR, 'savemat', savemat, 'savemat_prefix', savemat_prefix, 'savenc', 0, 'NumWorkers', NumWorkers);
                    
                    slct_data_new_hist = char(strcat(slct_data,'_',strrep(model,'-','_'),'_',strrep(exp_his,'-','_')));  %/ Always to remove the historical clim

                    
                    disp(SR_temp.(basin_name_strrep))
                    [slct_data_new_hist, str_grouping_mode, '_date', str_timescale_bc]

                    hist_period = SR_temp.(basin_name_strrep).([slct_data_new_hist, str_grouping_mode, '_date', str_timescale_bc]);
                    ind_dates = findismember_loop(hist_period, anom_base_period);
                    % size(hist_period)
                    % size(anom_base_period)
                    % size(ind_dates)
                    
                    if length(ind_dates) ~= length(anom_base_period)
                        hist_period
                        anom_base_period
                        error('Data period do not cover the ''anom_base_period''!');
                    end
                    grouped_regional_data_ts_anom = grouped_regional_data_ts - mean(SR_temp.(basin_name_strrep).([slct_data_new_hist, str_grouping_mode, str_timescale_bc])(:,ind_dates), 2, 'omitnan');
                    
                else
            
                    % if isempty(anom_base_period)
                    %     ind_dates = 1:length(date_ts);
                    % else
                    ind_dates = findismember_loop(date_ts, anom_base_period);
                    if length(ind_dates) ~= length(anom_base_period)
                        date_ts
                        anom_base_period
                        error('Data period do not cover the ''anom_base_period''!');
                    end
                    % end
                    grouped_regional_data_ts_anom = grouped_regional_data_ts - mean(grouped_regional_data_ts(:,ind_dates), 2, 'omitnan');
                end
            else
                grouped_regional_data_ts_anom = [];
            end

            %/ Compute CMF = 1/(1-rho_land)
            if ismember(slct_data_new, {'CMF'})  
                ind_land = findismember_loop(grouped_src, 'land');
                grouped_regional_data_daily = 1./(1-grouped_regional_data_daily(ind_land,:)/100);
                grouped_regional_data_ts    = 1./(1-grouped_regional_data_ts(ind_land,:)/100);
                grouped_regional_data_clim  = 1./(1-grouped_regional_data_clim(ind_land,:)/100);
            end

            %/ 8. Store into a temporary S struct (not SR, otherwise it will
            %/    store all irrelevant variablse in SR struct)
            S = [];
            S.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_daily'])                        = grouped_regional_data_daily;
            S.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_daily_unit'])                   = unit;     %/ original unit
            S.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_', date_fld])                   = final_data_dates;
            S.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, str_timescale])                   = grouped_regional_data_ts;
            S.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, str_timescale, '_clim'])          = grouped_regional_data_clim;
            S.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, str_timescale, '_sd'])            = grouped_regional_data_sd;
            S.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, str_timescale, '_unit'])          = unit_ts;  %/ unit after processing
            S.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_src_name'])                     = grouped_src;
            S.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_date', str_timescale])          = date_ts; 
            if ~isempty(anom_base_period)
                S.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, str_timescale, str_anom_base]) = grouped_regional_data_ts_anom;
                S.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_anom_base_period'])          = anom_base_period;
            end
            if ~isempty(model)
                if contains(model, 'MME')   
                    S.(basin_name_strrep).(strcat(slct_data_new{:}, '_model_list')) = dataset.(slct_data_raw{:}).('model_list');    %/ The 'model_list' field was only stored for MME data
                    S.(basin_name_strrep).(strcat(slct_data_new{:}, '_ens_list'))   = dataset.(slct_data_raw{:}).('ens_list');      %/ The 'ens_list'   field was only stored for MME data
                end
            end
            % Do not remove the data from dataset for it causes repreated data loading -> slow.
            % if ~isempty(select_data)
                % dataset = rmfield(dataset,dataname{select_data}); 
            % end

            if savemat
                fprintf('*** Saving SR time series into %s *** \n', SR_filename)
                save(SR_filename, 'S', '-v7.3');
            end
        end

        %/ Load all fields from S into SR.
        flds = fieldnames(S);    %/ loop over all fields from all basin fields
        for ii = 1:length(flds)
            flds2 = fieldnames(S.(flds{ii}));
            for jj = 1:length(flds2)
                SR.(flds{ii}).(flds2{jj}) = S.(flds{ii}).(flds2{jj});
            end
        end

        %/ [IMPORTANT] Calculate *annual* / *non-JJA* value for 'domain-wise-exact'!
        if contains(slct_data_new, {'Pm', 'Cf_map'}) && isequal(grouping_mode, 'domain-wise-exact') ...
           && (mth == 0 || mth == 20) ...              %/ original request 
           && mth_bc == mth_list_bc(end) ...           %/ hidden   request 

            %/ Double check the setting
            if ~isequal(area_mean_or_sum, 'sum')
                error_message = [ 'You must set "area_mean_or_sum" to "sum" for %s to get reasonable result!\n',...
                                      'Let say if you want to compute the "areal mean" contribution from westerly regime,\n'...
                                      'it will weigh the area of all the westerly grids globally.\n' , ...
                                      'This is why the value can end up to be very small!\n' ];
                error(error_message, slct_data_new{:});
            end
            if ~isequal(ins_or_acc, 'acc')
                error('You must set ''ins_or_acc'' to ''acc'' for %s to get correct result!', slct_data_new{:});
            end

            %/ Initialization
            target_src_list    = SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_src_name']);
            target_year_list   = SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_date_MAM']); 
            target_seasonal_ts = nan(length(mth_list_bc), length(year_list), length(target_src_list));
            P_seasonal_ts      = nan(length(mth_list_bc), length(year_list), 1);

            cnt = 0;
            for mm = mth_list_bc
                cnt = cnt + 1;
                str_timescale2 = strcat('_', str_mth{mm});

                %/ Check forcing_P if calculating fractional contributions
                if contains(slct_data_new, {'Pm_frac'}) 
%                     fprintf('*** [get_SR]: Perform recursion to get seasonal values of P (CERA)... ***\n');
                    [SR, ~] = get_SR('SR', SR, 'project_name', project_name, 'param', param, 'dataset', dataset, 'dataname', dataname, 'data_folder', data_folder, 'valid_WSV_dataname', valid_WSV_dataname,...
                                     'slct_data', forcing_P, 'area_mean_or_sum', area_mean_or_sum, 'ins_or_acc', ins_or_acc, 'grouping_mode', [], 'grouping_mode_labels', [],...
                                     'regime_ver', regime_ver, 'basin_name', basin_name, 'basin_catalog', basin_catalog, 'thres_agglomerate', thres_agglomerate,...
                                     'select_field', select_field, 'DivByAreaOf', DivByAreaOf, 'mth', mm, 'str_mth', str_mth, 'year_list', year_list, 'shared_year_list', shared_year_list, 'anom_base_period', anom_base_period,...
                                     'recompute_SR', recompute_SR, 'savemat', savemat, 'savemat_prefix', savemat_prefix);

                    %/ Store seasonal time series of P
                    P_seasonal_ts(cnt,:,:) = SR.(basin_name_strrep).([forcing_P, str_timescale2]);
                end

                %/ Store seasonal time series of Pm_frac
                ind = findismember_loop(year_list, target_year_list);
                target_seasonal_ts(cnt,ind,:) = SR.(basin_name_strrep).(strcat(slct_data_new{:}, str_grouping_mode, str_timescale2))';
            end

            %/ NOTE: - As domain-wise-exact has a seasonally varying regimes, we need to calculate the annual mean 
            %/          by weighting the seasonal value by seasonal P.
            %/       - Likewise, as Pm_frac_adj is adjusted by seasonal P, 
            %/          we also need to restore the weighting to compute annual value.
            if contains(slct_data_new, {'Pm_frac'}) || contains(slct_data_new, {'Cf_map_frac'}) 
                ts = squeeze(sum(target_seasonal_ts.*P_seasonal_ts, 1)./sum(P_seasonal_ts, 1))'; %/ squeeze, then transpose to src x years
                ts_unit = '%';
                
            elseif contains(slct_data_new, {'Pm', 'Cf_map'}) 
                ts = squeeze(sum(target_seasonal_ts, 1))'; %/ squeeze, then transpose to src x years
                if ~isempty(DivByAreaOf) && isequal(area_mean_or_sum, 'sum')
                    ts_unit = 'mm yr^{-1}'; 
                else
                    ts_unit = 'km^{3} yr^{-1}';
                end
            else
                error('code not set!');
            end
            ts_clim = mean(ts, 2, 'omitnan');
            ts_sd   = std(ts, 0, 2, 'omitnan');

            %/ Store into SR 
            SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, str_timescale_bc])                = ts;
            SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, str_timescale_bc, '_clim'])       = ts_clim;
            SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, str_timescale_bc, '_sd'])         = ts_sd;
            SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_date', str_timescale_bc])       = SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_date_MAM']);
            SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, str_timescale_bc, '_unit'])       = ts_unit;
            if ~isempty(anom_base_period)
                ts_anom = ts - ts_clim;
                SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, str_timescale_bc, str_anom_base]) = ts_anom;
            end
        end
    end
    
    if savenc
        %/ Date string (baed on the original request 'mth')
        if ~isempty(st_month) && ~isempty(st_day) && ~isempty(ed_month) && ~isempty(ed_day)
            str_timescale = '_interannual'; %/ To imply the mean/acc of the same specific period for different year.
            str_dates = sprintf('_%d%02d%02d-%d%02d%02d', year_list(1), st_month, st_day, year_list(end), ed_month, ed_day);
        else
            if mth == 0        
                str_timescale = '_annual';    
                str_dates     = str_years_meteo; %/ [CAVEAT]: Please do NOT change this line, otherwise it will affect the reading of the processed WSV data.
            else
                str_timescale = strcat('_', str_mth{mth});
                str_dates     = strcat(str_years_meteo, str_timescale); 
            end
        end

        if contains(slct_data_new, {'Pm', 'Cf_map'})
            SR_filename_suffix  = strcat(savemat_prefix, '_', slct_data_callfile, '_', area_mean_or_sum, '_', ins_or_acc, str_DivByAreaOf, '_', grouping_mode, str_regime_ver, str_thres_agllo, str_dates, str_anom_base);
        else
            SR_filename_suffix  = strcat(savemat_prefix, '_', slct_data_callfile, '_', area_mean_or_sum, '_', ins_or_acc, str_DivByAreaOf, str_thres_agllo, str_dates, str_anom_base);
        end
        ncfilename = char(strcat(data_folder, 'SR_', SR_filename_suffix, '_', basin_name, '.nc'));
        fprintf('*** Saving SR time series into %s *** \n', ncfilename)

        %/ Always save the daily SR time series (sources x dates)
        nc_data              = SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_daily']);
        data_units           = SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_daily_unit']);
        nc_date              = [];
        date_format          = [];
        othervar             = []; 
        othervar.dates       = SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_', date_fld]);

        % slct_data
        if contains(slct_data, 'Pm')
            % nc_data     = SR.(basin_name_strrep).([slct_data_new, str_grouping_mode, str_timescale]);
            % nc_data     = SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_daily']);
            data_shortname    = 'Cb';
            data_standardname = 'Cb';
            data_longname     = sprintf('Backward contributions to %s', basin_name);
            othervar.source_name = SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_src_name']);
            % nc_date           = [];  %/ Put dates in othervars, for Pm is sources x dates, not lon x lat x dates
            % date_format       = [];
            % othervar.source_name               = SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_src_name']);
            % othervar.dates                     = SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_', date_fld]);
            % othervar.(strcat(data_shortname, str_timescale_bc)) = SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, str_timescale]);
            % othervar.years = year_list;

        elseif contains(slct_data, 'Cf_map')
            % nc_data     = SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, str_timescale]);
            % nc_data     = SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_daily']);
            data_shortname    = 'Cf';
            data_standardname = 'Cf';
            data_longname     = sprintf('Forward contributions to %s', basin_name);
            othervar.source_name = SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_src_name']);
            % nc_date = [];
            % othervar.source_name               = SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_src_name']);
            % othervar.dates                     = SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, '_', date_fld]);
            % othervar.(strcat(data_shortname, str_timescale_bc)) = SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, str_timescale]);
            % othervar.years = year_list;
        
        elseif ismember(slct_data, {'P_LA'})
            data_shortname    = 'precip_estimate'; %/ Standard name requested by intercomparison community
            data_standardname = 'precip_estimate'; %/ Standard name requested by intercomparison community
            data_longname     = 'Precipitation estimated by FLEXPART-WaterSip forced by ERA5 model data';
            othervar.('precip_estimate_sum') = sum(nc_data, 'all');

        elseif ismember(slct_data, {'ERA5_P_05'})
            data_shortname    = 'precip_era5'; %/ Standard name requested by intercomparison community
            data_standardname = 'precip_era5'; %/ Standard name requested by intercomparison community
            data_longname     = 'precipitation from ERA5';
            othervar.('precip_era5_sum') = sum(nc_data, 'all');
        else
            data_shortname    = slct_data{:};
            data_standardname = slct_data{:};
            data_longname     = sprintf('%s in %s', slct_data{:}, basin_name);
        end

        if isempty(st_month) && isempty(st_day) && isempty(ed_month) && isempty(ed_day)
            othervar.(strcat(data_shortname, str_timescale_bc)) = SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, str_timescale]);
            othervar.(strcat(data_shortname, str_timescale_bc, '_units')) = SR.(basin_name_strrep).([slct_data_new{:}, str_grouping_mode, str_timescale, '_unit']);
            othervar.years = year_list;
        end

        if isempty(nc_remark)
            if ismember(slct_data, valid_WSV_dataname)
                nc_remark = strrep(str_expmntinfo, '_', ' ');  %/ Basic tracking settinng
            else
                nc_remark = [];
            end
        end

        write_nc('ncfilename', ncfilename, 'data', nc_data, 'data_shortname', data_shortname,...
                 'data_standardname', data_standardname, 'data_longname', data_longname, 'data_units', data_units,...
                 'lon', [], 'lat', [], 'plev', [], 'time', [], 'time_unit', [], ...
                 'date', nc_date, 'date_format', date_format, 'remark', nc_remark, 'othervar', othervar);
        % ncread(ncfilename, data_shortname)
        % ncread(ncfilename, 'source_name')
    end
end