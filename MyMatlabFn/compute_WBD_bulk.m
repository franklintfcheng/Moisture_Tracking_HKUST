function WBD = compute_WBD_bulk(varargin)

    pnames = {'WBD_mode', 'data_source', 'param', 'his_period', 'fut_period', 'dataset', 'mth', 'str_mth',...
              'stlon', 'edlon', 'stlat', 'edlat', 'stlevel', 'edlevel', 'level_unit',...
              'select_field',  'ins_or_acc', 'sp_mode', 'lnP_coor', 'interp_mode', 'nboot', 'nboot_seed', 'data_folder', 'savemat', ...
              'recompute_yrly_data', 'recompute_WBD', 'NumWorkers', 'NumWorkers_HPC'}; 
    dflts  = cell(length(pnames), 1);
    [          WBD_mode,   data_source,   param,  his_period,   fut_period,   dataset,   mth,   str_mth,...
               stlon,   edlon,   stlat,   edlat,   stlevel,   edlevel,   level_unit,...
               select_field,    ins_or_acc,   sp_mode,   lnP_coor,    interp_mode,  nboot,   nboot_seed,   data_folder,   savemat,...
                 ~,                   recompute_WBD,   NumWorkers,   NumWorkers_HPC] = internal.stats.parseArgs(pnames, dflts, varargin{:});
%%
    %/=====================================================================
    %/       Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Updated: 21 Feb 2024
    %/
    %/  Description: The function ouputs the bulk water budgets only,
    %/               e.g., dwdt, IVTdiv terms
    %/      
    %/  Data needed: 3D Data: q, U, V, W from stlevel to edlevel.
    %/               2D Data: sp, P, E
    %/=====================================================================
    
    if mth == 0
        str_mth_bc = '';
    else
        str_mth_bc = strcat('_',str_mth{mth});
    end
    
    if iscell(select_field)  select_field = char(select_field);  end
    
    if isempty(dataset)         dataset.placeholder = [];                         end
    if ~any(ismember(WBD_mode, {'bulk', 'bulk3D'}))
        error('Invalid input of ''WBD_mode''!');
    end
    if lnP_coor                 str_lnP_coor  = '_lnP';                      else  str_lnP_coor  = '';      end 
    if ~isempty(interp_mode)    str_interp_mode = strcat('_', interp_mode);  else  str_interp_mode = '';    end
    if isempty(NumWorkers) 
        NumWorkers = 10;     %/ For i/o
    end
    if isempty(NumWorkers_HPC)
        NumWorkers_HPC = 40; %/ For calculation (vertinte)
    end

    run(param);
    if isempty(level_unit)  level_unit = 'hPa';         end

    if isequal(WBD_mode, 'bulk3D')
        var_list        = strcat(data_source, '_', {'q', 'sp', 'W', 'U', 'V', 'P', 'E'});
    else
        var_list        = strcat(data_source, '_', {'q', 'sp', 'U', 'V', 'P', 'E'});
    end
    str_lonlat          = sprintf('%d-%dE_%d-%dN', stlon, edlon, stlat, edlat);
    str_plev            = sprintf('%d_%d', stlevel, edlevel);
    his_yrs             = his_period(1):his_period(2);
    fut_yrs             = fut_period(1):fut_period(2);
    str_fut_minus_his   = sprintf('%d-%d_minus_%d-%d',fut_period(1),fut_period(2),his_period(1),his_period(2));
    
    %/ Bootstrapping (if quiried)
    if ~isempty(nboot)       
        if nboot ~= round(nboot)   error('nboot must an integer! (Set it as [] if no intension to use bootstrapping)');     end
        if nboot <= 1              error('nboot must be > 1! (Set it as [] if no intension to use bootstrapping)');         end
        if isempty(nboot_seed)     nboot_seed = 'default';    end
        str_nboot      = sprintf('_nboot%d_%s', nboot, nboot_seed);
    else
        nboot          = 1;  %/ set 1 as a dummy dimension
        str_nboot      = [];
    end


    %/ Parameter Setting for Water Budget Decomposition (WBD)
    r         = 6371e3; 
    rho_w     = 997;                           %/ Water density (kg m^-3)
    WBD       = [];
    WBD.data_source = data_source;
    WBD.his_yrs = his_yrs;
    WBD.fut_yrs = fut_yrs;

    %/ WBD filename
    WBD_matname         = sprintf('WBD_%s_%s_%s_lv%s_%s_%s%s%s_%s%s%s.mat', WBD_mode, data_source, ins_or_acc, str_plev, str_lonlat,...
                                                                           str_fut_minus_his, str_mth_bc, str_lnP_coor, sp_mode, str_interp_mode, str_nboot);    
    WBD_matname         = strcat(data_folder, strrep(WBD_matname, ' ', '_'));
    disp(WBD_matname)
    if recompute_WBD == 0 && isfile(WBD_matname)
        fprintf('!!! mat file is found: %s. Loading from it... !!! \n', WBD_matname)
        load(WBD_matname, 'WBD');
    else
        %============== Load 4D/3D data for historical period ================%
        %/ Period Loop 
        q = []; dwdt = []; MFC = []; MFCv = [];
        
        for t = 1:2  %/ Period 1 and 2
            if t == 1
                yrs        = his_yrs;
                str_hisfut = '_his';
            elseif t == 2
                yrs        = fut_yrs;
                str_hisfut = '_fut';
            end
            str_yrs = sprintf('%d-%d',yrs(1),yrs(end));

            if mth == 0
                season = [];      st_month = [];  ed_month = []; 
            elseif ismember(mth, 1:12)
                season = [];      st_month = mth; ed_month = mth; 
            elseif mth > 12
                season = mth-12;  st_month = [];  ed_month = [];  
            end
            case_date   = date_array_gen('year_list', yrs, 'season', season, 'st_month', st_month, 'ed_month', ed_month, 'dt_slct_hr', 24, 'output_date_format', 'yyyymmdd');
            select_data = findismember_loop(dataname, var_list);   %/ we have executed run('param_flexpart') in this function
%                 dataname(select_data)
            if ~isequal(length(select_data), length(var_list))  
                error('Missing variables from dataname!'); 
            end

            for k = select_data
                %/ For example, dataname_callfile for q maybe q10-1000, while its dataname is q.
                disp(strcat(dataname{k}, str_hisfut))
                if ismember(dataname{k}, strcat(data_source,'_',{'sp', 'P', 'E'}))
                    time_dim        = 3;
                    mainfield_list  = [select_field, {'lon'}, {'lat'}, {'date_yyyymmdd_AllYr'}];
                    loadfile        = strcat(data_folder, dataname{k}, '_', str_lonlat, '_', select_field, '_', str_yrs, '.mat');
                else
                    time_dim        = 4;
                    mainfield_list  = [select_field, {'lon'}, {'lat'}, {'level'}, {'date_yyyymmdd_AllYr'}];
                    loadfile        = strcat(data_folder, dataname{k}, str_plev, '_', str_lonlat, '_', select_field, '_', str_yrs, '.mat');
                end
                disp(loadfile)

                %/==================== Data Reading =======================
                %/ 1. Check if the field has been loaded, 
                %/ 2. Check if it has been processed, 
                %/ 3. Otherwise will do data processing.
                %/=========================================================
                if isfield(dataset, strcat(dataname{k}, str_hisfut))
                    disp('The quried data has been loaded into ''dataset''.');
                else
                    if isfile(loadfile) %/ find the field exists/has been loaded.
                        disp(['Processed data is found. Loading ', loadfile, ' ...']);
                        load(loadfile, 'S');
                        fld = fieldnames(S);
                        for f = 1:length(fld)
                            dataset.(dataname{k}).(fld{f}) = S.(fld{f});
                        end
                        clear S;
                    else
                        if isequal(data_source, 'CERA')   %/ Use read_CERA
                            dataset.(dataname{k}) = read_CERA('time_dim', time_dim,...
                                                              'datafolder',   data_path_parts{whichfolder(k),1}, 'filename_prefix', data_path_parts{whichfolder(k),2},...
                                                              'filename_part',data_path_parts{whichfolder(k),3}, ...
                                                              'dataname', dataname{k},'dataname_callfile', dataname_callfile{k}, 'varname', varname{k},...
                                                              'dataunitconv',dataunitconv(k),'datatype',datatype{k}, 'allyears', yrs,...
                                                              'stlon', stlon, 'edlon', edlon, 'stlat', stlat, 'edlat', edlat,...
                                                              'stlevel', stlevel,'edlevel', edlevel, 'slct_level', [],...
                                                              'case_date', case_date, 'StepsInADay', StepsInADay(k), 'timeshift', [], 'select_field',mainfield_list,...
                                                              'fc_steps', fc_steps{k}, 'NumWorkers', NumWorkers);

                        elseif isequal(data_source, 'ERA5')
                            dataset.(dataname{k}) = read_ERA5('time_dim', time_dim,...
                                                              'datafolder',   data_path_parts{whichfolder(k),1}, 'filename_prefix', data_path_parts{whichfolder(k),2},...
                                                              'filename_part',data_path_parts{whichfolder(k),3}, ...
                                                              'dataname', dataname{k},'dataname_callfile', dataname_callfile{k}, 'varname', varname{k},...
                                                              'dataunitconv',dataunitconv(k),'datatype',datatype{k}, 'allyears', yrs,...
                                                              'stlon', stlon, 'edlon', edlon, 'stlat', stlat, 'edlat', edlat,...
                                                              'stlevel', stlevel,'edlevel', edlevel, 'slct_level', [],...
                                                              'case_date', case_date, 'StepsInADay', StepsInADay(k), 'timeshift', [], 'select_field', mainfield_list,...
                                                              'NoOfWorkers', NumWorkers); 

                        else
                            error('Code not set for the data source %s!', data_source);
                        end
                        if savemat
                            S = dataset.(dataname{k});
                            disp(['Writing into ', loadfile, ' ...'])
                            save(loadfile, 'S','-v7.3');
                            clear S;
                        end
                    end
                    %/ Label it with '_his' or '_fut'
                    dataset.(strcat(dataname{k}, str_hisfut)) = dataset.(dataname{k});
                    dataset = rmfield(dataset, dataname{k});  %/ spare memory
                end

                %/==================== Data Processing ====================
                %/ Process dwdt, MFC, MFCv 
                if isequal(dataname{k}, strcat(data_source, '_q'))
                    q.lon   = dataset.(strcat(dataname{k}, str_hisfut)).lon;   %/ as a reference lon
                    q.lat   = dataset.(strcat(dataname{k}, str_hisfut)).lat;   %/ as a reference lat
                    q.level = dataset.(strcat(dataname{k}, str_hisfut)).level; %/ as a reference p-level
                    q.date_yyyymmdd_AllYr = dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr;

                    if t == 1
                        WBD.lon   = q.lon;
                        WBD.lat   = q.lat;
                        WBD.level = q.level;
                        level     = WBD.level;
                        fprintf('*** Detected min level = %d, max level = %d ***\n', min(level), max(level));
                        [~, lat_2D] = meshgrid(WBD.lon, WBD.lat);
                        lat_2D      = lat_2D';
                        hx          = diff(WBD.lon(1:2))*pi/180;   %/ Change degree to radian
                        hy          = diff(WBD.lat(1:2))*pi/180;   %/ Change degree to radian
                    end
                else
                    %/ First, make lon/lat dim consistent 
                    if ~isequal(dataset.(strcat(dataname{k}, str_hisfut)).lon, q.lon)...
                        && ~isequal(dataset.(strcat(dataname{k}, str_hisfut)).lat, q.lat)
                        fprintf('** Detected that %s has a different lon/lat with q, now making them consistent... ***\n', dataname{k})

                        %/ Subsetting
                        ind_lon = findismember_loop(dataset.(strcat(dataname{k}, str_hisfut)).lon, q.lon);
                        ind_lat = findismember_loop(dataset.(strcat(dataname{k}, str_hisfut)).lat, q.lat);

                        %/ Checking
                        if length(ind_lon) ~= length(q.lon) error('It seems impossible to make lon of %s consistent with q. Check your data!\n', dataname{k}); end
                        if length(ind_lat) ~= length(q.lat) error('It seems impossible to make lat of %s consistent with q. Check your data!\n', dataname{k}); end

                        if time_dim == 3
                            dataset.(strcat(dataname{k}, str_hisfut)).daily = dataset.(strcat(dataname{k}, str_hisfut)).daily(ind_lon, ind_lat, :); 
                        elseif time_dim == 4
                            dataset.(strcat(dataname{k}, str_hisfut)).daily = dataset.(strcat(dataname{k}, str_hisfut)).daily(ind_lon, ind_lat, :, :); 
                        end
                        dataset.(strcat(dataname{k}, str_hisfut)).lon = dataset.(strcat(dataname{k}, str_hisfut)).lon(ind_lon);
                        dataset.(strcat(dataname{k}, str_hisfut)).lat = dataset.(strcat(dataname{k}, str_hisfut)).lat(ind_lat);
                    end

                    %/ Prepare data
                    if isequal(dataname{k}, strcat(data_source,'_sp'))

                        %/ With q and sp loaded, we can compute -dwdt
                        ht    = 1;
                        [~, ~, ~, dqdt_daily] = gradient(dataset.(strcat(data_source, '_q', str_hisfut)).daily, ht);  %/ [kg/kg day-1]
    
                        %/ dwdt
                        dwdt_daily = -1/rho_w*vertinte('data', dqdt_daily, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                                       'sp_mode', sp_mode,   'sp',    dataset.(strcat(data_source, '_sp', str_hisfut)).daily, 'topo_lon', WBD.lon,  'topo_lat',   WBD.lat,...
                                                       'interp_mode', interp_mode, 'NumWorkers', [])*1000;  %/ from m/day to mm/day
                        clear dqdt_daily;
                        [dwdt.(strcat('yrly', str_hisfut)), ~, ~, ~] = daily2any('data_daily', dwdt_daily,...
                                                                                 'dates',      dataset.(strcat(data_source, '_q', str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);
                        clear dwdt_daily;

                    elseif isequal(dataname{k}, strcat(data_source, '_U'))
                        uq_daily = dataset.(strcat(data_source, '_U', str_hisfut)).daily.*dataset.(strcat(data_source, '_q', str_hisfut)).daily;

                    elseif isequal(dataname{k}, strcat(data_source, '_V'))
                        %/ With q, sp, U, V loaded, we can compute MFC
                        vq_daily = dataset.(strcat(data_source, '_V', str_hisfut)).daily.*dataset.(strcat(data_source, '_q', str_hisfut)).daily;
                        
                        %/ Moisture flux divergence
                        [   ~, dAdx, ~, ~] = gradient(uq_daily, hx);
                        [dAdy,    ~, ~, ~] = gradient(vq_daily, hy);
                        
                        div_A = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D));  %/ == du/dx + dv/dy in cartesian coor.
                        
                        %/ Vertically Integrated Moisture Flux Convergence
                        MFC_daily = -1/rho_w*vertinte('data', div_A, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                                       'sp_mode', sp_mode,   'sp',    dataset.(strcat(data_source, '_sp', str_hisfut)).daily, 'topo_lon', WBD.lon,  'topo_lat',   WBD.lat,...
                                                       'interp_mode', interp_mode, 'NumWorkers', [])*1000*3600*24;  %/ from m/s to mm/day
                        %/ Spare memory
                        clear dAdx;
                        clear dAdy;
                        clear div_A;
                        clear uq_daily;
                        clear vq_daily;
                        [MFC.(strcat('yrly', str_hisfut)), ~, ~, ~] = daily2any('data_daily', MFC_daily,...
                                                                                 'dates',      dataset.(strcat(data_source, '_q', str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);
                        clear MFC_daily;

                    elseif isequal(dataname{k}, strcat(data_source, '_W'))
                        wq_daily = dataset.(strcat(data_source, '_W', str_hisfut)).daily.*dataset.(strcat(data_source, '_q', str_hisfut)).daily;

                        dAdp  = gradient_Pcoor('data', wq_daily,  'level',   level, 'level_unit', level_unit, 'lnP_coor', lnP_coor, 'lev_dim', 3);
                        
                        %/ Vertically Integrated Moisture Flux Convergence
                        MFCv_daily = -1/rho_w*vertinte('data', dAdp, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                                       'sp_mode', sp_mode,   'sp',    dataset.(strcat(data_source, '_sp', str_hisfut)).daily, 'topo_lon', WBD.lon,  'topo_lat',   WBD.lat,...
                                                       'interp_mode', interp_mode, 'NumWorkers', [])*1000*3600*24;  %/ from m/s to mm/day
                        clear dAdp;
                        clear wq_daily;

                        [MFCv.(strcat('yrly', str_hisfut)), ~, ~, ~] = daily2any('data_daily', MFCv_daily,...
                                                                                 'dates',      dataset.(strcat(data_source, '_q', str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);
                        clear MFCv_daily;

                    elseif isequal(dataname{k}, strcat(data_source, '_P'))
                        [P.(strcat('yrly', str_hisfut)), ~, ~, ~]  = daily2any('data_daily', dataset.(strcat(dataname{k}, str_hisfut)).daily,...
                                                                               'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);
                        P.date_yyyymmdd_AllYr = dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr;
                    
                    elseif isequal(dataname{k}, strcat(data_source, '_E'))
                        [E.(strcat('yrly', str_hisfut)), ~, ~, ~]  = daily2any('data_daily', dataset.(strcat(dataname{k}, str_hisfut)).daily,...
                                                                               'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);
                        E.date_yyyymmdd_AllYr = dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr;

                    else
                        error('Code not set for %s!', dataname{k});                        
                    end
                end

                %/ Spare memory (except q and sp)
                if ~isequal(dataname{k}, strcat(data_source, '_q')) && ~isequal(dataname{k}, strcat(data_source, '_sp'))
                    dataset = rmfield(dataset, strcat(dataname{k}, str_hisfut));
                end
            end
        end
        %/ Spare memory by removing q from dataset
        dataset = rmfield(dataset, strcat(data_source, '_q', str_hisfut));
        dataset = rmfield(dataset, strcat(data_source, '_sp', str_hisfut));

        %/ Bootstrapping (if quiried)
        if ~isempty(nboot)
            fprintf('*** Detected non-empty ''nboot''. bootstrapping over samples in the two periods will be performed. ***\n')
            rng(nboot_seed)  %/ for reproducibility
            [~, bootsam_his] = bootstrp(nboot, @mean, his_yrs); %/ bootsam_his is an index matrix (nyear * nboot) 
            [~, bootsam_fut] = bootstrp(nboot, @mean, fut_yrs); %/ bootsam_fut is an index matrix (nyear * nboot) 
            WBD.nboot      = nboot;       %/ for record
            WBD.nboot_seed = nboot_seed;  %/ for record
        end
        
        %/ Reallocation
        if isequal(WBD_mode, 'bulk')
            WBD_terms       = {'dP', 'dE', 'MFC', 'ST', 'Res', 'dP_minus_dE'}; 
        elseif isequal(WBD_mode, 'bulk3D')
            WBD_terms       = {'dP', 'dE', 'MFC', 'MFCv', 'ST', 'Res', 'dP_minus_dE'}; 
        end
        for i = 1:length(WBD_terms)
            WBD.(WBD_terms{i}) = nan(length(WBD.lon), length(WBD.lat), nboot);
        end
        
        if isempty(gcp('nocreate')) && ~isempty(NumWorkers_HPC) %/ if set worker number
            parpool('Threads', NumWorkers_HPC) %/ use process-based parpool (threads-based parpool is only available after R2020a :((
        end
    
        %/ Compute Water Budget Decomposition (WBD)
        WBD_MFC         = cell(nboot, 1); 
        WBD_MFCv        = cell(nboot, 1);
        WBD_dP_minus_dE = cell(nboot, 1); 
        WBD_dP          = cell(nboot, 1); 
        WBD_dE          = cell(nboot, 1); 
        WBD_ST          = cell(nboot, 1); 
        tic;
        % for n = 1:nboot
        parfor n = 1:nboot
            if nboot > 1   
                fprintf('*** Bootstrapping (%d/%d)... ***\n', n, nboot);
            end

            %/ Broadcase variable
            dwdt_bc = dwdt;
            MFC_bc = MFC;
            P_bc = P;
            E_bc = E;
            
            %/ Mean of yearly data in the historical period
            dwdt_his  = mean(dwdt_bc.yrly_his(:,:,bootsam_his(:,n)),    3, 'omitnan');
            MFC_his   = mean(MFC_bc.yrly_his(:,:,bootsam_his(:,n)),     3, 'omitnan');
            P_his     = mean(P_bc.yrly_his(:,:,bootsam_his(:,n)),       3, 'omitnan');
            E_his     = mean(E_bc.yrly_his(:,:,bootsam_his(:,n)),       3, 'omitnan');

            %/ Mean of yearly data in the historical period
            dwdt_fut  = mean(dwdt_bc.yrly_fut(:,:,bootsam_fut(:,n)),    3, 'omitnan');
            MFC_fut   = mean(MFC_bc.yrly_fut(:,:,bootsam_fut(:,n)),     3, 'omitnan');
            P_fut     = mean(P_bc.yrly_fut(:,:,bootsam_fut(:,n)),       3, 'omitnan');
            E_fut     = mean(E_bc.yrly_fut(:,:,bootsam_fut(:,n)),       3, 'omitnan');

            %/ Take difference in the period mean
            delta_dwdt  = dwdt_fut - dwdt_his;
            delta_MFC   = MFC_fut - MFC_his;
            delta_P     = P_fut  - P_his;
            delta_E     = E_fut  - E_his;
            
            WBD_ST{n}  = delta_dwdt;
            WBD_MFC{n} = delta_MFC;
            WBD_dP_minus_dE{n} = delta_P - delta_E;
            WBD_dP{n} = delta_P;
            WBD_dE{n} = delta_E;

            if isequal(WBD_mode, 'bulk3D')
                MFCv_bc  = MFCv;
                MFCv_his = mean(MFCv_bc.yrly_his(:,:,bootsam_his(:,n)),     3, 'omitnan');
                MFCv_fut = mean(MFCv_bc.yrly_fut(:,:,bootsam_fut(:,n)),     3, 'omitnan');
                delta_MFCv = MFCv_fut - MFCv_his;
                WBD_MFCv{n} = delta_MFCv;
            end
        end
        toc;

        %/ Concatenate
        if ~isempty(WBD_MFC)           WBD.('MFC') = cat(3, WBD_MFC{:});   end
        if ~isempty(WBD_MFCv)          WBD.('MFCv') = cat(3, WBD_MFCv{:}); end
        if ~isempty(WBD_dP_minus_dE)   WBD.('dP_minus_dE') = cat(3, WBD_dP_minus_dE{:});   end
        if ~isempty(WBD_dP)            WBD.('dP') = cat(3, WBD_dP{:});     end
        if ~isempty(WBD_dE)            WBD.('dE') = cat(3, WBD_dE{:});     end
        if ~isempty(WBD_ST)            WBD.('ST') = cat(3, WBD_ST{:});     end
        WBD_MFC = []; WBD_MFCv = []; WBD_dP_minus_dE = []; WBD_dP = []; WBD_dE = []; WBD_ST = [];

        %/ Residual (Res) term
        if isequal(WBD_mode, 'bulk')
            WBD.('Res') = WBD.('dP_minus_dE') - WBD.('MFC') - WBD.('ST');
        elseif isequal(WBD_mode, 'bulk3D')
            WBD.('Res') = WBD.('dP_minus_dE') - WBD.('MFC') - WBD.('MFCv') - WBD.('ST');
        end
        
%         %/ Check Res
%         max_res       = max(abs(WBD.('Res')), [], 'all');
%         max_res_prct  = max(abs(WBD.('Res')./WBD.('dP')), [], 'all');
%         mean_res      = mean(abs(WBD.('Res')), 'all');
%         mean_res_prct = mean(abs(WBD.('Res')./WBD.('dP')), 'all');
%         fprintf('*** Max. abs. residual = %.2f mm/day ***\n', max_res)
%         fprintf('*** Max. abs. ratio of residual to dP = %.2f ***\n', max_res_prct)
%         fprintf('*** Mean abs. residual = %.2f mm/day ***\n', mean_res)
%         fprintf('*** Mean abs. ratio of residual to dP = %.2f ***\n', mean_res_prct)
            
        if savemat
           save(WBD_matname, 'WBD', '-v7.3');
           fprintf('!!! Saved WBD data into %s. !!!\n', WBD_matname);
        end
        
        fprintf('!!! Water budget decomposition completed !!!\n')
    end
    
    %/ Close parpool when not using
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        delete(poolobj)
    end

end