 %%
function WBD = compute_WBD(varargin)

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
    %/  Description: The function ouputs the fully decomposed water budget 
    %/               following Dai et al. (2022).
    %/      
    %/               Dai, L., Cheng, T.F. & Lu, M. Anthropogenic warming disrupts 
    %/               intraseasonal monsoon stages and brings dry-get-wetter climate
    %/               in future East Asia. npj Clim Atmos Sci 5, 11 (2022). 
    %/               https://doi.org/10.1038/s41612-022-00235-9
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
    if isempty(WBD_mode)        
        WBD_mode = 'complete'; 
    elseif ~any(ismember(WBD_mode, {'complete', '2D', 'basic', 'basic3D'}))
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

    var_list            = strcat(data_source, '_', {'q', 'U', 'V', 'W', 'sp', 'P', 'E'});
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
        q = []; dqdt = []; u = []; v = []; w = []; P = []; E = []; sp = [];
        TE_qu = []; TE_qv = []; TE_qw = [];
        
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
                %/ 1. Check if the field has been loaded, 
                %/ 2. Check if it has been processed, 
                %/ 3. Otherwise will do data processing.
                if isfield(dataset, strcat(dataname{k}, str_hisfut))
                    disp('The quried data has been loaded in ''dataset''.');
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

                %/ Process q, dqdt, UVW, P, E, sp
                if isequal(dataname{k}, strcat(data_source, '_q'))
%                     if isequal(dataname{k}, strcat(data_source, '_q', str_plev))
                    q.lon   = dataset.(strcat(dataname{k}, str_hisfut)).lon;   %/ as a reference lon
                    q.lat   = dataset.(strcat(dataname{k}, str_hisfut)).lat;   %/ as a reference lat
                    q.level = dataset.(strcat(dataname{k}, str_hisfut)).level; %/ as a reference p-level

                    %/ Prepare data
                    ht    = 1;
                    [~, ~, ~, dqdt_daily] = gradient(dataset.(strcat(dataname{k}, str_hisfut)).daily, ht);  %/ [kg/kg day-1]

                    [q.(strcat('yrly', str_hisfut)), ~, ~, q.(strcat('prime_daily', str_hisfut))] = daily2any('data_daily', dataset.(strcat(dataname{k}, str_hisfut)).daily,...
                                                                                                              'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);
                    [dqdt.(strcat('yrly', str_hisfut)), ~, ~, ~] = daily2any('data_daily', dqdt_daily,...
                                                                             'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);
                    q.date_yyyymmdd_AllYr = dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr;

                    %/ Spare memory
                    clear dqdt_daily; 
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
                    if isequal(dataname{k}, strcat(data_source, '_U'))
%                         if isequal(dataname{k}, strcat('U', str_plev))
                        [u.(strcat('yrly', str_hisfut)), ~, ~, u.(strcat('prime_daily', str_hisfut))] = daily2any('data_daily', dataset.(strcat(dataname{k}, str_hisfut)).daily,...
                                                                                                                  'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);
                        
                        %/ Transient Eddies (TE) q'u'
                        X = q.(strcat('prime_daily', str_hisfut)).*u.(strcat('prime_daily', str_hisfut));   
                        [TE_qu.(strcat('yrly', str_hisfut)), ~, ~, ~] = daily2any('data_daily', X,...
                                                                                  'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);
                        u.date_yyyymmdd_AllYr = dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr;

                        %/ Spare memory
                        u = rmfield(u, strcat('prime_daily', str_hisfut));      
                        clear X;                                                  

                    elseif isequal(dataname{k}, strcat(data_source, '_V'))
                        [v.(strcat('yrly', str_hisfut)), ~, ~, v.(strcat('prime_daily', str_hisfut))] = daily2any('data_daily', dataset.(strcat(dataname{k}, str_hisfut)).daily,...
                                                                                                                  'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);

                        %/ Transient Eddies (TE) q'v'
                        X = q.(strcat('prime_daily', str_hisfut)).*v.(strcat('prime_daily', str_hisfut));   
                        [TE_qv.(strcat('yrly', str_hisfut)), ~, ~, ~] = daily2any('data_daily', X,...
                                                                                  'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);
                        v.date_yyyymmdd_AllYr = dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr;

                        %/ Spare memory                         
                        v = rmfield(v, strcat('prime_daily', str_hisfut));       
                        clear X;   

                    elseif isequal(dataname{k}, strcat(data_source, '_W'))
                        [w.(strcat('yrly', str_hisfut)), ~, ~, w.(strcat('prime_daily', str_hisfut))] = daily2any('data_daily', dataset.(strcat(dataname{k}, str_hisfut)).daily,...
                                                                                                                  'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);

                        %/ Transient Eddies (TE) q'w'
                        X = q.(strcat('prime_daily', str_hisfut)).*w.(strcat('prime_daily', str_hisfut));   
                        [TE_qw.(strcat('yrly', str_hisfut)), ~, ~, ~] = daily2any('data_daily', X,...
                                                                                  'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);
                        w.date_yyyymmdd_AllYr = dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr;

                        %/ Spare memory     
                        w = rmfield(w, strcat('prime_daily', str_hisfut));      
                        clear X;   

                    elseif isequal(dataname{k}, strcat(data_source, '_P'))
                        [P.(strcat('yrly', str_hisfut)), ~, ~, ~]  = daily2any('data_daily', dataset.(strcat(dataname{k}, str_hisfut)).daily,...
                                                                               'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);
                        P.date_yyyymmdd_AllYr = dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr;
                    elseif isequal(dataname{k}, strcat(data_source, '_E'))
                        [E.(strcat('yrly', str_hisfut)), ~, ~, ~]  = daily2any('data_daily', dataset.(strcat(dataname{k}, str_hisfut)).daily,...
                                                                               'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);
                        E.date_yyyymmdd_AllYr = dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr;
                    elseif isequal(dataname{k}, strcat(data_source, '_sp'))
                        [sp.(strcat('yrly', str_hisfut)), ~, ~, ~] = daily2any('data_daily', dataset.(strcat(dataname{k}, str_hisfut)).daily,...
                                                                               'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);
                        sp.date_yyyymmdd_AllYr = dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr;
                    else
                        error('Code not set for %s!', dataname{k});                        
                    end
                end
                %/ Spare memory
                dataset = rmfield(dataset, strcat(dataname{k}, str_hisfut));
            end
        end

        %/ Parameter Setting for Water Budget Decomposition (WBD)
        r         = 6371e3; 
        rho_w     = 997;                           %/ Water density (kg m^-3)
        WBD       = [];
        WBD.lon   = q.lon;
        WBD.lat   = q.lat;
        WBD.level = q.level;
        level     = WBD.level;
        WBD.data_source = data_source;
        WBD.his_yrs = his_yrs;
        WBD.fut_yrs = fut_yrs;
        fprintf('*** Detected min level = %d, max level = %d ***\n', min(level), max(level));
        [~, lat_2D] = meshgrid(WBD.lon, WBD.lat);
        lat_2D      = lat_2D';
        hx          = diff(WBD.lon(1:2))*pi/180;   %/ Change degree to radian
        hy          = diff(WBD.lat(1:2))*pi/180;   %/ Change degree to radian
        
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
        if isequal(WBD_mode, 'complete')
            WBD_terms       = {'dP', 'dE', 'HEh', 'HEv', 'DYh', 'DYv', 'TE', 'NL', 'ST', 'Res', 'dP_minus_dE'}; 
        elseif isequal(WBD_mode, '2D')
            WBD_terms       = {'dP', 'dE', 'HEh',  'DYh', 'TE', 'NL', 'ST', 'Res', 'dP_minus_dE'}; 
        elseif isequal(WBD_mode, 'basic')
            WBD_terms       = {'dP', 'dE', 'MC', 'ST', 'Res', 'dP_minus_dE'}; 
        elseif isequal(WBD_mode, 'basic3D')
            WBD_terms       = {'dP', 'dE', 'MC3D', 'ST', 'Res', 'dP_minus_dE'}; 
        end
        for i = 1:length(WBD_terms)
            WBD.(WBD_terms{i}) = nan(length(WBD.lon), length(WBD.lat), nboot);
        end
        
        if isempty(gcp('nocreate')) && ~isempty(NumWorkers_HPC) %/ if set worker number
            parpool('Threads', NumWorkers_HPC) %/ use process-based parpool (threads-based parpool is only available after R2020a :((
        end
    
        %/ Compute Water Budget Decomposition (WBD)
        WBD_MC          = cell(nboot, 1); 
        WBD_MC3D        = cell(nboot, 1);
        WBD_HEh         = cell(nboot, 1); 
        WBD_DYh         = cell(nboot, 1); 
        WBD_HEv         = cell(nboot, 1); 
        WBD_DYv         = cell(nboot, 1); 
        WBD_NL          = cell(nboot, 1); 
        WBD_TE          = cell(nboot, 1); 
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
            q_bc = q;
            u_bc = u;
            v_bc = v;
            w_bc = w;
            dqdt_bc = dqdt;
            P_bc = P;
            E_bc = E;
            TE_qu_bc = TE_qu;
            TE_qv_bc = TE_qv;
            TE_qw_bc = TE_qw;
            sp_bc = sp;
            WBD_bc = WBD;
            
            %/ Mean of yearly data in the historical period
            q_his     = mean(q_bc.yrly_his(:,:,:,bootsam_his(:,n)),     4, 'omitnan');   %/ The lower level data may contain NaN (likely due to levels below topo in some weather conditions)
            u_his     = mean(u_bc.yrly_his(:,:,:,bootsam_his(:,n)),     4, 'omitnan');
            v_his     = mean(v_bc.yrly_his(:,:,:,bootsam_his(:,n)),     4, 'omitnan');
            w_his     = mean(w_bc.yrly_his(:,:,:,bootsam_his(:,n)),     4, 'omitnan');
            dqdt_his  = mean(dqdt_bc.yrly_his(:,:,:,bootsam_his(:,n)),  4, 'omitnan');
            P_his     = mean(P_bc.yrly_his(:,:,bootsam_his(:,n)),       3, 'omitnan');
            E_his     = mean(E_bc.yrly_his(:,:,bootsam_his(:,n)),       3, 'omitnan');
            TE_qu_his = mean(TE_qu_bc.yrly_his(:,:,:,bootsam_his(:,n)), 4, 'omitnan');   %/ Transient Eddies (TE) q'u'
            TE_qv_his = mean(TE_qv_bc.yrly_his(:,:,:,bootsam_his(:,n)), 4, 'omitnan');   %/ Transient Eddies (TE) q'v'
            TE_qw_his = mean(TE_qw_bc.yrly_his(:,:,:,bootsam_his(:,n)), 4, 'omitnan');   %/ Transient Eddies (TE) q'w'
            sp_his    = mean(sp_bc.yrly_his(:,:,bootsam_his(:,n)),      3, 'omitnan');   %/ Following Moon and Ha (2020)

            %/ Mean of yearly data in the historical period
            q_fut     = mean(q_bc.yrly_fut(:,:,:,bootsam_fut(:,n)),     4, 'omitnan');   %/ The lower level data may contain NaN (likely due to levels below topo in some weather conditions)
            u_fut     = mean(u_bc.yrly_fut(:,:,:,bootsam_fut(:,n)),     4, 'omitnan');
            v_fut     = mean(v_bc.yrly_fut(:,:,:,bootsam_fut(:,n)),     4, 'omitnan');
            w_fut     = mean(w_bc.yrly_fut(:,:,:,bootsam_fut(:,n)),     4, 'omitnan');
            dqdt_fut  = mean(dqdt_bc.yrly_fut(:,:,:,bootsam_fut(:,n)),  4, 'omitnan');
            P_fut     = mean(P_bc.yrly_fut(:,:,bootsam_fut(:,n)),       3, 'omitnan');
            E_fut     = mean(E_bc.yrly_fut(:,:,bootsam_fut(:,n)),       3, 'omitnan');
            TE_qu_fut = mean(TE_qu_bc.yrly_fut(:,:,:,bootsam_fut(:,n)), 4, 'omitnan');   %/ Transient Eddies (TE) q'u'
            TE_qv_fut = mean(TE_qv_bc.yrly_fut(:,:,:,bootsam_fut(:,n)), 4, 'omitnan');   %/ Transient Eddies (TE) q'v'
            TE_qw_fut = mean(TE_qw_bc.yrly_fut(:,:,:,bootsam_fut(:,n)), 4, 'omitnan');   %/ Transient Eddies (TE) q'w'
            sp_fut    = mean(sp_bc.yrly_fut(:,:,bootsam_fut(:,n)),      3, 'omitnan');   

            %/ Take difference in the period mean
            delta_q     = q_fut  - q_his;
            delta_u     = u_fut  - u_his;
            delta_v     = v_fut  - v_his;
            delta_w     = w_fut  - w_his;
            delta_dqdt  = dqdt_fut - dqdt_his;
            delta_P     = P_fut  - P_his;
            delta_E     = E_fut  - E_his;
            delta_TE_qu = TE_qu_fut - TE_qu_his;
            delta_TE_qv = TE_qv_fut - TE_qv_his;
            delta_TE_qw = TE_qw_fut - TE_qw_his;

            delta_qu    = q_fut.*u_fut - q_his.*u_his;
            delta_qv    = q_fut.*v_fut - q_his.*v_his;
            delta_qw    = q_fut.*w_fut - q_his.*w_his;
            
            % delta_sp    = sp_fut - sp_his;
            % max(abs(delta_sp), [], 'all')
            % mean(abs(delta_sp), 'all')

            if ismember(WBD_mode, {'basic'})
                %/ Moisture convergence (MC) term
                [   ~, dAdx, ~] = gradient(delta_qu, hx);
                [dAdy,    ~, ~] = gradient(delta_qv, hy);
                
                div_A           = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D));  %/ == du/dx + dv/dy in cartesian coor.
                
                WBD_MC{n} = -1/rho_w*vertinte('data', div_A, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                             'sp_mode', sp_mode,  'sp',   sp_his,   'topo_lon', WBD_bc.lon, 'topo_lat',   WBD_bc.lat,...
                                             'interp_mode', interp_mode, 'NumWorkers', [])*1000*3600*24;  %/ from m/s to mm/day

            elseif ismember(WBD_mode, {'basic3D'})
                %/ Moisture convergence (MC) term
                [   ~, dAdx, ~] = gradient(delta_qu, hx);
                [dAdy,    ~, ~] = gradient(delta_qv, hy);
                dAdp            = gradient_Pcoor('data', delta_qw,  'level',   level, 'level_unit', level_unit, 'lnP_coor', lnP_coor, 'lev_dim', 3);

                div_A           = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D)) + dAdp;     %/ == du/dx + dv/dy  + dw/dp

                WBD_MC3D{n} = -1/rho_w*vertinte('data', div_A, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                             'sp_mode', sp_mode,  'sp',   sp_his,   'topo_lon', WBD_bc.lon, 'topo_lat',   WBD_bc.lat,...
                                             'interp_mode', interp_mode, 'NumWorkers', [])*1000*3600*24;  %/ from m/s to mm/day

            else
                %/ HEh term
                [   ~, dAdx, ~] = gradient(u_his.*delta_q, hx);
                [dAdy,    ~, ~] = gradient(v_his.*delta_q, hy);
                div_A           = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D));  %/ == du/dx + dv/dy in cartesian coor.

                WBD_HEh{n} = -1/rho_w*vertinte('data', div_A, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                             'sp_mode', sp_mode,  'sp',   sp_his,   'topo_lon', WBD_bc.lon, 'topo_lat',   WBD_bc.lat,...
                                             'interp_mode', interp_mode, 'NumWorkers', [])*1000*3600*24;  %/ from m/s to mm/day
    
                %/ DYh term
                [   ~, dAdx, ~] = gradient(q_his.*delta_u, hx);
                [dAdy,    ~, ~] = gradient(q_his.*delta_v, hy);
                div_A           = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D));  %/ == du/dx + dv/dy in cartesian coor.

                WBD_DYh{n}  = -1/rho_w*vertinte('data',  div_A, 'level',      level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                              'sp_mode', sp_mode,  'sp',    sp_his,    'topo_lon',   WBD_bc.lon,   'topo_lat',   WBD_bc.lat,...
                                              'interp_mode', interp_mode, 'NumWorkers', [])*1000*3600*24;  %/ from m/s to mm/day
                       
                %/ If WBD_mode is 'complete', we include the vertical terms (CAUTION: it may produce large residue!!)
                if isequal(WBD_mode, 'complete')
                    %/ HEv term
                    dAdp        = gradient_Pcoor('data', w_his.*delta_q,  'level',   level, 'level_unit', level_unit, 'lnP_coor', lnP_coor, 'lev_dim', 3);
                    WBD_HEv{n} = -1/rho_w*vertinte('data',  dAdp, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                                 'sp_mode', sp_mode,  'sp',    sp_his,   'topo_lon', WBD_bc.lon,   'topo_lat',   WBD_bc.lat,...
                                                 'interp_mode', interp_mode, 'NumWorkers', [])*1000*3600*24;  %/ from m/s to mm/day
                                             
                    %/ DYv term
                    dAdp         = gradient_Pcoor('data', q_his.*delta_w,  'level',   level, 'level_unit', level_unit, 'lnP_coor', lnP_coor, 'lev_dim', 3);
                    WBD_DYv{n}  = -1/rho_w*vertinte('data', dAdp, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                                  'sp_mode', sp_mode,  'sp',    sp_his,    'topo_lon', WBD_bc.lon,   'topo_lat',   WBD_bc.lat,...
                                                  'interp_mode', interp_mode, 'NumWorkers', [])*1000*3600*24;  %/ from m/s to mm/day
                                  
                    %/ NL term (3D)
                    [   ~, dAdx, ~] = gradient(delta_q.*delta_u, hx);
                    [dAdy,    ~, ~] = gradient(delta_q.*delta_v, hy);
                    dAdp            = gradient_Pcoor('data', delta_q.*delta_w,  'level',   level, 'level_unit', level_unit, 'lnP_coor', lnP_coor, 'lev_dim', 3);
                    div_A           = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D)) + dAdp;     %/ == du/dx + dv/dy  + dw/dp
                    WBD_NL{n}   = -1/rho_w*vertinte('data', div_A, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                                  'sp_mode', sp_mode,   'sp',    sp_his,   'topo_lon', WBD_bc.lon,   'topo_lat',   WBD_bc.lat,...
                                                  'interp_mode', interp_mode, 'NumWorkers', [])*1000*3600*24;  %/ from m/s to mm/day
                                             
                    %/ TE term (3D)
                    [   ~, dAdx, ~] = gradient(delta_TE_qu, hx);
                    [dAdy,    ~, ~] = gradient(delta_TE_qv, hy);
                    dAdp            = gradient_Pcoor('data', delta_TE_qw, 'level', level, 'level_unit', level_unit, 'lnP_coor', lnP_coor, 'lev_dim', 3);
                    div_A           = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D)) + dAdp;     %/ == du/dx + dv/dy  + dw/dp
                    WBD_TE{n}   = -1/rho_w*vertinte('data', div_A, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                                  'sp_mode', sp_mode,   'sp',    sp_his,    'topo_lon', WBD_bc.lon,  'topo_lat',   WBD_bc.lat,...
                                                  'interp_mode', interp_mode, 'NumWorkers', [])*1000*3600*24;  %/ from m/s to mm/day
                                              
                elseif isequal(WBD_mode, '2D')
                    %/ NL term (2D)
                    [   ~, dAdx, ~] = gradient(delta_q.*delta_u, hx);
                    [dAdy,    ~, ~] = gradient(delta_q.*delta_v, hy);
                    div_A           = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D));     %/ == du/dx + dv/dy  + dw/dp
                    WBD_NL{n}   = -1/rho_w*vertinte('data', div_A, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                                  'sp_mode', sp_mode,   'sp',    sp_his,   'topo_lon', WBD_bc.lon,   'topo_lat',   WBD_bc.lat,...
                                                  'interp_mode', interp_mode, 'NumWorkers', [])*1000*3600*24;  %/ from m/s to mm/day

                    %/ TE term (2D)
                    [   ~, dAdx, ~] = gradient(delta_TE_qu, hx);
                    [dAdy,    ~, ~] = gradient(delta_TE_qv, hy);
                    div_A           = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D));     %/ == du/dx + dv/dy  + dw/dp
                    WBD_TE{n}   = -1/rho_w*vertinte('data', div_A, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                                  'sp_mode', sp_mode,   'sp',    sp_his,    'topo_lon', WBD_bc.lon,  'topo_lat',   WBD_bc.lat,...
                                                  'interp_mode', interp_mode, 'NumWorkers', [])*1000*3600*24;  %/ from m/s to mm/day
                end      
            end
            WBD_dP_minus_dE{n} = delta_P - delta_E;
            WBD_dP{n} = delta_P;
            WBD_dE{n} = delta_E;

            %/ Storage (ST) term
            WBD_ST{n}   = -1/rho_w*vertinte('data', delta_dqdt, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                          'sp_mode', sp_mode,   'sp',    sp_his,   'topo_lon', WBD_bc.lon,  'topo_lat',   WBD_bc.lat,...
                                          'interp_mode', interp_mode, 'NumWorkers', [])*1000;  %/ from m/day to mm/day
        end
        toc;

        %/ Concatenate
        if ~isempty(WBD_MC)            WBD.('MC') = cat(3, WBD_MC{:});     end
        if ~isempty(WBD_MC3D)          WBD.('MC3D') = cat(3, WBD_MC3D{:}); end
        if ~isempty(WBD_HEh)           WBD.('HEh') = cat(3, WBD_HEh{:});   end
        if ~isempty(WBD_DYh)           WBD.('DYh') = cat(3, WBD_DYh{:});   end
        if ~isempty(WBD_HEv)           WBD.('HEv') = cat(3, WBD_HEv{:});   end
        if ~isempty(WBD_DYv)           WBD.('DYv') = cat(3, WBD_DYv{:});   end
        if ~isempty(WBD_NL)            WBD.('NL') = cat(3, WBD_NL{:});     end
        if ~isempty(WBD_TE)            WBD.('TE') = cat(3, WBD_TE{:});     end
        if ~isempty(WBD_dP_minus_dE)   WBD.('dP_minus_dE') = cat(3, WBD_dP_minus_dE{:});   end
        if ~isempty(WBD_dP)            WBD.('dP') = cat(3, WBD_dP{:});     end
        if ~isempty(WBD_dE)            WBD.('dE') = cat(3, WBD_dE{:});     end
        if ~isempty(WBD_ST)            WBD.('ST') = cat(3, WBD_ST{:});     end
        WBD_MC = []; WBD_HEh = []; WBD_DYh = []; WBD_HEv = []; WBD_DYv = [];
        WBD_NL = []; WBD_TE = []; WBD_dP_minus_dE = []; WBD_dP = []; WBD_dE = []; WBD_ST = [];

        %/ Residual (Res) term
        if isequal(WBD_mode, 'complete')
            WBD.('Res') = WBD.('dP_minus_dE') - WBD.('HEh') - WBD.('HEv') - WBD.('DYh') - WBD.('DYv') ...
                         - WBD.('TE') - WBD.('NL') - WBD.('ST');
                     
        elseif isequal(WBD_mode, '2D')
            WBD.('Res') = WBD.('dP_minus_dE') - WBD.('HEh') - WBD.('DYh') - WBD.('TE') - WBD.('NL') - WBD.('ST');

        elseif isequal(WBD_mode, 'basic')
            WBD.('Res') = WBD.('dP_minus_dE') - WBD.('MC') - WBD.('ST');

        elseif isequal(WBD_mode, 'basic3D')
            WBD.('Res') = WBD.('dP_minus_dE') - WBD.('MC3D') - WBD.('ST');
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