%%
function WBD = WBD_v2_wrapper(varargin)

    pnames = {'data_source', 'his_period', 'fut_period', 'dataset', 'mth', 'str_mth', 'stlevel', 'edlevel', 'level_unit',...
              'select_field',  'ins_or_acc', 'sp_mode', 'lnP_coor', 'data_folder', 'savemat', ...
              'recompute_yrly_data', 'recompute_WBD', 'NumWorkers'}; 
    dflts  = cell(length(pnames), 1);
    [          data_source,   his_period,   fut_period,   dataset,   mth,   str_mth,   stlevel,   edlevel,   level_unit,...
               select_field,    ins_or_acc,   sp_mode,   lnP_coor,   data_folder,   savemat,...
               recompute_yrly_data,   recompute_WBD,   NumWorkers] = internal.stats.parseArgs(pnames, dflts, varargin{:});
%%
    %/ Author:               Fandy Cheng
    %/ Date of modification: 15 May 2023
    
    %/ Description: The function ouputs decomposed water budget following
    %/              Dai et al. (2022)
    %/
    %/ NOTE:    This is the version 2.
    error('code deprecated')
    if mth == 0
        str_mth_bc = '';
    else
        str_mth_bc = strcat('_',str_mth{mth});
    end
    str_WBD_remark = 'complete';
    domain_str = '_global'; 
    if isempty(dataset)     dataset.placefolder = [];   end
    if lnP_coor  str_lnP_coor  = '_lnP';  else  str_lnP_coor  = '';   end 
    if isequal(data_source, 'CERA-20C')   %/ Use read_CERA
        mode = 1;
        addpath(genpath('/home/tfchengac/Flexpart_Codes'));
        run('param_flexpart'); %/ load dataname, data_path_parts, etc.
        
        if isempty(level_unit)  level_unit = 'hPa';         end
        str_plev = sprintf('%d_%d',stlevel, edlevel);
        var_list = [strcat({'q', 'U', 'V', 'W'}, str_plev), 'sp', 'P', 'E'];  %/ load q before all others
        
    else
        error('Code not set for the data source %s!', data_source);
    end
    his_yrs = his_period(1):his_period(2);
    fut_yrs = fut_period(1):fut_period(2);
    str_fut_minus_his = sprintf('%d-%d_minus_%d-%d',fut_period(1),fut_period(2),his_period(1),his_period(2));

    matname = sprintf('WBD_v2_%s_%s%s%s_%s.mat', str_WBD_remark, str_fut_minus_his, str_mth_bc, str_lnP_coor, sp_mode);    
    matname = strcat(data_folder, strrep(matname, ' ', '_'));
    
    if recompute_WBD == 0 && isfile(matname)
        fprintf('!!! mat file is found: %s. Loading from it... !!! \n', matname)
        load(matname, 'WBD');
    else
        %============== Load 4D/3D data for historical period ================%
        %/ Period Loop 
        q = []; dqdt = []; u = []; v = []; w = []; P = []; E = []; sp = [];
        TE_qu = []; TE_qv = []; TE_qw = [];
        
        yrly_data_matname = sprintf('yrly_data_for_WBD_%s%s.mat', str_fut_minus_his, str_mth_bc);    
        yrly_data_matname = strcat(data_folder, strrep(yrly_data_matname, ' ', '_'));
        if recompute_yrly_data == 0 && isfile(yrly_data_matname)
            fprintf('!!! mat file is found: %s. Loading from it... !!! \n', yrly_data_matname)
            load(yrly_data_matname, 'q', 'dqdt', 'u', 'v', 'w', 'P', 'E', 'sp', 'TE_qu', 'TE_qv', 'TE_qw');
        else
            for pd = 1:2
                if pd == 1
                    yrs        = his_yrs;
                    str_hisfut = '_his';
                elseif pd == 2
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
                case_date   = date_array_gen('years', yrs, 'season', season, 'st_month', st_month, 'ed_month', ed_month, 'dt_slct_hr', 24, 'output_date_format', 'yyyymmdd');
                select_data = findismember_loop(dataname, var_list);  
                dataname(select_data)
                if ~isequal(length(select_data), length(var_list))  error('Check your dataname!'); end
                for k = select_data
                    disp(strcat(dataname{k}, str_hisfut))
                    if ismember(dataname{k}, {'sp', 'P', 'E'})
                        time_dim        = 3;
                        mainfield_list  = [select_field(:)', {'lon'}, {'lat'}, {'date_yyyymmdd_AllYr'}];
                    else
                        time_dim        = 4;
                        mainfield_list  = [select_field(:)', {'lon'}, {'lat'}, {'level'}, {'date_yyyymmdd_AllYr'}];
                    end
                    loadfile = strcat(data_folder, dataname{k}, '_', select_field{:}, domain_str, '_', str_yrs, '.mat');

                    %/ 1. Check if the field has been loaded, 
                    %/ 2. Check if it has been processed, 
                    %/ 3. Otherwise will do data processing.
                    if isfield(dataset, strcat(dataname{k}, str_hisfut))
                        disp(['The quried data has been loaded in ''dataset''.']);
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
                            if isequal(data_source, 'CERA-20C')   %/ Use read_CERA
                                dataset.(dataname{k}) = read_CERA('mode', mode, 'time_dim', time_dim, ...
                                                                  'datafolder',   data_path_parts{whichfolder(k),1}, 'filename_prefix', data_path_parts{whichfolder(k),2},...
                                                                  'filename_part',data_path_parts{whichfolder(k),3}, ...
                                                                  'dataname', dataname{k},'dataname_callfile', dataname_callfile{k}, 'varname', varname{k},...
                                                                  'dataunitconv',dataunitconv(k),'datatype',datatype{k}, 'allyears', yrs,...
                                                                  'stlon', stlon,'edlon', edlon, 'stlat', stlat, 'edlat', edlat,...
                                                                  'stlevel', stlevel,'edlevel', edlevel, 'slct_level', [],...
                                                                  'case_date', case_date, 'StepsInADay', StepsInADay(k), 'timeshift', timeshift, 'selectfield',mainfield_list,...
                                                                  'fc_steps', fc_steps{k}, 'NumWorkers', NumWorkers);
                            else
                                error('Code not set for the data source %s!', data_source);
                            end
                            if savemat
                                S = dataset.(dataname{k});
                                filename = strcat(data_folder,dataname{k},'_',mainfield_list{1}, domain_str, '_', str_yrs, '.mat');
                                disp(['Writing into ', filename, ' ...'])
                                save(filename, 'S','-v7.3');
                                clear S;
                            end
                        end
                        %/ Label it with '_his' or '_fut'
                        dataset.(strcat(dataname{k}, str_hisfut)) = dataset.(dataname{k});
                        dataset = rmfield(dataset, dataname{k});  %/ spare memory
                    end

                    %/ Process q, dqdt, UVW, P, E, sp
                    if isequal(dataname{k}, strcat('q', str_plev))
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
                        if isequal(dataname{k}, strcat('U', str_plev))
                            [u.(strcat('yrly', str_hisfut)), ~, ~, u.(strcat('prime_daily', str_hisfut))] = daily2any('data_daily', dataset.(strcat(dataname{k}, str_hisfut)).daily,...
                                                                                                                      'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);

                            X = q.(strcat('prime_daily', str_hisfut)).*u.(strcat('prime_daily', str_hisfut));   %/ Transient Eddies (TE) q'u'
                            [TE_qu.(strcat('yrly', str_hisfut)), ~, ~, ~] = daily2any('data_daily', X,...
                                                                                      'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);
                            u.date_yyyymmdd_AllYr = dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr;

                            %/ Spare memory
                            u = rmfield(u, strcat('prime_daily', str_hisfut));      
                            clear X;                                                  

                        elseif isequal(dataname{k}, strcat('V', str_plev))
                            [v.(strcat('yrly', str_hisfut)), ~, ~, v.(strcat('prime_daily', str_hisfut))] = daily2any('data_daily', dataset.(strcat(dataname{k}, str_hisfut)).daily,...
                                                                                                                      'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);

                            X = q.(strcat('prime_daily', str_hisfut)).*v.(strcat('prime_daily', str_hisfut));   %/ Transient Eddies (TE) q'v'
                            [TE_qv.(strcat('yrly', str_hisfut)), ~, ~, ~] = daily2any('data_daily', X,...
                                                                                      'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);
                            v.date_yyyymmdd_AllYr = dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr;

                            %/ Spare memory                         
                            v = rmfield(v, strcat('prime_daily', str_hisfut));       
                            clear X;   

                        elseif isequal(dataname{k}, strcat('W', str_plev))
                            [w.(strcat('yrly', str_hisfut)), ~, ~, w.(strcat('prime_daily', str_hisfut))] = daily2any('data_daily', dataset.(strcat(dataname{k}, str_hisfut)).daily,...
                                                                                                                      'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);

                            X = q.(strcat('prime_daily', str_hisfut)).*w.(strcat('prime_daily', str_hisfut));   %/ Transient Eddies (TE) q'w'
                            [TE_qw.(strcat('yrly', str_hisfut)), ~, ~, ~] = daily2any('data_daily', X,...
                                                                                      'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);
                            w.date_yyyymmdd_AllYr = dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr;

                            %/ Spare memory     
                            w       = rmfield(w, strcat('prime_daily', str_hisfut));      
                            clear X;   

                        elseif isequal(dataname{k}, 'P')
                            [P.(strcat('yrly', str_hisfut)), ~, ~, ~]  = daily2any('data_daily', dataset.(strcat(dataname{k}, str_hisfut)).daily,...
                                                                                   'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);
                            P.date_yyyymmdd_AllYr = dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr;
                        elseif isequal(dataname{k}, 'E')
                            [E.(strcat('yrly', str_hisfut)), ~, ~, ~]  = daily2any('data_daily', dataset.(strcat(dataname{k}, str_hisfut)).daily,...
                                                                                   'dates',      dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr, 'mth', mth, 'ins_or_acc', ins_or_acc);
                            E.date_yyyymmdd_AllYr = dataset.(strcat(dataname{k}, str_hisfut)).date_yyyymmdd_AllYr;
                        elseif isequal(dataname{k}, 'sp')
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
                %/ Finally, remove q prime to save memory
                q = rmfield(q, strcat('prime_daily', str_hisfut));
            end
            if savemat
                save(yrly_data_matname, 'q', 'dqdt', 'u', 'v', 'w', 'P', 'E', 'sp', 'TE_qu', 'TE_qv', 'TE_qw', '-v7.3');
            end
        end
        %/ Historical period
        q_his     = mean(q.yrly_his,    4, 'omitnan');    %/ The lower level data may contain NaN (likely due to levels below topo in some weather conditions)
        u_his     = mean(u.yrly_his,    4, 'omitnan');
        v_his     = mean(v.yrly_his,    4, 'omitnan');
        w_his     = mean(w.yrly_his,    4, 'omitnan');
        dqdt_his  = mean(dqdt.yrly_his, 4, 'omitnan');
        P_his     = mean(P.yrly_his,    3, 'omitnan');
        E_his     = mean(E.yrly_his,    3, 'omitnan');
        TE_qu_his = mean(TE_qu.yrly_his, 4, 'omitnan');   %/ Transient Eddies (TE) q'u'
        TE_qv_his = mean(TE_qv.yrly_his, 4, 'omitnan');   %/ Transient Eddies (TE) q'v'
        TE_qw_his = mean(TE_qw.yrly_his, 4, 'omitnan');   %/ Transient Eddies (TE) q'w'
        sp_his    = mean(sp.yrly_his, 3, 'omitnan');      %/ Following Moon and Ha (2020)
        
        %/ Future period
        q_fut     = mean(q.yrly_fut,    4, 'omitnan');    %/ the lower level data may contain NaN (likely due to levels below topo in some weather conditions)
        u_fut     = mean(u.yrly_fut,    4, 'omitnan');
        v_fut     = mean(v.yrly_fut,    4, 'omitnan');
        w_fut     = mean(w.yrly_fut,    4, 'omitnan');
        dqdt_fut  = mean(dqdt.yrly_fut, 4, 'omitnan');
        P_fut     = mean(P.yrly_fut,    3, 'omitnan');
        E_fut     = mean(E.yrly_fut,    3, 'omitnan');
        TE_qu_fut = mean(TE_qu.yrly_fut, 4, 'omitnan');   %/ Transient Eddies (TE) q'u'
        TE_qv_fut = mean(TE_qv.yrly_fut, 4, 'omitnan');   %/ Transient Eddies (TE) q'v'
        TE_qw_fut = mean(TE_qw.yrly_fut, 4, 'omitnan');   %/ Transient Eddies (TE) q'w'
        
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

        %/ Compute Water Budget Decomposition (WBD)
        r       = 6371e3; 
        wd      = 997;                       %/ water density (kg m^-3)
        WBD     = [];
        WBD.lon = q.lon;
        WBD.lat = q.lat;
        level   = q.level;
        [lon_D, lat_2D] = meshgrid(WBD.lon, WBD.lat);
        lon_D = lon_D'; lat_2D = lat_2D';
        hx = diff(WBD.lon(1:2))*pi/180;   %/ rmb to change degree to radian!!!!
        hy = diff(WBD.lat(1:2))*pi/180;   %/ rmb to change degree to radian!!!!
        
        %/ THh term
%         [   ~, dAdx, ~] = gradient(u_his.*delta_q, hx);
%         [dAdy,    ~, ~] = gradient(v_his.*delta_q, hy);
%         div_A           = 1./(r*cosd(lat_q_2D)).*(dAdx + dAdy.*cosd(lat_q_2D));  %/ == du/dx + dv/dy in cartesian coor.
%         WBD.('THh')  = -1/wd*vertinte('data', div_A, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
%                                      'sp_mode', sp_mode,  'sp',   sp_his,   'topo_lon', WBD.lon, 'topo_lat',   WBD.lat)*1000*3600*24;  %/ from m/s to mm/day

        %/ THh term
        uq_vertinte  = -1/wd*vertinte('data', u_his.*delta_q, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                     'sp_mode', sp_mode,  'sp',   sp_his,   'topo_lon', WBD.lon, 'topo_lat',   WBD.lat)*1000*3600*24;  %/ from m/s to mm/day
        vq_vertinte  = -1/wd*vertinte('data', v_his.*delta_q, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                     'sp_mode', sp_mode,  'sp',   sp_his,   'topo_lon', WBD.lon, 'topo_lat',   WBD.lat)*1000*3600*24;  %/ from m/s to mm/day
        [   ~, dAdx] = gradient(uq_vertinte, hx);
        [dAdy,    ~] = gradient(vq_vertinte, hy);
        WBD.('THh')     = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D));  %/ == du/dx + dv/dy in cartesian coor.
        
        
        %/ THv term
        dAdp        = gradient_Pcoor('data',  w_his.*delta_q,  'level',   level, 'level_unit', level_unit, 'lnP_coor', lnP_coor, 'lev_dim', 3);
        WBD.('THv') = -1/wd*vertinte('data',  dAdp, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                    'sp_mode', sp_mode,  'sp',    sp_his,   'topo_lon', WBD.lon,   'topo_lat',   WBD.lat)*1000*3600*24;  %/ from m/s to mm/day
        
        %/ DYh term
%         [   ~, dAdx, ~] = gradient(q_his.*delta_u, hx);
%         [dAdy,    ~, ~] = gradient(q_his.*delta_v, hy);
%         div_A           = 1./(r*cosd(lat_q_2D)).*(dAdx + dAdy.*cosd(lat_q_2D));  %/ == du/dx + dv/dy in cartesian coor.
%         WBD.('DYh')  = -1/wd*vertinte('data',  div_A, 'level',      level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
%                                      'sp_mode', sp_mode,  'sp',    sp_his,    'topo_lon',   WBD.lon,   'topo_lat',   WBD.lat)*1000*3600*24;  %/ from m/s to mm/day

        uq_vertinte  = -1/wd*vertinte('data', delta_u.*q_his, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                      'sp_mode', sp_mode,  'sp',   sp_his,   'topo_lon', WBD.lon, 'topo_lat',   WBD.lat)*1000*3600*24;  %/ from m/s to mm/day
        vq_vertinte  = -1/wd*vertinte('data', delta_v.*q_his, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                      'sp_mode', sp_mode,  'sp',   sp_his,   'topo_lon', WBD.lon, 'topo_lat',   WBD.lat)*1000*3600*24;  %/ from m/s to mm/day
        [   ~, dAdx] = gradient(uq_vertinte, hx);
        [dAdy,    ~] = gradient(vq_vertinte, hy);
        WBD.('DYh')     = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D));  %/ == du/dx + dv/dy in cartesian coor.
        
        
        %/ DYv term
        dAdp         = gradient_Pcoor('data', delta_w.*q_his,  'level',   level, 'level_unit', level_unit, 'lnP_coor', lnP_coor, 'lev_dim', 3);
        WBD.('DYv')  = -1/wd*vertinte('data', dAdp, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                     'sp_mode', sp_mode,  'sp',    sp_his,    'topo_lon', WBD.lon,   'topo_lat',   WBD.lat)*1000*3600*24;  %/ from m/s to mm/day

        %/ NL term
%         [   ~, dAdx, ~] = gradient(delta_q.*delta_u, hx);
%         [dAdy,    ~, ~] = gradient(delta_q.*delta_v, hy);
%         dAdp            = gradient_Pcoor('data', delta_q.*delta_w,  'level',   level, 'level_unit', level_unit, 'lnP_coor', lnP_coor, 'lev_dim', 3);
%         div_A           = 1./(r*cosd(lat_q_2D)).*(dAdx + dAdy.*cosd(lat_q_2D)) + dAdp;     %/ == du/dx + dv/dy  + dw/dp
%         WBD.('NL')  = -1/wd*vertinte('data', div_A, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
%                                     'sp_mode', sp_mode,   'sp',    sp_his,   'topo_lon', WBD.lon,   'topo_lat',   WBD.lat)*1000*3600*24;  %/ from m/s to mm/day

        uq_vertinte  = -1/wd*vertinte('data', delta_u.*delta_q, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                      'sp_mode', sp_mode,  'sp',   sp_his,   'topo_lon', WBD.lon, 'topo_lat',   WBD.lat)*1000*3600*24;  %/ from m/s to mm/day
        vq_vertinte  = -1/wd*vertinte('data', delta_v.*delta_q, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                      'sp_mode', sp_mode,  'sp',   sp_his,   'topo_lon', WBD.lon, 'topo_lat',   WBD.lat)*1000*3600*24;  %/ from m/s to mm/day
        
        [   ~, dAdx] = gradient(uq_vertinte, hx);
        [dAdy,    ~] = gradient(vq_vertinte, hy);
        dAdp            = gradient_Pcoor('data', delta_w.*delta_q, 'level', level, 'level_unit', level_unit, 'lnP_coor', lnP_coor, 'lev_dim', 3);
        VertTerm        = -1/wd*vertinte('data', dAdp, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                      'sp_mode', sp_mode,  'sp',   sp_his,   'topo_lon', WBD.lon, 'topo_lat',   WBD.lat)*1000*3600*24;  %/ from m/s to mm/day
        WBD.('NL')      = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D)) + VertTerm;  %/ == du/dx + dv/dy in cartesian coor.
        
        %/ TE term
%         [   ~, dAdx, ~] = gradient(delta_TE_qu, hx);
%         [dAdy,    ~, ~] = gradient(delta_TE_qv, hy);
%         dAdp            = gradient_Pcoor('data', delta_TE_qw, 'level', level, 'level_unit', level_unit, 'lnP_coor', lnP_coor, 'lev_dim', 3);
%         div_A           = 1./(r*cosd(lat_q_2D)).*(dAdx + dAdy.*cosd(lat_q_2D)) + dAdp;     %/ == du/dx + dv/dy  + dw/dp
%         WBD.('TE')  = -1/wd*vertinte('data', div_A, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
%                                     'sp_mode', sp_mode,   'sp',    sp_his,    'topo_lon', WBD.lon,  'topo_lat',   WBD.lat)*1000*3600*24;  %/ from m/s to mm/day

        uq_vertinte  = -1/wd*vertinte('data', delta_TE_qu, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                      'sp_mode', sp_mode,  'sp',   sp_his,   'topo_lon', WBD.lon, 'topo_lat',   WBD.lat)*1000*3600*24;  %/ from m/s to mm/day
        vq_vertinte  = -1/wd*vertinte('data', delta_TE_qv, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                      'sp_mode', sp_mode,  'sp',   sp_his,   'topo_lon', WBD.lon, 'topo_lat',   WBD.lat)*1000*3600*24;  %/ from m/s to mm/day
        dAdp            = gradient_Pcoor('data', delta_TE_qw, 'level', level, 'level_unit', level_unit, 'lnP_coor', lnP_coor, 'lev_dim', 3);
        VertTerm        = -1/wd*vertinte('data', dAdp, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                      'sp_mode', sp_mode,  'sp',   sp_his,   'topo_lon', WBD.lon, 'topo_lat',   WBD.lat)*1000*3600*24;  %/ from m/s to mm/day
        [   ~, dAdx] = gradient(uq_vertinte, hx);
        [dAdy,    ~] = gradient(vq_vertinte, hy);
        WBD.('TE')      = 1./(r*cosd(lat_2D)).*(dAdx + dAdy.*cosd(lat_2D)) + VertTerm;  %/ == du/dx + dv/dy in cartesian coor.
        
        
        %/ P-E term
        WBD.('dP_minus_dE') = delta_P - delta_E;
        WBD.('dP') = delta_P;
        WBD.('dE') = delta_E;

        %/ Storage (ST) term
        WBD.('ST')  = -1/wd*vertinte('data', delta_dqdt, 'level',    level, 'level_unit', level_unit, 'lnP_coor', lnP_coor,...
                                     'sp_mode', sp_mode,   'sp',    sp_his,   'topo_lon', WBD.lon,  'topo_lat',   WBD.lat)*1000;  %/ from m/day to mm/day

        %/ Residual (Res) term
        WBD.('Res') = WBD.('dP_minus_dE') - WBD.('THh') - WBD.('THv') - WBD.('DYh') - WBD.('DYv') ...
                     - WBD.('TE') - WBD.('NL') - WBD.('ST');
        
        max_res       = max(abs(WBD.('Res')), [], 'all');
        max_res_prct  = max(abs(WBD.('Res')./WBD.('dP')), [], 'all');
        mean_res      = mean(abs(WBD.('Res')), 'all');
        mean_res_prct = mean(abs(WBD.('Res')./WBD.('dP')), 'all');
        fprintf('*** Max. abs. residual = %.2f mm/day ***\n', max_res)
        fprintf('*** Max. abs. ratio of residual to dP = %.2f ***\n', max_res_prct)
        fprintf('*** Mean abs. residual = %.2f mm/day ***\n', mean_res)
        fprintf('*** Mean abs. ratio of residual to dP = %.2f ***\n', mean_res_prct)
        if savemat
           save(matname, 'WBD', '-v7.3');
           fprintf('!!! Saved WBD data into %s. !!!\n', matname);
        end
    end


end