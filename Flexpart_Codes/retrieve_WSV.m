function [WSV_data_daily, WSV_dates, rr, traj_catalog] = retrieve_WSV(varargin)

%/ create a set of valid parameters and their default value
pnames = {'project_name', 'WSV_name', 'WSV_dir', 'year_list', 'mth', 'st_month', 'st_day', 'ed_month', 'ed_day', 'str_RHc_dqc', 'ldirect', 'from_basin', 'maxtraj_day',...
          'str_optimal', 'dt_slct', 'str_traj_rm_jump', 'str_BLH_factor', 'str_remark', 'str_src_z', 'str_domain', 'str_domain_trajfile', 'str_sharpcut', 'NumWorkers', ...
          'basin_name', 'basin_catalog', 'compute_anombyAM', 'forcing'};  
dflts  = cell(1, length(pnames));

[          project_name,   WSV_name,   WSV_dir,   year_list,   mth,    st_month,  st_day,   ed_month,   ed_day,   str_RHc_dqc,   ldirect,   from_basin,   maxtraj_day,...
           str_optimal,   dt_slct,    str_traj_rm_jump,   str_BLH_factor,  str_remark,   str_src_z,   str_domain,  str_domain_trajfile,  str_sharpcut,  NumWorkers,...
           basin_name,   basin_catalog,  compute_anombyAM,   forcing] ...
               = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments
%%
    %===============================================================================
    %/      Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 15 Feb 2024
    %/
    %/ Description: This function is designd to retrieve WaterSip/WaterDrip products   
    %/
    %===============================================================================
    
    if isempty(forcing)   
        forcing = 'CERA';   %/ By default. CERA-forced FLEXPART simulation is only up to 2010-12-30
    end

    if ~isempty(str_optimal)
        str_optimal_traj = '_optimal';     
    else
        str_optimal_traj = [];  
    end

    %/ Change WSV_dir to char type
    if iscell(WSV_dir) && length(WSV_dir) == 1
        WSV_dir = WSV_dir{:};
        % else
            % error('retrieve_WSV cannot handle ''WSV_dir'' with multiple cells!');
        % end
    end

    % basin_name
    % project_name
    if from_basin
        if isequal(project_name, 'TP') && isequal(basin_name, 'TP')
            ind_basin = 1;
        else
            basin_namelist = [basin_catalog.name];
            ind_basin = findismember_loop(basin_namelist, basin_name);
    
            if isempty(ind_basin)  
                error('Detected that ''ind_basin'' is empty when from_basin is not 0! Check your code!');
            end
        end
    else
        ind_basin = [];
    end
    
    if isempty(str_remark)       
        str_remark = '_';    %/ then replace it with '_' to avoid reading a wrong file.
    end   

    %/ Set dates
    mth_bc            = mth;
    st_month_bc       = st_month;
    st_day_bc         = st_day;
    ed_month_bc       = ed_month;
    ed_day_bc         = ed_day;
    
    %/ Retireve data based on 'mth' if the period is not given
    if ~(isempty(st_month_bc) && isempty(st_day_bc) && isempty(ed_month_bc) && isempty(ed_day_bc))
        dt_target_dates = date_array_gen('year_list', year_list, 'st_month', st_month_bc, 'st_day', st_day_bc,...
                                   'ed_month', ed_month_bc, 'ed_day', ed_day_bc,...   
                                   'dt_slct_hr', dt_slct, 'output_date_format', 'datetime', 'skip_the_incomplete', 1);
    else
        if mth == 0
            season = [];
        else
            season = mth_bc - 12;
        end
        dt_target_dates = date_array_gen('year_list', year_list, 'season', season,...   
                                   'dt_slct_hr', dt_slct, 'output_date_format', 'datetime', 'skip_the_incomplete', 1);
    end
    disp(dt_target_dates([1:3,end-2:end]))

    %/ Intialization 
    WSV_data_daily     = []; 
    WSV_dates          = [];
    rr                 = [];
    traj_catalog       = []; 
    % traj_catalog.traj  = [];
    % traj_catalog.q     = [];
    % traj_catalog.T     = [];
    % traj_catalog.date  = [];
    % traj_catalog.ntraj = [];
    
    traj_real = cell(length(dt_target_dates),1);
    q         = cell(length(dt_target_dates),1);
    T         = cell(length(dt_target_dates),1);
    date      = cell(length(dt_target_dates),1);
    ntraj     = cell(length(dt_target_dates),1);

    %/ Date loop 
    %/  We can't run parallel because of the way we accumulate subdaily to daily.
    %/  This loop should be more robust and flexible than a year loop
    tic
    if isequal(WSV_name, 'traj')
        if isempty(NumWorkers)
            NumWorkers = 5;  
        end
        parpool_core = 'Processes';
        if isempty(gcp('nocreate'))  %/ if set worker number
            parpool(parpool_core, NumWorkers) %/ CAVEAT: Loading V7.3 MAT-files on threads is not supported.
        end
        tic;
        parfor t = 1:length(dt_target_dates) 
        % for t = 1:length(dt_target_dates) 
            dt_target_dates_t = dt_target_dates(t,:);
            int_target_date_t = datetime2int(dt_target_dates_t, 'yyyymmddHHMM');
            
            traj_suffix = strcat('_', str_RHc_dqc, '_', ldirect,'_',...
                                 num2str(maxtraj_day), 'd', str_sharpcut, str_optimal_traj, '_', num2str(dt_slct),...
                                 'h_', str_domain_trajfile,'_', str_src_z, str_BLH_factor, str_remark,...
                                 datestr(dt_target_dates_t, 'yyyymmddHHMM'), '.mat');


            if iscell(WSV_dir) && length(WSV_dir) > 1
                %/ [IMPORTANT] Since DJF data will use two-year data due to the 'date_array_gen' function. Need to correctly locate the dir.
                b = cellfun(@(x) strsplit(x, {'_', '/'}), {WSV_dir}, 'UniformOutput', false);  %/ cellfun to strplit each cell.
                b = cellfun(@(x)                x{end-1},         b, 'UniformOutput', false);  %/ cellfun to select year component from each cell.
                ind = findismember_loop(b, {dt_target_dates_t(1:4)});
                
                %/ Load the target data
                traj_fullpath = strcat(WSV_dir{ind}, WSV_name, traj_suffix);
            else
                traj_fullpath = strcat(WSV_dir, WSV_name, traj_suffix);
            end
            
            if ismember(ldirect, {'bwd'})
                traj_var = 'domfill';

            elseif ismember(ldirect, {'fwd'})
                traj_var = 'fwd_final';
            end

            if isfile(traj_fullpath)
                fprintf('*** [t = %d/%d] Loading traj data: %s... ***\n', t, length(dt_target_dates), traj_fullpath)
                domfill = par_load(traj_fullpath, traj_var);
            else
                error('!!! Traj data not found: %s !!!\nHave you run moisture_tracking.m?', traj_fullpath);
            end

            if ismember(ldirect, {'fwd'})  %/ For fwd traj file, it contains trajs from all basins. Subsetting is thus needed.
                domfill = domfill{ind_basin};
            end
            
            %/ 1. Correct the lon range to strictly in (-180, 180]) (sometimes there are x = 182.xx)
            x = domfill.traj(:,1,:);
            x(x > 180) = x(x > 180) - 360;
            domfill.traj(:,1,:) = x;
            
            %/ Reshape (231,161) to (231,1,161), for example.
            domfill.q    = reshape(domfill.q,    size(domfill.q,1),    1, size(domfill.q,2));
            domfill.T    = reshape(domfill.T,    size(domfill.T,1),    1, size(domfill.T,2));
            domfill.topo = reshape(domfill.topo, size(domfill.topo,1), 1, size(domfill.topo,2));
            
            %/ Restore the true altitude of the particle.
            domfill.traj_real        = domfill.traj;
            domfill.traj_real(:,3,:) = domfill.traj(:,3,:) + domfill.topo;
            sz = size(domfill.traj_real);

            %/ NOTE: For fwd trajs, we may have trajs that were sharp cut
            %/       on a certain date, as described by 'str_sharpcut'
            %/       We need to expand the trajtime dimension to ensure it is fixed. 
            max_trajtime = maxtraj_day*24/dt_slct + 1;
            append_n     = max_trajtime - sz(3);
            domfill.traj_real(:,:,end+1:end+append_n) = nan(sz(1), sz(2), append_n);
            domfill.q(:,:,end+1:end+append_n)         = nan(sz(1),     1, append_n);
            domfill.T(:,:,end+1:end+append_n)         = nan(sz(1),     1, append_n);

            %/ Store in cells
            traj_real{t} = domfill.traj_real;
            q{t}         = domfill.q;
            T{t}         = domfill.T;
            date{t}      = int_target_date_t;
            ntraj{t}     = sz(1);

            %/ Concatenate into trajs_catalog fields
            % disp(size(domfill.traj_real))
            % traj_catalog.traj  = cat(1, traj_catalog.traj,  domfill.traj_real);
            % traj_catalog.q     = cat(1, traj_catalog.q,     domfill.q);
            % traj_catalog.T     = cat(1, traj_catalog.T,     domfill.T);
            % traj_catalog.date  = cat(1, traj_catalog.date,  int_target_date_t);
            % traj_catalog.ntraj = cat(1, traj_catalog.ntraj, sz(1));
        end
        traj_catalog.traj  = cat(1, traj_real{:});
        traj_catalog.q     = cat(1, q{:});
        traj_catalog.T     = cat(1, T{:});
        traj_catalog.date  = cat(1, date{:});
        traj_catalog.ntraj = cat(1, ntraj{:});
        % disp(traj_catalog)
        toc;
    else
        %/ Compute the cumulated fwd contr from a region to elsewhere on a certain date
        if ismember(ldirect, {'fwd'})
            %/ Set indices of the end of the days -> to calc daily uptake value.
            ind_EOD = [find(diff(datetime2int(dt_target_dates, 'yyyymmdd')) ~= 0); length(dt_target_dates)]; 
            nday    = length(ind_EOD);
            WSV_data = 0; freq_map = 0; 

            cnt = 1;
            for t = 1:length(dt_target_dates) 
                dt_target_dates_t = dt_target_dates(t,:);
                int_target_date_t = datetime2int(dt_target_dates_t, 'yyyymmddHHMM');
        
                if isempty(NumWorkers)
                    NumWorkers = 60;  
                end
                parpool_core = 'Processes';
                if isempty(gcp('nocreate'))  %/ if set worker number
                    parpool(parpool_core, NumWorkers)   %/ CAVEAT: Loading V7.3 MAT-files on threads is not supported.
                end
    
                %/ A second loop to sum over the fwd contr
                maxtraj_hr = maxtraj_day*24 - dt_slct; %/ Maximum traj time steps. Why -dt_slct? Because we need at least 3 time steps to compute forward contr.
                ntt        = maxtraj_hr/dt_slct;
                parfor tt = 1:ntt
                    dt_src_date  = dt_target_dates_t - hours(tt*dt_slct);
                    str_target_date = datestr(dt_src_date, 'yyyymmddHHMM');
    
                    %/ Construct the suffix to read files
                    str_sharpcut_bc = '';
                    if ~isempty(str_sharpcut) 
                        date_sharpcut    = strrep(str_sharpcut,'_sharpcut','');
                        date_sharpcut_dt = datetime(date_sharpcut,"InputFormat",'yyyyMMddHHmm'); %/ String to datetime
                        date_begins_sharpcut = date_sharpcut_dt - days(maxtraj_day); 

                        if dt_src_date > date_begins_sharpcut
                            str_sharpcut_bc = str_sharpcut;
                        end
                    end

                    suffix = strcat('_', str_RHc_dqc, '_', ldirect,'_',...
                                                   num2str(maxtraj_day), 'd', str_sharpcut_bc, str_optimal, '_', num2str(dt_slct), 'h_', str_domain,... 
                                                   '_', str_traj_rm_jump, str_src_z, str_BLH_factor, str_remark, str_target_date, '.mat');

                    basin_suffix = strcat('_', str_RHc_dqc, '_', ldirect,'_',...
                                                   num2str(maxtraj_day), 'd', str_sharpcut_bc, str_optimal, '_', num2str(dt_slct), 'h_', basin_name,... 
                                                   '_', str_traj_rm_jump, str_src_z, str_BLH_factor, str_remark, str_target_date, '.mat');

                    %/ First, load the dates corresponding to forward trajtime 
                    if iscell(WSV_dir) && length(WSV_dir) > 1
                        WSV_dir_bc = WSV_dir{t};
                    else
                        WSV_dir_bc = WSV_dir;
                    end
                    fullpath = strcat(WSV_dir_bc, 'Cf_map_dates', suffix);

                    if isfile(fullpath)         
                        % fprintf('*** [t = %d, tt = %d] Loading %s *** \n', t, tt, fullpath)
                        % dates = load(fullpath);
                        % Cf_map_dates = dates.Cf_map_dates;
                        Cf_map_dates = par_load(fullpath, 'Cf_map_dates');
                        
                    else
                        error('!!! Cf_map_dates NOT found: %s !!!',fullpath);
                    end
    
                    %/ Second, set the index of the target date
                    int_Cf_map_dates = datetime2int(Cf_map_dates, 'yyyymmddHHMM');
                    ind_date = find(ismember(int_Cf_map_dates, int_target_date_t));
                    if isempty(ind_date)
                        error('ind_date is empty!');
                    end
    
                    %/ Third, load the target data
                    fullpath = strcat(WSV_dir_bc, WSV_name, basin_suffix);
                    if isfile(fullpath)         
                        fprintf('*** [t = %d, tt = %d] Loading %s *** \n', t, tt, fullpath)
                        L = par_load(fullpath, WSV_name);
                    else
                        error('!!! %s NOT found: %s !!!', WSV_name, fullpath);
                    end
                    %/ Subset the target date and do summation by reduction
                    L = L(:,:,ind_date);
    
                    a = L;  a(~isnan(a)) = 1;  a(isnan(a))  = 0;
                    freq_map = freq_map + a;
    
                    nonnan_data = L; nonnan_data(isnan(nonnan_data)) = 0; 
                    WSV_data = WSV_data + nonnan_data;            %/ summing if it is uptake (g/kg per 3h) or Pm (mm per 3h)  
                end
    
                if ~isempty(find(isnan(WSV_data), 1))  
                    error('reduction variable contains nan. modify the code!');   
                end
    
                if t == 1                                                      %/ initialize
                    WSV_data_daily = zeros(size(WSV_data,1), size(WSV_data,2), nday);
                end
    
                if any(ismember(t, ind_EOD))                                   %/ store into daily data at the end of day
                    %/ Only the below WSV will have to take average.
                    if ismember(WSV_name, {'optimal_trajtime', 'CWRT', 'RH2', 'rr_tot_L', 'rr_tot_NLL', 'rr_tot_NLO'})
                        WSV_data = WSV_data./freq_map;                         %/ divide by # of non-nans
                    end
                    WSV_data_daily(:,:,cnt) = WSV_data;
                    freq_map = 0;   WSV_data = 0;                              %/ reset
                    cnt = cnt + 1;                                             %/ since we will skip the last time step (or should we?), cnt = 356 or 366
                end
            end
        end
    
        %/ Compute the cumulated bwd contr from a region
        if ismember(ldirect, {'bwd'})
            if isempty(NumWorkers)
                NumWorkers = 40;  
            end

            parpool_core = 'Processes';  %/ NOTE: threads cannot load V7.3 MAT-files
            if isempty(gcp('nocreate'))  
                parpool(parpool_core, NumWorkers) 
            end 

            WSV_data_daily_EachYr = cell(length(year_list), 1);
            rr_EachYr             = cell(length(year_list), 1);

            for y = 1:length(year_list)
            % parfor y = 1:length(year_list)

                %/ Use 'date_array_gen' to allow reading seasons across two years (e.g., DJF)
                if ~(isempty(st_month_bc) && isempty(st_day_bc) && isempty(ed_month_bc) && isempty(ed_day_bc))
                    dt_target_dates_eachYr = date_array_gen('year_list', year_list(y), 'st_month', st_month_bc, 'st_day', st_day_bc,...
                                               'ed_month', ed_month_bc, 'ed_day', ed_day_bc,...   
                                               'dt_slct_hr', dt_slct, 'output_date_format', 'datetime', 'skip_the_incomplete', 0);
                else
                    if mth == 0
                        season = [];
                    else
                        season = mth_bc - 12;
                    end

                    if y == length(year_list) && isequal(season, 4)
                        continue;
                    end

                    dt_target_dates_eachYr = date_array_gen('year_list', year_list(y), 'season', season,...   
                                               'dt_slct_hr', dt_slct, 'output_date_format', 'datetime', 'skip_the_incomplete', 0);
                end  
                ind_EOD_eachYr = [find(diff(datetime2int(dt_target_dates_eachYr, 'yyyymmdd')) ~= 0); length(dt_target_dates_eachYr)]; 
                nday_eachYr    = length(ind_EOD_eachYr);

                %/ Initialization
                cnt = 1; WSV_data = 0; freq_map = 0; WSV_data_daily = 0;
                rr = [];
                rr_pfor = cell(length(dt_target_dates_eachYr), 1); 
                
                for t = 1:length(dt_target_dates_eachYr) 
                    % disp(t)
                    dt_target_dates_t = dt_target_dates_eachYr(t,:);
                    % int_target_date_t = datetime2int(dt_target_dates_t, 'yyyymmddHHMM');
            
                    if isequal(forcing, 'EA')
                        str_target_date = datestr(dt_target_dates_t,   'yyyymmddHHMM');
          
                    elseif isequal(forcing, 'CERA')
                        if isequal(datestr(dt_target_dates_t, 'yyyymmdd'), '20101231')   %/ No Dec 31 in year 2010. Skip it.
                            if any(ismember(t, ind_EOD_eachYr))
                                cnt = cnt + 1;
                            end
                            continue;

                        elseif isequal(datestr(dt_target_dates_t, 'mmddHHMM'), '12312100')
                            str_target_date = datestr(dt_target_dates_eachYr(t-1,:), 'yyyymmddHHMM'); %/ No 1231 21:00 data from flexpart expmnt, so we copy the 1231 18:00 instead.

                        elseif isequal(datestr(dt_target_dates_t, 'yyyymmddHHMM'), '201012302100') 
                            str_target_date = datestr(dt_target_dates_eachYr(t-1,:), 'yyyymmddHHMM'); %/ No Dec 30, 21:00 in year 2010, here we borrow 12301800 to compute daily value.

                        else
                            str_target_date = datestr(dt_target_dates_t,   'yyyymmddHHMM');
                        end
                    else
                        error('str_target_date not set for the forcing ''%s''!', forcing);
                    end
            
                    suffix = strcat('_', str_RHc_dqc, '_', ldirect,'_',...
                                    num2str(maxtraj_day), 'd', str_sharpcut, str_optimal, '_', num2str(dt_slct), 'h_', str_domain,... 
                                    '_', str_traj_rm_jump, str_src_z, str_BLH_factor, str_remark, str_target_date, '.mat');
            
                    %/ Load the target data
                    if iscell(WSV_dir) && length(WSV_dir) > 1
                        WSV_dir_bc = WSV_dir{y};
                    else
                        WSV_dir_bc = WSV_dir;
                    end
                    fullpath = strcat(WSV_dir_bc, WSV_name, suffix);

                    if isfile(fullpath)         
                        fprintf('*** [t = %d] Loading %s *** \n', t, fullpath)
                        L = load(fullpath);
                    else
                        error('!!! File NOT found: %s !!!',fullpath);
                    end
                    fld = fieldnames(L);                                               %/ extract all the fieldname
        
                    if isequal(WSV_name, 'watersip')                                   %/ only for bwd
                        if from_basin
                            rr_pfor{t} = mean(table2array(L.(fld{:}){ind_basin}), 1, 'omitnan'); %/ take the mean recycling ratios of all trajs from the hs at each time step.
                        else
                            if isa(L.(fld{:}), 'double')
                                rr_pfor{t} = mean(L.(fld{:}), 1, 'omitnan');
                            elseif isa(L.(fld{:}), 'table')
                                rr_pfor{t} = mean(table2array(L.(fld{:})), 1, 'omitnan');
                            else
                                error('What type is it??')
                            end
                        end
                    else
                        if from_basin
                            a = L.(fld{:})(:,:,ind_basin);  a(~isnan(a)) = 1;  a(isnan(a))  = 0;
                            freq_map = freq_map + a;
                            nonnan_data = L.(fld{:})(:,:,ind_basin);
                            nonnan_data(isnan(nonnan_data)) = 0; 
                            WSV_data = WSV_data + nonnan_data;                     %/ summing if it is uptake (g/kg per 3h) or Pm (mm per 3h)
                        else
                            a = L.(fld{:});  a(~isnan(a)) = 1;  a(isnan(a))  = 0;
                            freq_map = freq_map + a;
                            nonnan_data = L.(fld{:});
                            nonnan_data(isnan(nonnan_data)) = 0; 
                            WSV_data = WSV_data + nonnan_data;                     %/ summing if it is uptake (g/kg per 3h) or Pm (mm per 3h)
                        end
                    
                        if ~isempty(find(isnan(WSV_data), 1))  
                            error('reduction variable contains nan. modify the code!');   
                        end
                        if t == 1
                            WSV_data_daily = zeros(size(WSV_data,1), size(WSV_data,2), nday_eachYr);
                        end
                        if any(ismember(t, ind_EOD_eachYr))                            %/ store into daily data at the end of day
                            %/ Only the below WSVs will have to take average.
                            if ismember(WSV_name, {'optimal_trajtime', 'CWRT', 'RH2', 'rr_tot_L', 'rr_tot_NLL', 'rr_tot_NLO'})
                                WSV_data = WSV_data./freq_map;                         %/ divide by # of non-nans
                            end
                            WSV_data_daily(:,:,cnt) = WSV_data;
                            freq_map = 0;   WSV_data = 0;                              %/ reset
                            cnt = cnt + 1;                                             %/ since we will skip the last time step (or should we?), cnt = 356 or 366
                        end
                    end
                end
                
                if isequal(WSV_name, 'watersip')                                       %/ need not to take daily data when retrieving rr.
                    rr             = cat(1, rr_pfor{:});                               %/ since there are cell arrays in cell arrays
                    rr             = cat(1, rr{:});
                else
                    if cnt ~= nday_eachYr + 1
                        disp(cnt)
                        disp(nday_eachYr + 1)
                        error('check if the quired dates have data!!')
                    end
                end
    
                %/ Store results from the date loop into a cell array
                WSV_data_daily_EachYr{y} = WSV_data_daily;
                rr_EachYr{y}             = rr;
            end

            %/ Concatenate the cells
            WSV_data_daily = cat(3, WSV_data_daily_EachYr{:});
            rr             = cat(1, rr_EachYr{:});
        end
        toc
    
        % A              = cellfun(@(x) x.', slct_dates_dt, 'UniformOutput',false);  %/ convert to row vector in each cell
        % WSV_dates      = cat(2, unique(datetime2int([A{:}], 'yyyymmdd')))';
        WSV_dates = unique(datetime2int(dt_target_dates, 'yyyymmdd'));
    
        if compute_anombyAM
            %/ Incurrsion
            [WSV_data_daily_all, ~, ~] = retrieve_WSV('WSV_name', WSV_name, 'WSV_dir', WSV_dir, 'year_list', year_list, 'mth', 0, 'str_RHc_dqc', str_RHc_dqc, 'ldirect', ldirect, 'from_basin', from_basin, 'maxtraj_day', maxtraj_day,...
                                                         'str_optimal', str_optimal, 'dt_slct', dt_slct, 'str_traj_rm_jump', str_traj_rm_jump, 'str_BLH_factor', str_BLH_factor, 'str_remark', str_remark, 'str_src_z', str_src_z,...
                                                         'str_domain', str_domain, 'str_sharpcut', str_sharpcut, 'NumWorkers', NumWorkers,...
                                                         'basin_name', basin_name, 'basin_catalog', basin_catalog, 'forcing', forcing);
    
            WSV_data_AM = mean(WSV_data_daily_all, 3, 'omitnan');          %/ AM = annual mean
            WSV_data_daily = WSV_data_daily - WSV_data_AM;                 %/ anomalies w.r.t. annual mean
        end
    
        %/ Quick check
        fprintf('*** Max(WSV_data_daily) = %.2f ***\n', max(WSV_data_daily, [], 'all'));
        fprintf('*** Mean(WSV_data_daily) = %.2f ***\n', mean(WSV_data_daily, 'all'));
        fprintf('*** Min(WSV_data_daily) = %.2f ***\n', min(WSV_data_daily, [], 'all'));
    end

    %/ Shut down parpool (only for Threads; avoid shutting down for
    %/ Processes, which would spend more time in restarting.
    if isequal(parpool_core, 'Threads')
        poolobj = gcp('nocreate');
        delete(poolobj);
    end
end