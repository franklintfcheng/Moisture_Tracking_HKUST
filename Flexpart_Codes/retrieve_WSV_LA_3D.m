function [WSV_data_daily_zm, WSV_data_daily_mm, WSV_dates,...
          basin_lon_range, basin_lat_range, basin_z_range, zm_range, mm_range] = retrieve_WSV_LA_3D(varargin)

    %/ create a set of valid parameters and their default value
    pnames = {'WSV_name', 'WSV_dir', 'year_list', 'mth', 'str_RHc_dqc', 'ldirect', 'from_basin', 'maxtraj_day',...
              'str_optimal', 'dt_slct', 'str_traj_rm_jump', 'str_BLH_factor', 'str_remark', 'str_src_z', 'str_domain', 'str_sharpcut', 'NumWorkers',...
              'slct_reg', 'basin_name', 'basin_bndry', 'basin_catalog'};  
    dflts  = cell(1, length(pnames));

    [          WSV_name, WSV_dir, year_list, mth, str_RHc_dqc, ldirect, from_basin, maxtraj_day,...
               str_optimal, dt_slct, str_traj_rm_jump, str_BLH_factor, str_remark, str_src_z, str_domain, str_sharpcut, NumWorkers,...
               slct_reg,   basin_name,   basin_bndry,   basin_catalog] ...
                   = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments

%%
    if (from_basin == 0 && isequal(ldirect, 'bwd')) || isequal(ldirect, 'fwd')
        error('This function is not yet designed for from_basin = %d and ldirect = %s', from_basin, ldirect);
    end

    WSV_data_daily_zm  = cell(length(year_list), 1);
    WSV_data_daily_mm  = cell(length(year_list), 1);
    slct_dates_dt      = cell(length(year_list), 1);
    
    if ~isempty(slct_reg)
        slct_bndry = box_region(slct_reg);
    else
        slct_bndry = [];
    end
    
    %/ NOTE: 'basin_catalog' can be useful if later wanted to compute freq_3D 
    %/       for basin_name = 'TP' from the 15 TP basins.
    if isequal(basin_name, 'TP')
        basin_namelist  = [basin_catalog.name];
    else
        basin_namelist  = {basin_name};
    end
    
    %/============ Retrieve basic info ============%      
    a = date_array_gen('year_list', year_list(1), 'st_month', 1, 'st_day', 1,...
                       'ed_month', 12, 'ed_day', 31, 'dt_slct_hr', dt_slct, 'output_date_format', 'datetime');

    str_target_date = datestr(a(1), 'yyyymmddHHMM');

    suffix = strcat('_', str_RHc_dqc, '_', ldirect,'_',...
                    num2str(maxtraj_day), 'd', str_sharpcut, str_optimal, '_', num2str(dt_slct), 'h_', str_domain,... 
                    '_', str_traj_rm_jump, str_src_z, str_BLH_factor, str_remark, str_target_date, '.mat');

    % suffix = strcat('_', str_RHc_dqc, '_', ldirect,'_',...
    %                                num2str(maxtraj_day), 'd', str_optimal, '_', num2str(dt_slct), 'h_', str_domain,... 
    %                              '_', str_traj_rm_jump, str_src_z, str_BLH_factor, str_remark, str_target_date, '.mat');

    %/ Since DJF data will use two-year data due to the 'date_array_gen' function. Need to correctly locate the dir.
    b = cellfun(@(x) strsplit(x, {'_', '/'}), WSV_dir, 'UniformOutput', false);  %/ cellfun to strplit each cell.
    b = cellfun(@(x)                x{end-1},       b, 'UniformOutput', false);  %/ cellfun to select year component from each cell.
    ind = findismember_loop(b, {str_target_date(1:4)});

    %/ Since temp vas in the parfor loop cannot be saved. 
    %/ we read an arbitrary file to get the basic info about (lon, lat, z).
    fullpath = strcat(WSV_dir{ind}, 'LA_3D_', basin_namelist{1}, suffix);
    load(fullpath, 'LA_3D');  %/ output is LA_3D

    %/ NOTE: I set the z range the same for all basins, but not for lon and
    %/       lat range.. Hence, we need a mother domain to integrate the data.
    res        = 1;  %/ Here we assume the horizontal resolution is 1 deg. Adapt it if not the case!
    mother_lon = 0:res:359;  
    mother_lat = -90:res:90;
    mother_z   = LA_3D.basin_z_range;  
    
    %/ Basin Loop (if basin_name == 'TP')
    %/ Year Loop
    if isempty(gcp('nocreate')) && ~isempty(NumWorkers) %/ if set worker number
        % parpool('Threads', NumWorkers) %/ Loading V7.3 MAT-files on threads is not supported.
        parpool('Processes', NumWorkers) %/ use process-based parpool
    end
    parfor y = 1:length(year_list)
        season = []; st_month = []; ed_month = []; st_day = []; ed_day = [];

        %/ Broadcast vars
        year_list_bc    = year_list;
        basin_bndry_bc  = basin_bndry; 
        slct_bndry_bc   = slct_bndry;

%         if test_mode        %/ only consider Jan 1.
%             fprintf('==============================================================================\n')
%             fprintf('============= Test Mode On (Only run Fed 1, 8 time steps) ===============\n')
%             fprintf('==============================================================================\n')
%             st_month = 2; st_day = 1;   
%             ed_month = 2; ed_day = 10;
%         else
            if  mth == 0    %/ mth=0 --> whole period
                st_month = 1;     st_day = 1;   ed_month  = 12; 
                if year_list_bc(y) == 2010 && ismember(ldirect, {'bwd'})  
                    ed_day = 30;
                elseif year_list_bc(y) == 2010 && ismember(ldirect, {'fwd'})  
                    ed_day = 15;  
                else
                    ed_day = 31;
                end

            elseif mth >= 13
                season = mth - 12;    %/ 13: MAM, 14: JJA, 15: SON, 16: DJF, 17: Apr-Sep, 18: Oct-Mar
            else
                st_month  = mth; st_day = 1;
                ed_month  = st_month; ed_day = eomday(year_list_bc(y),ed_month);
            end
%         end
        if y == length(year_list_bc) && isequal(season, 4)    continue;   end           %/ skip it since no complete DJF for the last year. Do NOT write season == 4, since season can be [].

        slct_dates_dt{y} = date_array_gen('year_list', year_list_bc(y), 'st_month', st_month, 'st_day', st_day,...
                                          'ed_month', ed_month, 'ed_day', ed_day, 'season', season,...   %/ omit 20101231 if using season function.
                                          'dt_slct_hr', dt_slct, 'output_date_format', 'datetime');

        ind_EOD = [find(diff(yyyymmdd(slct_dates_dt{y})) ~= 0); length(slct_dates_dt{y})]; %/ index of the end of the day -> to compute daily value.
        nday = length(unique(yyyymmdd(slct_dates_dt{y})));

        if ismember(WSV_name, {'BLH_2D'})
            WSV_data     = zeros(length(mother_lon), length(mother_lat)); 
            No_Of_NonNaN = zeros(length(mother_lon), length(mother_lat)); 
        else
            WSV_data     = zeros(length(mother_lon), length(mother_lat), length(mother_z)); 
            No_Of_NonNaN = zeros(length(mother_lon), length(mother_lat), length(mother_z)); 
        end
        cnt = 1; 
        for t = 1:length(slct_dates_dt{y})  %/ this loop can NOT run parallel, since we need to accum subdaily to daily.
            if isequal(datestr(slct_dates_dt{y}(t,:), 'mmddHHMM'), '12312100')
                %/ NOTE: No 1231 21:00 data from flexpart expmnt, so we copy the 1231 18:00 instead.
                str_target_date = datestr(slct_dates_dt{y}(t-1,:), 'yyyymmddHHMM'); 
            elseif year_list_bc(y) == 2010 && isequal(datestr(slct_dates_dt{y}(t,:), 'mmddHHMM'), '12302100') %/ only till 20101230
                str_target_date = datestr(slct_dates_dt{y}(t-1,:), 'yyyymmddHHMM'); 
            else
                str_target_date = datestr(slct_dates_dt{y}(t,:),   'yyyymmddHHMM');
            end

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

            % suffix = strcat('_', str_RHc_dqc, '_', ldirect,'_',...
            %                      num2str(maxtraj_day), 'd', str_optimal, '_', num2str(dt_slct), 'h_', str_domain,... 
            %                      str_traj_rm_jump, str_src_z, str_BLH_factor, str_remark, str_target_date, '.mat');

            %/ Since DJF data will use two-year data due to the 'date_array_gen' function. Need to correctly locate the dir.
            b = cellfun(@(x) strsplit(x, {'_', '/'}), WSV_dir, 'UniformOutput', false);  %/ cellfun to strplit each cell.
            b = cellfun(@(x)                x{end-1},       b, 'UniformOutput', false);  %/ cellfun to select year component from each cell.
            ind = findismember_loop(b, {str_target_date(1:4)});

            %/ Loading LA_3D file for the basin(s)
            LA_3D_all = []; 
            if ismember(WSV_name, {'BLH_2D'})
                LA_3D_all.(WSV_name) = zeros(length(mother_lon), length(mother_lat)); 
                counting             = zeros(length(mother_lon), length(mother_lat));
            else
                LA_3D_all.(WSV_name) = zeros(length(mother_lon), length(mother_lat), length(mother_z)); 
                counting             = zeros(length(mother_lon), length(mother_lat), length(mother_z));
            end
            
            for top = 1:length(basin_namelist)
                basin_name_each = basin_namelist{top};

                fullpath = strcat(WSV_dir{ind}, 'LA_3D_', basin_name_each, suffix);
                fprintf('*** [t = %d] Loading from it... %s *** \n', t, fullpath)
                LA_3D            = par_load(fullpath, 'LA_3D');   
                fld              = fieldnames(LA_3D);
                fld_to_remove    = setdiff(fld, [WSV_name, {'basin_lon_range', 'basin_lat_range', 'basin_z_range'}]);
%                 fld_to_remove    = setdiff(fld, [WSV_name, {'freq_3D', 'basin_lon_range', 'basin_lat_range', 'basin_z_range'}]);
                LA_3D            = rmfield(LA_3D, fld_to_remove);
                
                if abs(unique(diff(LA_3D.basin_lon_range))) ~= res
                    disp(res)
                    disp(abs(unique(diff(LA_3D.basin_lon_range))))
                    error('Please modify ''res'' in the function to fit the actual res of the data!');
                end
                ind_lon = findismember_loop(mother_lon, LA_3D.basin_lon_range);
                ind_lat = findismember_loop(mother_lat, LA_3D.basin_lat_range);
                ind_z   = findismember_loop(mother_z,   LA_3D.basin_z_range);
                
%                 size(LA_3D.(WSV_name))
                %/ Always set nan to 0 prior to reduction
                nonnan_data = LA_3D.(WSV_name);
                nonnan_data(isnan(nonnan_data)) = 0;  
%                 size(nonnan_data)
                
                if ismember(WSV_name, {'BLH_2D'})
                    LA_3D_all.(WSV_name)(ind_lon, ind_lat)        = LA_3D_all.(WSV_name)(ind_lon, ind_lat)        + nonnan_data;           %/ Do summation first
                    counting(ind_lon, ind_lat)                    = counting(ind_lon, ind_lat)                    + logical(nonnan_data);  %/ Do not use freq_3D, as ins vars like u_3D have already been computed based on freq_3D.
                else
                    LA_3D_all.(WSV_name)(ind_lon, ind_lat, ind_z) = LA_3D_all.(WSV_name)(ind_lon, ind_lat, ind_z) + nonnan_data;           %/ Do summation first
                    counting(ind_lon, ind_lat, ind_z)             = counting(ind_lon, ind_lat, ind_z)             + logical(nonnan_data);  %/ Do not use freq_3D, as ins vars like u_3D have already been computed based on freq_3D.
                end
%                 size(LA_3D_all.(WSV_name))
            end

            if ismember(WSV_name, {'BLH_2D', 'rr_3D', 'u_3D', 'v_3D', 'w_3D'})
                LA_3D_all.(WSV_name) = LA_3D_all.(WSV_name)./counting; %/ get the mean u_3D, etc. from the basin(s)
            end
            
            %/ Hence, by default, we strictly subset the basin boundary
            %/   range for zonal/meridional averaging
            if ~isempty(slct_bndry_bc)
                zm_range = mother_lon(mother_lon >= min(slct_bndry_bc(:,1)) & mother_lon <= max(slct_bndry_bc(:,1)));  
                mm_range = mother_lat(mother_lat >= min(slct_bndry_bc(:,2)) & mother_lat <= max(slct_bndry_bc(:,2))); 
            else
                zm_range = mother_lon(mother_lon >= min(basin_bndry_bc(:,1)) & mother_lon <= max(basin_bndry_bc(:,1)));  
                mm_range = mother_lat(mother_lat >= min(basin_bndry_bc(:,2)) & mother_lat <= max(basin_bndry_bc(:,2))); 
            end

            %/ Always set nan to 0 prior to reduction
            nonnan_data = LA_3D_all.(WSV_name);
            nonnan_data(isnan(nonnan_data)) = 0;

            if numel(nonnan_data) == 1 && nonnan_data == 0
                warning('Data is missing (due likely to no traj in this basin at this time step.')
                %/ To avoid bug, create a zero 2D or 3D data.
                if ismember(WSV_name, {'BLH_2D'})
                    nonnan_data = zeros(length(mother_lon), length(mother_lat));
                else
                    nonnan_data = zeros(length(mother_lon), length(mother_lat), length(mother_z));
                end
            end
%             size(nonnan_data)
            No_Of_NonNaN = No_Of_NonNaN + logical(nonnan_data); %/ using logical to count the # of Non-NaN in a day at each grid.
            WSV_data     = WSV_data + nonnan_data;                                 

            if ~isempty(find(isnan(WSV_data), 1))  
                error('reduction variable contains nan. modify the code!');   
            end

            if t == 1  %/ initialize
                WSV_data_daily_zm{y} = squeeze(zeros(size(WSV_data,2), size(WSV_data,3), nday));
                WSV_data_daily_mm{y} = squeeze(zeros(size(WSV_data,1), size(WSV_data,3), nday));
            end

            %/ store into daily data at the end of day
            if any(ismember(t, ind_EOD))                
                %/ Perform zonal mean (zm) / meridional mean (mm) on
                %/ each day to save memory space.
                if ismember(WSV_name, {'BLH_2D'})
                    WSV_data = WSV_data./No_Of_NonNaN;

                    %/ NOTE: Due to the original design of 'compute_zm_mm',
                    %/       Here when zm_or_mm == 1, we specify the zm_range 
                    %/       but in the largest lat range availble. 
                    %/       Similarly for zm_or_mm == 2.
                    mean_or_sum = 1;
                    [WSV_data_daily_zm{y}(:,cnt), ~] = compute_zm_mm('geo_data', WSV_data, 'lon',  mother_lon, 'lat', mother_lat,...
                                  'lon_range', zm_range, 'lat_range', mother_lat, 'zm_or_mm', 1, 'mean_or_sum', mean_or_sum);

                    [WSV_data_daily_mm{y}(:,cnt), ~] = compute_zm_mm('geo_data', WSV_data, 'lon',  mother_lon, 'lat', mother_lat,...
                                  'lon_range', mother_lon, 'lat_range', mm_range, 'zm_or_mm', 2, 'mean_or_sum', mean_or_sum);

                else
                    %/ NOTE: Due to the original design of 'compute_zm_mm',
                    %/       Here when zm_or_mm == 1, we specify the zm_range 
                    %/       but in the largest lat range availble. 
                    %/       Similarly for zm_or_mm == 2.
                    if ismember(WSV_name, {'freq_3D'})
                        %/ DO nothing on WSV_data -> get the daily total freq_3D.
                        mean_or_sum = 2;   %/ Only for freq_3D we take zonal / meridional sum.
                    else
                        WSV_data = WSV_data./No_Of_NonNaN;
                        mean_or_sum = 1;
                    end
                    [WSV_data_daily_zm{y}(:,:,cnt), ~] = compute_zm_mm('geo_data', WSV_data, 'lon',  mother_lon, 'lat', mother_lat,...
                                  'lon_range', zm_range, 'lat_range', mother_lat, 'zm_or_mm', 1, 'mean_or_sum', mean_or_sum);

                    [WSV_data_daily_mm{y}(:,:,cnt), ~] = compute_zm_mm('geo_data', WSV_data, 'lon',  mother_lon, 'lat', mother_lat,...
                                  'lon_range', mother_lon, 'lat_range', mm_range, 'zm_or_mm', 2, 'mean_or_sum', mean_or_sum);
                end

    %             freq_map = 0;   
                WSV_data = 0;      %/ reset
                No_Of_NonNaN = 0;  %/ reset
                cnt = cnt + 1;       %/ since we will skip the last time step (or should we?), cnt = 356 or 366
            end
        end
    end

    if ismember(WSV_name, {'BLH_2D'})
        WSV_data_daily_zm = cat(2, WSV_data_daily_zm{:});              %/ 2nd dim = dates
        WSV_data_daily_mm = cat(2, WSV_data_daily_mm{:});              %/ 2nd dim = dates
    else
        WSV_data_daily_zm = cat(3, WSV_data_daily_zm{:});              %/ 3rd dim = dates
        WSV_data_daily_mm = cat(3, WSV_data_daily_mm{:});              %/ 3rd dim = dates
    end

    A = cellfun(@(x) x.', slct_dates_dt, 'UniformOutput',false);       %/ convert to row vector in each cell
    WSV_dates = cat(2, unique(yyyymmdd([A{:}])))';

    %/ update
    basin_lon_range = mother_lon;
    basin_lat_range = mother_lat;
    basin_z_range   = mother_z/1e3; %/ at last, convert from m to km.
    
    %/ Output the range for zonal/meridional averaging
    if ~isempty(slct_bndry)
        zm_range = mother_lon(mother_lon >= min(slct_bndry(:,1)) & mother_lon <= max(slct_bndry(:,1)));  
        mm_range = mother_lat(mother_lat >= min(slct_bndry(:,2)) & mother_lat <= max(slct_bndry(:,2))); 
    else
        zm_range = mother_lon(mother_lon >= min(basin_bndry(:,1)) & mother_lon <= max(basin_bndry(:,1)));  
        mm_range = mother_lat(mother_lat >= min(basin_bndry(:,2)) & mother_lat <= max(basin_bndry(:,2)));
    end
    
end