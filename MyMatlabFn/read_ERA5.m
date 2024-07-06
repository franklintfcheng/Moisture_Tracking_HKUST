function [output, fld_data, fld_lon, fld_lat, fld_date_yyyymmdd_AllYr] = read_ERA5(varargin)

    %/ Create a set of valid parameters and their default value
    pnames = {'datafolder',         'filename_prefix',  'filename_part', 'dataname',...
              'dataname_callfile',  'varname',          'dataunitconv',  'datatype',     'dataunit',...
              'stlon',              'edlon',            'stlat',         'edlat',...
              'stlevel',            'edlevel',          'slct_level',    'case_date', 'StepsInADay',...
              'timeshift',          'select_field',     'NumWorkers',    'NonStruct'};

    dflts  = cell(length(pnames), 1);

    %/ Parse function arguments
             [ datafolder,          filename_prefix,   filename_part,   dataname,...
               dataname_callfile,   varname,           dataunitconv,    datatype,         dataunit,...
               stlon,               edlon,             stlat,           edlat,...
               stlevel,             edlevel,           slct_level,      case_date,     StepsInADay,...
               timeshift,           select_field,      NumWorkers,      ~ ] ...
                   = internal.stats.parseArgs(pnames, dflts, varargin{:});
%%
    %================================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: 7 Feb 2024
    %/
    %/ This is a function to read *all* or *partial* continuous time entries in ERA5 3D/4D data /% 
    %/
    %/ Available output fields: 
    %/                  lon, lat, time, daily
    %/                  date_yyyymmdd_AllYr, date_mmdd_AllYr, date_mmdd_OneYr
    %/                  daily_clim, PMC, dailyAnom_PMC
    %/                  monthly, monthlyAnom
    %/
    %/          mode = 1: Normal data extraction. Compute Pentad-day moving mean c
    %/                    limatology based on given time period.
    %/
    %/          mode = 2: Similar to mode 1, but when the given time period
    %/                    is short or/and memory in server is limited (Deprecated)
    %================================================================================
    %%
    % tic
    output = [];
    if numel(num2str(case_date(1))) == 8
        year = unique(round(case_date/1e4));
    else
        error('incorrect format of case_date!');
    end
    
    if ~(isequal(datatype, 'acc') || isequal(datatype, 'ins'))
        disp(datatype)
        error('Specify ''datatype'' properly! ''ins'' (instantaneous) or ''acc'' (accumulated)?');
    end

    % if mode == 2 && isempty(allyears)   error('Specify ''allyears''!');              end %/ allyears is only used in mode2. Depricated.
    % disp(allyears)
    % if isempty(NumWorkers) NumWorkers = Inf; end  %/ if NumWorkers = 0 --> no parallel computing.

    %/ Broadcast var
    dataname_callfile_bc = dataname_callfile;
    varname_bc           = varname;
     
    if isempty(stlon)     stlon = 0;        end 
    if isempty(edlon)     edlon = 360;      end
    if isempty(stlat)     stlat = -90;      end
    if isempty(edlat)     edlat = 90;       end

    %/ set timeshift = 0 if not specified or datatype = 'ins'.
    if isempty(timeshift)
        timeshift = 0;
    else
        if isequal(datatype, 'ins')
            timeshift = 0;
        end
    end
    % disp(['NOTE: timeshift = ', num2str(timeshift), ' days.']);
    dn1900 = datenum(1900, 1, 1, 0, 0, 0);  %/ in days

    %/ Read basic info (e.g., lat, lon, level, time) 
    styear = year(1);
    if isempty(filename_prefix)    
        error('filename_prefix is not specified! Note: it reads "filename_prefix + dataname_callfile + filename_part + year(t) + .nc"');
    elseif isempty(filename_part)  
        error('filename_part is not specified!   Note: it reads "filename_prefix + dataname_callfile + filename_part + year(t) + .nc"');
    else                           
        filename_full = strcat(filename_prefix, dataname_callfile_bc, filename_part);        
    end

    %/ Search for complete data. If not existed, search for incomplete data
    datapath_st = [];
    suffix = {'', '_incomplete'}; 
    for i = 1:length(suffix)
        filename_st = strcat(filename_full, num2str(styear),suffix{i},'.nc');
        datapath_st = strcat(datafolder,filename_st);
        if isfile(datapath_st)
            break;
        end
    end
    if isempty(datapath_st)
        error('File %s not exist!')
    end

    ncdisp(datapath_st);
    lon = double(ncread(datapath_st,'longitude'));
    lat = double(ncread(datapath_st,'latitude'));
    lon = lon(stlon <= lon & lon <= edlon);
    lat = lat(stlat <= lat & lat <= edlat);
    nlon = length(lon);
    nlat = length(lat);
    output.lon = lon; %/ update 
    output.lat = lat; %/ update
    units      = ncreadatt(datapath_st, varname_bc, 'units');
    if isequal(datatype, 'acc')
        if contains(select_field{1}, 'daily')
            output.units = strcat(units, {' day**-1'});
        elseif contains(select_field{1}, 'monthly')
            output.units = strcat(units, {' month**-1'});
        else
            error('code not set for ''%s''!', select_field{1});
        end
    else
        if ~isempty(dataunit)
            output.units = dataunit;  %/ If given (useful when dataunitconv is not 1)
        else
            output.units = units;     %/ Will save the original unit for the queried field
        end
    end
    
    %/ Check if the 'level variable exists
    try
        output.level = double(ncread(datapath_st,'level'))';
        if ~isempty(slct_level) %/ Based on slct_level if given, otherwise stlevel and edlevel
            ind_level = find(ismember(output.level, slct_level));
        else
            if stlevel > edlevel
                %/ Switch the two
                ind_level = find(edlevel <= output.level & output.level <= stlevel);
            else
                ind_level = find(stlevel <= output.level & output.level <= edlevel);
            end
        end
        nlevel = length(ind_level);  
        output.level = output.level(ind_level);
        flag_4D = 1; time_dim = 4; %/ Assume
    catch
        ind_level = []; nlevel = []; %/ still need to specify them as broadcast vars
        flag_4D = 0; time_dim = 3; %/ Assume
    end
    
    if flag_4D == 1 && isempty(ind_level)
        error('ind_level is empty! Check your input!');
    end

    subdaily_1D_list        = cell(length(year),1);
    daily_1D_list           = cell(length(year),1);
    yearly_sum_1D_list      = cell(length(year),1);
    time_list               = cell(length(year),1);
    
    if isempty(gcp('nocreate')) && ~isempty(NumWorkers)   
        parpool('Processes', NumWorkers);    
        % parpool('Threads', NumWorkers);   %/ Threads cannot work with ncread!
    end

    filename = []; subdaily = []; daily = [];
    parfor t = 1:length(year)
    % for t = 1:length(year)
    
        %/ KEEP the initialization to avoid warnings for temporary var generated in the parfor loop.
        filename        = [];
        subdaily        = [];
        ind_level_bc    = ind_level;
        % mmdd            = [];
        % dailyAnom_PMC   = [];

        %/ Search for complete data. If not existed, search for incomplete data
        datapath = []; flag_incomplete = 0;
        suffix = {'', '_incomplete'}; 
        for i = 1:length(suffix)
            filename = strcat(filename_full, num2str(year(t)), suffix{i}, '.nc');
            datapath = strcat(datafolder,filename);
            if isfile(datapath)
                if i == 2
                    flag_incomplete = 1;
                end
                break;
            end
        end
        if isempty(datapath)
            error('File %s not exist!')
        end

        %/ Check if 'expver' exists (i.e., contains preliminary data)
        try
            ncread(datapath, 'expver');
            flag_expver = 1;
            disp(['*** reading ', filename, ' (contains expver 5) ***']);
        catch 
            flag_expver = 0;
            disp(['*** reading ', filename, ' ***']);
        end

        %/ as lon lat somewhow may be different among each year's nc data (due to human error)
        lon = double(ncread(datapath,'longitude'));
        lat = double(ncread(datapath,'latitude'));
        ind_lon = find(stlon <= lon & lon <= edlon);
        ind_lat = find(stlat <= lat & lat <= edlat);
%         nlon = length(ind_lon);
%         nlat = length(ind_lat);

        time_ori          = double(ncread(datapath,'time'));
        date_yyyymmdd     = int64(str2num(datestr(dn1900 + time_ori/24+timeshift,'yyyymmdd')));
        date_yyyymmddHHMM = int64(str2num(datestr(dn1900 + time_ori/24+timeshift,'yyyymmddHHMM')));
        ind_time = find(ismember(date_yyyymmdd, case_date)); % let it auto sort. case_date in the output shall be sorted
        if isempty(ind_time) %/ e.g., a leap day in a non-leap year or incomplete nc data.
            warning(['[read_ERA5]: No data can be retrieved given the input "case_date", skip the year ', num2str(year(t)),'.'])
            continue
        end

        %/ If not all the subdaily steps are available on the last date, do not read them.
        ntime            = length(ind_time);
        if flag_incomplete && mod(ntime, StepsInADay) ~= 0 
            ind_rm           = (ntime-mod(ntime, StepsInADay)+1):ntime;
            ind_time(ind_rm) = [];  
        end
        ntime            = length(ind_time);  %/ Update
        case_date_eachYr = unique(date_yyyymmdd(ind_time)); %/ unique() to convert subdaily dates into daily dates
        ndays            = length(case_date_eachYr); %NOTE: case_date_eachYr ~= case_date cos we split the case year by year


        %/ List of acc fields excluded from the data length checking (for only 17 hrly time steps for acc fields on 19790101)
        exclusion_list = {'tp', 'e', 'ttr', 'str', 'tsr', 'tisr', 'ssr', 'slhf', 'sshf'};  

        %/ Check if the input StepsInADay match the data 
        %/ (if the data is not preliminary)
        if StepsInADay ~= ntime/ndays 
            if any(contains(varname_bc, exclusion_list)) && ismember(year(t), [1940, 1979])
                warning('Known issue: The first seven hrly fields are missing in ''acc'' field like %s on %d0101. Skip calling errors.', varname_bc, year(t));
            else 
                ncdisp(datapath);
                disp(date_yyyymmdd([1:4,end-3:end]));
                disp(date_yyyymmddHHMM([1:4,end-3:end]));
                error('The input StepsInADay (%.2f) ~= ntime/ndays (%.2f) in year %d. Consider Add %s to the exclusion list for acc fields!', StepsInADay, ntime/ndays, year(t), varname_bc);
            end
        end

        if flag_4D
            if StepsInADay == 1
                if flag_expver
                    disp(datapath)
                    error('[read_ERA5]: ''expver'' detected! Code not ready for processindg 4D preliminary data!')
                else
                    daily = double(ncread(datapath,varname_bc, ...
                                            [ind_lon(1)  ind_lat(1)  ind_level_bc(1) ind_time(1)],... % lon, lat, time
                                            [nlon        nlat        nlevel          ntime]));
                end
            else
                if flag_expver
                    disp(datapath)
                    error('[read_ERA5]: ''expver'' detected! Code not ready for processindg 4D preliminary data!')
                else
                    subdaily = double(ncread(datapath,varname_bc, ...
                                            [ind_lon(1)  ind_lat(1)  ind_level_bc(1) ind_time(1)],... % lon, lat, time
                                            [nlon        nlat        nlevel          ntime]));
                end

                if time_dim ~= length(size(subdaily))
                    error('[read_ERA5]: time_dim ~= length(size(subdaily))!');
                end
                daily = nan(nlon, nlat, nlevel, ndays);
                for i = 1:ndays
                    if isequal(datatype, 'acc')
                        daily(:,:,:,i) = sum(subdaily(:,:,:,(1:StepsInADay)+ (i-1)*StepsInADay),time_dim, 'omitnan');
                    else
                        daily(:,:,:,i) = mean(subdaily(:,:,:,(1:StepsInADay)+(i-1)*StepsInADay),time_dim, 'omitnan');
                    end
                end
            end
        else
            if StepsInADay == 1
                if flag_expver
                    daily = double(ncread(datapath,varname_bc, ...
                                        [ind_lon(1)  ind_lat(1)  1   ind_time(1)],... % lon, lat, expver, time
                                        [nlon        nlat        1   ntime]));

                    daily_expver5 = double(ncread(datapath,varname_bc, ...
                                        [ind_lon(1)  ind_lat(1)  2   ind_time(1)],... % lon, lat, expver, time
                                        [nlon        nlat        1   ntime]));
                    
                    %/ Replace NaN in ERA5 (expver 1) with data from ERA5T (expver 5) 
                    logi_nan = isnan(daily);
                    daily(logi_nan) = daily_expver5(logi_nan);
                    
                    if ~isempty(find(isnan(daily), 1))
                        warning('[read_ERA5]: daily still contains %d NaNs after replacement of preliminary data (expver 5). Check if this is expected!', length(find(isnan(daily))));
                    end
                else
                    daily = double(ncread(datapath,varname_bc, ...
                                        [ind_lon(1)  ind_lat(1)  ind_time(1)],... % lon, lat, time
                                        [nlon        nlat        ntime]));
                end
            else
                if flag_expver
                    subdaily = squeeze(double(ncread(datapath,varname_bc, ...
                                        [ind_lon(1)  ind_lat(1)  1   ind_time(1)],... % lon, lat, expver, time
                                        [nlon        nlat        1   ntime])));

                    subdaily_expver5 = squeeze(double(ncread(datapath,varname_bc, ...
                                        [ind_lon(1)  ind_lat(1)  2   ind_time(1)],... % lon, lat, expver, time
                                        [nlon        nlat        1   ntime])));


                    %/ Replace NaN in ERA5 (expver 1) with data from ERA5T (expver 5) 
                    logi_nan = isnan(subdaily);
                    subdaily(logi_nan) = subdaily_expver5(logi_nan);
                    
                    if ~isempty(find(isnan(subdaily), 1))
                        warning('[read_ERA5]: subdaily still contains %d NaNs after replacement of preliminary data (expver 5). Check if this is expected!', length(find(isnan(subdaily))));
                    end
                else
                    subdaily = double(ncread(datapath,varname_bc, ...
                                        [ind_lon(1)  ind_lat(1)  ind_time(1)],... % lon, lat, time
                                        [nlon        nlat        ntime]));
                end

                if time_dim ~= length(size(subdaily))
                    error('[read_ERA5]: time_dim ~= length(size(subdaily))!');
                end

                daily = nan(nlon, nlat, ndays);
                for i = 1:ndays
                    if any(strcmp(varname_bc, exclusion_list)) && ismember(year(t), [1940, 1979]) %/ 19790101 (1st date of satellite-era) or 19400101 (1st date of ERA5 extended data) has 17 time steps for acc. field
                        if i == 1
                            if isequal(datatype, 'acc')
                                daily(:,:, i) = sum(subdaily(:,:,1:17),time_dim, 'omitnan');
                            else
                                daily(:,:, i) = mean(subdaily(:,:,1:17),time_dim, 'omitnan');
                            end
                        else
                            if isequal(datatype, 'acc')
                                daily(:,:, i) = sum(subdaily(:,:,(1:StepsInADay)+ 17 + (i-2)*StepsInADay),time_dim, 'omitnan');
                            else
                                daily(:,:, i) = mean(subdaily(:,:,(1:StepsInADay)+ 17 + (i-2)*StepsInADay),time_dim, 'omitnan');
                            end
                        end
                    else
                        if isequal(datatype, 'acc')
                            daily(:,:, i) = sum(subdaily(:,:,(1:StepsInADay)+(i-1)*StepsInADay),time_dim, 'omitnan');
                        else
                            daily(:,:, i) = mean(subdaily(:,:,(1:StepsInADay)+(i-1)*StepsInADay),time_dim, 'omitnan');
                        end
                    end 
                end
            end
        end
        daily = daily*dataunitconv;
        subdaily = subdaily*dataunitconv;

        if any(ismember(select_field, {'subdaily'}))
            subdaily = reshape(subdaily, [], 1);
            subdaily_1D_list{t} = subdaily;
        end

        if any(ismember(select_field, {'yearly_sum'}))
            yearly_sum = sum(daily,time_dim, 'omitnan');
            yearly_sum = reshape(yearly_sum, [], 1); %convert to 1D, auto: lon*lat*time, 1
            yearly_sum_1D_list{t} = yearly_sum; % store into a cell array since sizes of col arrays are different.
        end

        if any(ismember(select_field, {'daily', 'daily_clim', 'PMC', 'dailyAnom_PMC', 'monthly', 'montlyAnom'})) %/ all these fields will need daily field.
            daily = reshape(daily, [], 1); %convert to 1D, auto: lon*lat*time, 1
            daily_1D_list{t} = daily; % store into a cell array since sizes of col arrays are different.
        end

        time_list{t} = time_ori(ind_time);

%         %/ clean up
%         daily       = [];
%         subdaily    = [];
%         yearly_sum  = [];
    end

    %/ Close parpool when not using
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        delete(poolobj)
    end

    output.subdaily     = [];
    output.daily        = [];
    output.yearly_sum   = [];
    output.time         = [];
    for t = 1:length(year)
        %/ output subdaily only upon request.
        if any(ismember(select_field, {'subdaily'}))
            if flag_4D
                output.subdaily = cat(time_dim, output.subdaily, reshape(subdaily_1D_list{t}, nlon, nlat, nlevel, []));
            else
                output.subdaily = cat(time_dim, output.subdaily, reshape(subdaily_1D_list{t}, nlon, nlat, []));
            end
        end
        if any(ismember(select_field, {'daily', 'daily_clim', 'PMC', 'dailyAnom_PMC', 'monthly', 'montlyAnom'}))
            if flag_4D
                output.daily = cat(time_dim, output.daily, reshape(daily_1D_list{t}, nlon, nlat, nlevel, []));
            else
                output.daily = cat(time_dim, output.daily, reshape(daily_1D_list{t}, nlon, nlat, []));
            end
        end
        if any(ismember(select_field, {'yearly_sum'}))
            if flag_4D
                output.yearly_sum = cat(time_dim, output.yearly_sum, reshape(yearly_sum_1D_list{t}, nlon, nlat, nlevel, []));
            else
                output.yearly_sum = cat(time_dim, output.yearly_sum, reshape(yearly_sum_1D_list{t}, nlon, nlat, []));
            end
        end
        output.time = cat(1, output.time, time_list{t});
    end

    %/ clean up
    clear subdaily_1D_list 
    clear daily_1D_list
    clear yearly_sum_1D_list

    %/ NOTE: For ERA5, 2nd January 2017 time = 00 will give you total precipitation data to cover 23 - 24 UTC for 1st January 2017
    %/       https://confluence.ecmwf.int/display/CKB/ERA5%3A+How+to+calculate+daily+total+precipitation
    output.date_UTC_AllYr = datetime(dn1900 + output.time/24+timeshift, 'ConvertFrom','datenum', 'Format', 'yyyyMMddHHmm', 'TimeZone', 'UTC'); %/ UTC 
    output.date_HKT_AllYr = datetime(output.date_UTC_AllYr, 'TimeZone', 'Asia/Hong_Kong'); %/ If the input data are character vectors or strings that include a time zone, then the datetime function converts all values to the specified time zone.

    output.date_yyyymmdd_AllYr = unique(str2num(datestr(dn1900 + output.time/24+timeshift,'yyyymmdd'))); %/ double
    output.date_mmdd_AllYr = str2num(datestr(dn1900 + output.time/24+timeshift,'mmdd'));
    output.date_mmdd_AllYr = output.date_mmdd_AllYr(1:StepsInADay:end);
    output.date_mmdd_OneYr = unique(output.date_mmdd_AllYr);
    nday = length(output.date_mmdd_OneYr);

%     disp(length(output.date_mmdd_AllYr))
%     disp(length(output.date_mmdd_OneYr))
%     disp(output.date_mmdd_OneYr)

    if any(ismember(select_field, {'daily_clim', 'PMC', 'dailyAnom_PMC'}))
        if flag_4D
            for i = 1:nday
               day_ind = find(output.date_mmdd_AllYr == output.date_mmdd_OneYr(i));
    %            disp(day_ind)
               output.daily_clim(:,:,:,i) = mean(output.daily(:,:,:,day_ind),time_dim, 'omitnan');
            end
        else
            for i = 1:nday
                day_ind = find(output.date_mmdd_AllYr == output.date_mmdd_OneYr(i));
                output.daily_clim(:,:,i) = mean(output.daily(:,:,day_ind),time_dim, 'omitnan');
            end
        end

        % Pentad Moving Mean Clim (PMC)
        if output.date_mmdd_OneYr(1) == 101 && output.date_mmdd_OneYr(end) == 1231
            %then we can append 2 days to each of the two ends for 5-day moving mean
            if flag_4D
                a = cat(time_dim, output.daily_clim(:,:,:,end-1:end),...
                                                          output.daily_clim(:,:,:,:),...
                                                          output.daily_clim(:,:,:,1:2));
                output.PMC = movmean(a,5,time_dim); % movmean(X,movingdays,dim)
                output.PMC = output.PMC(:,:,:,3:end-2);
            else
                a = cat(time_dim, output.daily_clim(:,:,end-1:end),...
                                                          output.daily_clim(:,:,:),...
                                                          output.daily_clim(:,:,1:2));
                output.PMC = movmean(a,5,time_dim); % movmean(X,movingdays,dim)
                output.PMC = output.PMC(:,:,3:end-2);
            end
        else
            output.PMC = movmean(output.daily_clim,5,time_dim); % movmean(X,movingdays,dim)
        end
    end

    % daily anom on PMC
    if any(ismember(select_field, {'dailyAnom_PMC'}))
        if flag_4D
            [n1, n2, n3, n4] = size(output.daily);
            output.dailyAnom_PMC = nan(n1, n2, n3, n4);
            for i = 1:nday
                day_ind = find(output.date_mmdd_AllYr == output.date_mmdd_OneYr(i));
                output.dailyAnom_PMC(:,:,:,day_ind) = ...
                                output.daily(:,:,:,day_ind) - output.PMC(:,:,:,i);
            end
        else
            [n1, n2, n3] = size(output.daily);
            output.dailyAnom_PMC = nan(n1, n2, n3);
            for i = 1:nday
                day_ind = find(output.date_mmdd_AllYr == output.date_mmdd_OneYr(i));
                output.dailyAnom_PMC(:,:,day_ind) = ...
                                output.daily(:,:,day_ind) - output.PMC(:,:,i);
            end
        end
        output.daily = []; %/ spare memory
    end

    %/ monthly
    if any(ismember(select_field, {'monthly', 'monthlyAnom'}))
        output.date_yyyymm_AllYr = str2num(datestr(dn1900 + output.time/24+timeshift,'yyyymm')); % 14610x1, do NOT use str2double -> it causes bug!
        output.date_yyyymm_AllYr = output.date_yyyymm_AllYr(1:StepsInADay:end); % 14610x1
%         disp(output.date_yyyymm_AllYr)
    %     output.date_yyyymm_AllYr = [reshape(repmat(year, 12, 1), 1, [])*100 + repmat(1:12, 1, length(year))]'; % 480x1
        output.date_mm_AllYr = repmat(1:12, 1, length(year));  % 480x1
        
        [n1, n2, n3, ~] = size(output.daily);
        if flag_4D
            output.monthly = nan(n1, n2, n3, length(year)*12);
            for t = 1:length(year)
                for m = 1:12
                    day_ind = find(output.date_yyyymm_AllYr == year(t)*100 + m);
                    if isempty(day_ind)  
                        disp(output.date_yyyymm_AllYr)
                        disp(year(t)*100 + m)
                        error('empty day_ind!'); 
                    end
                    if isequal(datatype, 'acc')
                        output.monthly(:,:,:,m + 12*(t-1)) = sum(output.daily(:,:,:,day_ind),time_dim, 'omitnan');
                    else      
                        output.monthly(:,:,:,m + 12*(t-1)) = mean(output.daily(:,:,:,day_ind),time_dim, 'omitnan');
                    end
                end
            end
        else
            output.monthly = nan(n1, n2, length(year)*12);
            for t = 1:length(year)  
                for m = 1:12
                    day_ind = find(output.date_yyyymm_AllYr == year(t)*100 + m);
                    if isempty(day_ind)  
                        disp(output.date_yyyymm_AllYr)
                        disp(year(t)*100 + m)
                        error('empty day_ind!'); 
                    end
                    if isequal(datatype, 'acc')
                        output.monthly(:,:,m + 12*(t-1)) = sum(output.daily(:,:,day_ind),time_dim, 'omitnan');
                    else     
                        output.monthly(:,:,m + 12*(t-1)) = mean(output.daily(:,:,day_ind),time_dim, 'omitnan');
                    end
                end
            end
        end
        output.daily = [];  %/ spare memory

        %/ monthly clim
        [n1, n2, n3, ~] = size(output.monthly);
        if flag_4D
            output.monthly_clim = nan(n1, n2, n3, 12);
            for m = 1:12
               ind = find(output.date_mm_AllYr == m);
               output.monthly_clim(:,:,:,m) = mean(output.monthly(:,:,:,ind),time_dim, 'omitnan');
            end
        else
            output.monthly_clim = nan(n1, n2, 12);
            for m = 1:12
               ind = find(output.date_mm_AllYr == m);
               output.monthly_clim(:,:,m) = mean(output.monthly(:,:,ind),time_dim, 'omitnan');
            end
        end

        %/ monthly Anom
        output.monthlyAnom = nan(size(output.monthly));
        if flag_4D
            output.monthlyAnom = nan(size(output.monthly));
            for m = 1:12
               day_ind = find(output.date_mm_AllYr == m);
               if isempty(day_ind)  error('empty day_ind!'); end
               output.monthlyAnom(:,:,:,day_ind) = output.monthly(:,:,:,day_ind) - output.monthly_clim(:,:,:,m);
            end
        else
            output.monthlyAnom = nan(size(output.monthly));
            for m = 1:12
               day_ind = find(output.date_mm_AllYr == m);
               if isempty(day_ind)  error('empty day_ind!'); end
               output.monthlyAnom(:,:,day_ind) = output.monthly(:,:,day_ind) - output.monthly_clim(:,:,m);
            end
        end
        output.date_yyyymmdd_AllYr = int64(unique(floor(output.date_yyyymmdd_AllYr/1e2))*1e2+1); %/ old, but keep it for now
        output.date_yyyymm         = int64(unique(floor(output.date_yyyymmdd_AllYr/1e2)));       %/ new 
    end

    if flag_4D
        %/ Finally, squeeze out any singluar level dim (when selecting a single level)
        output.(select_field{1}) = squeeze(output.(select_field{1}));
    end

    %/ remove all other fields in the struct data (i.e. output) that are not selected.
    flds = fieldnames(output);
    omitfields = setdiff(flds, select_field); 
    for r = 1:length(omitfields)
        if isfield(output, omitfields{r})
            output = rmfield(output, omitfields{r});
        end
    end

    % if NonStruct  %/ Then do not output the struct data (i.e. output)
    %     %/ only returns 4 fields.
    %     if ~any(ismember(select_field{1}, {'subdaily', 'daily', 'yearly_sum', 'daily_clim', 'PMC', 'dailyAnom_PMC', 'monthly', 'monthlyAnom'}))
    %         error('Rewrite your list of "select_field", put the desired field in the 1st position!')
    %     else
    %         fld_data = output.(select_field{1});
    %     end
    %     fld_lon = output.lon;
    %     fld_lat = output.lat;
    %     fld_date_yyyymmdd_AllYr = int64(output.date_yyyymmdd_AllYr);
    %     output = [];
    % end
    % toc
end