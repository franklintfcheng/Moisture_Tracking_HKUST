function [output, fld_data, fld_lon, fld_lat, fld_date] = read_CERA(varargin)

% create a set of valid parameters and their default value
pnames = {'time_dim',          'datafolder',  'filename_prefix', 'filename_part',...
          'dataname',  'dataname_callfile', 'varname',     'dataunitconv',    'datatype',...
          'allyears',  'stlon',             'edlon',       'stlat',           'edlat',...
          'stlevel',   'edlevel',           'slct_level',  'case_date',       'StepsInADay',...
          'timeshift', 'select_field',       'NumWorkers',  'NonStruct',       'fc_steps'};

dflts  = cell(length(pnames), 1);

%/ parse function arguments
         [ time_dim,          datafolder,  filename_prefix, filename_part,...
           dataname,  dataname_callfile, varname,     dataunitconv,    datatype,...
           allyears,  stlon,             edlon,       stlat,           edlat,...
           stlevel,   edlevel,           slct_level,  case_date,       StepsInADay,...
           timeshift, select_field,       NumWorkers,  NonStruct,       fc_steps ] ...
               = internal.stats.parseArgs(pnames, dflts, varargin{:});
    %%
    %/=====================================================================
    %/       Author: Franklin Cheng
    %/ Last Updated: 12 Feb 2024
    %/
    %/ This is a function to read *all* or *partial* continuous time entries in CERA-20C data /% 
    %/
    %/ Available output fields: 
    %/             lon, lat, time, daily
    %/             date_yyyymmdd_AllYr, date_mmdd_AllYr, date_mmdd_OneYr
    %/             daily_clim, PMC, dailyAnom_PMC
    %/             monthly, monthlyAnom
    %/
    %/ *** WARNING ***:
    %/ CERA fc data starts from 18 UTC. To obtain the daily data for day t:
    %/      daily(t) = subdaily(t, 0~18 UTC) + subdaily(t, 18~24 UTC)
    %/      daily(t) = subdaily(t-1, +24h step) - subdaily(t-1, +6h step) + subdaily(t, +6 step)
    %/  ==> The first day of each year will be missing for an acc variable!
    %/
    %/=====================================================================
    
    output = [];
    if numel(num2str(case_date(1))) == 8
        year = unique(round(case_date/1e4));
    else
        error('incorrect format of case_date!');
    end
    
    if isempty(stlon)     stlon = 0;        end 
    if isempty(edlon)     edlon = 360;      end
    if isempty(stlat)     stlat = -90;      end
    if isempty(edlat)     edlat = 90;       end
    if isempty(timeshift) timeshift = 0;    end
    % disp(year);
    if time_dim ~= 3 && time_dim ~= 4             error('Only time_dim = 3 or 4 is allowed.'); end
    if ~any(ismember(datatype, {'ins', 'acc'}))   error('Specify datatype! ins or acc?');      end
    % if isempty(allyears) allyears = 1971:2010; end %/ allyears is only used in mode2.
    % disp(allyears)
    % if isempty(NoOfWorkers) NoOfWorkers = Inf; end  %/ if NoOfWorkers = 0 --> no parallel computing.
    
    %/ Broadcast var /%
    dataname_callfile_bc = dataname_callfile;
    % dataname_bc          = dataname;
    varname_bc           = varname;
    
    %/ set timeshift = 0 if not specified or datatype = 'ins'.
    if isempty(timeshift)
        timeshift = 0;
    else
        if ismember(datatype, {'ins'})
            timeshift = 0;
        end
    end
    % disp(['NOTE: timeshift = ', num2str(timeshift), ' days.']);
    dn1900 = datenum(1900,1,1);
    
    %/ Read basic info (e.g., lat, lon, level, time) /%
    styear = year(1);
    if      isempty(filename_prefix)  error('ERROR: filename_prefix is not specified! Note: it reads "filename_prefix + dataname_callfile + filename_part + year(t) + .nc"');
    else if isempty(filename_part)    error('ERROR: filename_part is not specified!   Note: it reads "filename_prefix + dataname_callfile + filename_part + year(t) + .nc"');
    else
        filename_full = strcat(filename_prefix, dataname_callfile_bc,filename_part);
    end
    
    filename = strcat(filename_full, num2str(styear),'.nc');
    lon = double(ncread(strcat(datafolder,filename),'longitude'));
    lat = double(ncread(strcat(datafolder,filename),'latitude'));
    ind_lon = find(stlon <= lon & lon <= edlon);
    ind_lat = find(stlat <= lat & lat <= edlat);
    nlon = length(ind_lon);
    nlat = length(ind_lat);
    
    output.lon = lon(ind_lon); %update 
    output.lat = lat(ind_lat); %update
    
    if time_dim == 4
        output.level = double(ncread(strcat(datafolder,filename),'level'));  
        if ~isempty(slct_level) %let slct_level overwirtes stlevel and edlevel if it's given.
            ind_level = find(ismember(output.level, slct_level));
        else
            ind_level = find(stlevel <= output.level & output.level <= edlevel);
        end
        
        nlevel = length(ind_level);  
        output.level = output.level(ind_level);
    else
        ind_level = []; nlevel = []; %/ still need to specify them as broadcast vars
    end
    
    subdaily_1D_list = cell(length(year),1);
    daily_1D_list = cell(length(year),1);
    yearly_sum_1D_list = cell(length(year),1);
    % dailyAnom_PMC_1D_list = cell(length(year),1);
    time_list = cell(length(year),1);
    
    if isempty(gcp('nocreate')) && ~isempty(NumWorkers) %/ if set worker number
        parpool('local', NumWorkers) %/ use process-based parpool (threads-based parpool is only available after R2020a :((
    end
    
    parfor t = 1:length(year)
    % for t = 1:length(year)
    %     disp(t)
        daily = []; %/ KEEP THIS. Initialize it to avoid warnings for temporary var generated in the parfor loop.
        subdaily = [];
        % mmdd = [];
        % dailyAnom_PMC = [];
    
        filename = strcat(filename_full, num2str(year(t)),'.nc');
    %     ncdisp(strcat(datafolder,filename))
        disp(['*** reading ', filename, ' ***']);
            
        time_ori          = double(ncread(strcat(datafolder,filename),'time'));
        date_yyyymmdd     = str2num(datestr(dn1900 + time_ori/24+timeshift,'yyyymmdd'));
        ind_date_yyyymmdd = find(ismember(date_yyyymmdd, case_date)); % let it auto sort. case_date in the output shall be sorted
        
        if isempty(ind_date_yyyymmdd) %/ e.g., a leap day in a non-leap year or incomplete nc data.
            warning(['(read_CERA): No data can be retrived given the input "case_date", skip the year ', num2str(year(t)),'.'])
            continue
        end
        
        case_date_eachYr = unique(date_yyyymmdd(ind_date_yyyymmdd)); %unique() for houly data on the same day
        ndays            = length(case_date_eachYr); %NOTE: case_date_eachYr ~= case_date cos we split the case year by year
    
        ndate = length(ind_date_yyyymmdd);
        if StepsInADay ~= ndate/ndays % check if the input StepsInADay match the data
            fprintf('ndate/ndays=%.5f \n', ndate/ndays)
            error('StepsInADay does not match ndate/ndays, check if the data is complete.');
        end
    
        if time_dim == 4
            if StepsInADay == 1
                daily = double(ncread(strcat(datafolder,filename),varname_bc, ...
                                        [ind_lon(1)  ind_lat(1)  ind_level(1) ind_date_yyyymmdd(1)],... % lon, lat, time
                                        [nlon        nlat        nlevel       ndate]));
            else
                subdaily = double(ncread(strcat(datafolder,filename),varname_bc, ...
                                        [ind_lon(1)  ind_lat(1)  ind_level(1) ind_date_yyyymmdd(1)],... % lon, lat, time
                                        [nlon        nlat        nlevel       ndate]));
    
                daily = nan(nlon, nlat,nlevel, ndays);
                for i = 1:ndays     %/ normally no 4d data are of accum type (?)
                    if any(strcmp(datatype, {'acc'})) 
                        error('Check if the 4D data is really an accumulate var!!!\n')
                    else
                        daily(:,:,:, i) = mean(subdaily(:,:,:,(1:StepsInADay)+(i-1)*StepsInADay),time_dim);
                    end
                end
            end
        else
            if StepsInADay == 1
                daily = double(ncread(strcat(datafolder,filename),varname_bc, ...
                                    [ind_lon(1)  ind_lat(1)  ind_date_yyyymmdd(1)],... % lon, lat, time
                                    [nlon        nlat        ndate]));
            else
                subdaily = double(ncread(strcat(datafolder,filename),varname_bc, ...
                                    [ind_lon(1)  ind_lat(1)  ind_date_yyyymmdd(1)],... % lon, lat, time
                                    [nlon        nlat        ndate]));               %/ get one more time step for computation.
                                
                if any(strcmp(datatype, {'acc'}))
                    if isequal(fc_steps, [6,24])   
                        daily = nan(nlon, nlat, ndays-1); %/ ndays-1 cos no data to compute daily data for the last day.
                        for i = 1:ndays-1
                            %/ NOTE: CERA fc data starts from 18 UTC. To obtain the daily data for day t:
                            %/ daily(t) = subdaily(t, 0~18 UTC) + subdaily(t, 18~24 UTC)
                            %/ daily(t) = subdaily(t-1, +24h step) - subdaily(t-1, +6h step) + subdaily(t, +6 step)
    
                            daily(:,:, i) =   subdaily(:,:, StepsInADay*i) - subdaily(:,:, StepsInADay*i-1)...
                                            + subdaily(:,:, StepsInADay*i+1);
                        end
                        
                    elseif isequal(fc_steps, [3:3:24])
                        subdaily_correct = nan(nlon, nlat, ndate);  
                        for i = 1:ndate
                            if mod(i-1, 8) == 0
                                subdaily_correct(:,:,i) = subdaily(:,:,i);
                            else
                                subdaily_correct(:,:,i) = subdaily(:,:,i) - subdaily(:,:,i-1);
                            end
                        end
                        subdaily = subdaily_correct;   %/ now the subdaily prcp/evap is correct.
                    else
                        error('!!! Check the fc_steps !!!'); 
                    end
                else
                    for i = 1:ndays
                        daily(:,:,i) = mean(subdaily(:,:,(1:StepsInADay)+(i-1)*StepsInADay),time_dim);
                    end
                end
            end
        end
        daily    = daily*dataunitconv;
        subdaily = subdaily*dataunitconv;
    
        if any(ismember(select_field, {'subdaily'}))
            subdaily = reshape(subdaily, [], 1);
            subdaily_1D_list{t} = subdaily;
        end
    
        if any(ismember(select_field, {'yearly_sum'}))
            yearly_sum = sum(daily,time_dim);
            yearly_sum = reshape(yearly_sum, [], 1); %convert to 1D, auto: lon*lat*time, 1
            yearly_sum_1D_list{t} = yearly_sum; % store into a cell array since sizes of col arrays are different.
        end
    
        if any(ismember(select_field, {'daily', 'daily_clim', 'PMC', 'dailyAnom_PMC', 'monthly', 'montlyAnom'})) %/ all these fields will need daily field.
            daily = reshape(daily, [], 1); %convert to 1D, auto: lon*lat*time, 1
            daily_1D_list{t} = daily; % store into a cell array since sizes of col arrays are different.
        end
    
        if any(strcmp(datatype, {'acc'})) && isequal(fc_steps, [6,24])
            time_list{t} = time_ori(ind_date_yyyymmdd(1:end-length(fc_steps))); %/ remove the last day (2 time steps) if handling fc data.
        else
            time_list{t} = time_ori(ind_date_yyyymmdd);
        end
    end

    %/ Close parpool when not using
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        delete(poolobj)
    end
     
    output.subdaily = [];
    output.daily = [];
    output.yearly_sum = [];
    output.time = [];
    for t = 1:length(year)
        %/ output subdaily only upon request.
        if any(ismember(select_field, {'subdaily'}))
            if time_dim == 4
                output.subdaily = cat(time_dim, output.subdaily, reshape(subdaily_1D_list{t}, nlon, nlat, nlevel, []));
            else
                output.subdaily = cat(time_dim, output.subdaily, reshape(subdaily_1D_list{t}, nlon, nlat, []));
            end
        end
        
        if any(ismember(select_field, {'daily', 'daily_clim', 'PMC', 'dailyAnom_PMC', 'monthly', 'montlyAnom'}))
            if time_dim == 4
                output.daily = cat(time_dim, output.daily, reshape(daily_1D_list{t}, nlon, nlat, nlevel, []));
            else
                output.daily = cat(time_dim, output.daily, reshape(daily_1D_list{t}, nlon, nlat, []));
            end
        end
        
        if any(ismember(select_field, {'yearly_sum'}))
            if time_dim == 4
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
    
    
    if ~isempty(fc_steps)
        output.date_yyyymmddHHMMSS = datestr(dn1900 + output.time/24+timeshift,'yyyymmdd HH:MM:SS');
    %         disp(output.date_yyyymmddHHMMSS(1:10,:))
    %         disp(output.date_yyyymmddHHMMSS(end-10:end, :))
    end
    output.date_yyyymmdd_AllYr = unique(str2num(datestr(dn1900 + output.time/24+timeshift,'yyyymmdd')));
    output.date_mmdd_AllYr     = str2num(datestr(dn1900 + output.time/24+timeshift,'mmdd'));
    output.date_mmdd_AllYr     = output.date_mmdd_AllYr(1:StepsInADay:end);
    output.date_mmdd_OneYr     = unique(output.date_mmdd_AllYr);
    %     disp(length(output.date_mmdd_AllYr))
    %     disp(length(output.date_mmdd_OneYr))
    %     disp(output.date_mmdd_OneYr)
    
    nday = length(output.date_mmdd_OneYr);
    if any(ismember(select_field, {'daily_clim', 'PMC', 'dailyAnom_PMC'}))
        if time_dim == 4
            for i = 1:nday
               day_ind = find(output.date_mmdd_AllYr == output.date_mmdd_OneYr(i));
    %            disp(day_ind)
               output.daily_clim(:,:,:,i) = mean(output.daily(:,:,:,day_ind),time_dim);
            end
        else
            for i = 1:nday
                day_ind = find(output.date_mmdd_AllYr == output.date_mmdd_OneYr(i));
                output.daily_clim(:,:,i) = mean(output.daily(:,:,day_ind),time_dim);
            end
        end
    
        % Pentad Moving Mean Clim (PMC)
        if output.date_mmdd_OneYr(1) == 101 && output.date_mmdd_OneYr(end) == 1231
            %then we can append 2 days to each of the two ends for 5-day moving mean
            if time_dim == 4
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
        if time_dim == 4
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
    end
    
    %/ monthly
    if any(ismember(select_field, {'monthly', 'monthlyAnom'}))
       
        output.date_yyyymm_AllYr = str2num(datestr(dn1900 + output.time/24+timeshift,'yyyymm')); % 14610x1, do NOT use str2double -> it causes bug!
        output.date_yyyymm_AllYr = output.date_yyyymm_AllYr(1:StepsInADay:end); % 14610x1
    %         disp(output.date_yyyymm_AllYr)
        
    %     output.date_yyyymm_AllYr = [reshape(repmat(year, 12, 1), 1, [])*100 + repmat(1:12, 1, length(year))]'; % 480x1
        output.date_mm_AllYr = repmat(1:12, 1, length(year));  % 480x1
        [n1, n2, n3, ~] = size(output.daily);
        
        if time_dim == 4
            output.monthly = nan(n1, n2, n3, length(year)*12);
            for t = 1:length(year)
                for m = 1:12
                    day_ind = find(output.date_yyyymm_AllYr == year(t)*100 + m);
        %            disp(day_ind)
                    if any(strcmp(datatype, {'acc'}))
                        output.monthly(:,:,:,m + 12*(t-1)) = sum(output.daily(:,:,:,day_ind),time_dim);
                    elseif any(strcmp(datatype, {'ins'}))      
                        output.monthly(:,:,:,m + 12*(t-1)) = mean(output.daily(:,:,:,day_ind),time_dim);
                    end
                end
            end
        else
            output.monthly = nan(n1, n2, length(year)*12);
            for t = 1:length(year)
                for m = 1:12
                    day_ind = find(output.date_yyyymm_AllYr == year(t)*100 + m);
    %                 disp(day_ind)
                    if any(strcmp(datatype, {'acc'})) 
                        output.monthly(:,:,m + 12*(t-1)) = sum(output.daily(:,:,day_ind),time_dim);
                    elseif any(strcmp(datatype, {'ins'}))          
                        output.monthly(:,:,m + 12*(t-1)) = mean(output.daily(:,:,day_ind),time_dim);
                    end
                end
            end
        end
        
        %/ monthly clim
        [n1, n2, n3, ~] = size(output.monthly);
        if time_dim == 4
            output.monthly_clim = nan(n1, n2, n3, 12);
            for m = 1:12
               ind = find(output.date_mm_AllYr == m);
               output.monthly_clim(:,:,:,m) = mean(output.monthly(:,:,:,ind),time_dim);
            end
        else
            output.monthly_clim = nan(n1, n2, 12);
            for m = 1:12
               ind = find(output.date_mm_AllYr == m);
               output.monthly_clim(:,:,m) = mean(output.monthly(:,:,ind),time_dim);
            end
        end
        
        %/ monthly Anom
        output.monthlyAnom = nan(size(output.monthly));
        if time_dim == 4
            output.monthlyAnom = nan(size(output.monthly));
            for m = 1:12
               ind = find(output.date_mm_AllYr == m);
               output.monthlyAnom(:,:,:,ind) = output.monthly(:,:,:,ind) - output.monthly_clim(:,:,:,m);
            end
        else
            output.monthlyAnom = nan(size(output.monthly));
            for m = 1:12
               ind = find(output.date_mm_AllYr == m);
               output.monthlyAnom(:,:,ind) = output.monthly(:,:,ind) - output.monthly_clim(:,:,m);
            end
        end
        
        output.date_yyyymm = int64(unique(floor(output.date_yyyymmdd_AllYr/1e2))*1e2+1); 
        output = rmfield(output, 'date_yyyymmdd_AllYr');
    end
    
    %/ remove all other fields in the struct data (i.e. output) that are not selected.
    flds = fieldnames(output);
    omitfields = setdiff(flds, select_field); 
    for r = 1:length(omitfields)
        if isfield(output, omitfields{r})
            output = rmfield(output, omitfields{r});
        end
    end
    
    if NonStruct == 1 %/ Then do not output the struct data (i.e. output)
        %/ only returns 4 fields.
        if ~any(ismember(select_field{1}, {'subdaily', 'daily', 'yearly_sum', 'daily_clim', 'PMC', 'dailyAnom_PMC', 'monthly', 'monthlyAnom'}))
            error('Rewrite your list of "select_field", put the wanted data field in the 1st position!')
        else
            fld_data = output.(select_field{1});
        end
        fld_lon = output.lon;
        fld_lat = output.lat;
        if any(ismember(select_field, {'monthly', 'monthlyAnom'}))
            fld_date = int64(output.date_yyyymmdd);
        else
            fld_date = int64(output.date_yyyymmdd_AllYr);
        end
        output = [];
    end

end