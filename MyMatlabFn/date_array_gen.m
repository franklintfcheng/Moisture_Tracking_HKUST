function [output_dates, str_dates] = date_array_gen(varargin)
    
    %/ create a set of valid parameters and their default value
    pnames = {'year_list', 'st_year', 'st_month', 'st_day', 'ed_year', 'ed_month', 'ed_day',...
              'st_hr', 'st_min', 'ed_hr', 'ed_min',...
              'season', 'dt_slct_mo', 'dt_slct_hr', 'dt_slct_min', 'output_date_format',...
              'noleap', 'skip_the_incomplete', 'calendar'}; 
    
    dflts  = cell(length(pnames), 1);
    
    [          year_list,   st_year,   st_month,   st_day,   ed_year,  ed_month,   ed_day,...
               st_hr,   st_min,   ed_hr,   ed_min,...
               season,  dt_slct_mo,   dt_slct_hr,   dt_slct_min,   output_date_format,...
               noleap,   skip_the_incomplete,   calendar] ...
                   = internal.stats.parseArgs(pnames, dflts, varargin{:});

%%
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 20 Feb 2023
    %/
    %/ Description: This function generate date arrays
    %/         'season': 
    %/                   1 to 4: MAM / JJA / SON / D(0)JF(1) 
    %/                        5: Apr-Sep
    %/                        6: Oct(0)-Mar(1)
    %/                        7: JFD(0)
    %/                        8: non-JJA
    %/                        9: May-Oct
    %/                       10: Nov(0)-Apr(1)
    %/                       11: May-Sep
    %/
    %/         'dt_slct': time internal (hr), if not given, assume 24 (daily)
    %/     'dt_slct_min': time internal (min)
    %/=====================================================================
    

    % str_mth = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',...
               % 'Sep', 'Oct', 'Nov', 'Dec'};

    str_season = {'MAM', 'JJA', 'SON', 'DJF',...
                 'AMJJAS','ONDJFM','JFD0','nonJJA','MJJASO','NDJFMA','MJJAS'};

    if ~isempty(season)
        str_dates = strcat('_', str_season{season});
    elseif isempty(st_month) && isempty(ed_month)
        str_dates = '';  %/ asumme to extract all dates
    else
        str_dates = sprintf('_%d%d-%d%d', st_month, st_day, ed_month, ed_day);
    end

    if isempty(output_date_format)  output_date_format  = 'yyyymmdd'; end
    if isempty(skip_the_incomplete) skip_the_incomplete = 0;          end  %/ Turn it off by default, since it depends on the upper-level coding script. Turn it on may cause bugs.
    
    % %/ Full half-hourly datetime array
    % date_st_dt = datetime(2009,   3,     1, 'format','yyyyMMdd HH:mm:ss');
    % date_ed_dt = datetime(2009,   5,  31+1, 'format','yyyyMMdd HH:mm:ss'); %/ use eomday() to find the last day of the month
    % dates_dt = [date_st_dt:hours(dt_slct):date_ed_dt]';  
    % dates_dt(end) = [];  
    
    if isempty(dt_slct_hr)  dt_slct_hr = 24;       end
    if isempty(st_day)      st_day     = 1;        end
    if isempty(ed_day)      ed_day     = 31;       end
    if isempty(st_hr)       st_hr      = 0;        end
    if isempty(st_min)      st_min     = 0;        end
    if isempty(ed_hr)       ed_hr      = 0;        end
    if isempty(ed_min)      ed_min     = 0;        end

    %/ If the period is not specified, retrieve the entire period of the year
    if isempty(season) && isempty(st_month) && isempty(ed_month)
        st_month =  1; 
        ed_month = 12; 
    end

    %/ Retrieve dates based on 'st_year', 'ed_year' directly if given 
    if ~isempty(st_year) && ~isempty(ed_year) && ~isempty(st_month) && ...
       ~isempty(ed_month) && ~isempty(st_day) && ~isempty(ed_day)

        flag_contin_dates = 1;
        years_unique = 1;  %/ dummy
    else
        flag_contin_dates = 0;
        if isempty(year_list)
            error('year_list is missing!');
        end
        years_unique = unique(year_list, 'stable');
        if ~isequal(year_list, years_unique)
            warning('[date_array_gen]: Duplicated year_list detected and removed!')
        end
    end
    
    %/ Special treatment of julian calendar 
    if isequal(calendar, 'julian')
        fprintf('*** As per the requirement, generating julian dates... ***')
    
        if flag_contin_dates
            error('code not ready!')
        end

        %/ Generate a 365-day datetime first
        t1 = datetime('20230101', 'InputFormat','yyyyMMdd');
        t2 = datetime('20231231', 'InputFormat','yyyyMMdd');
        tv = (t1:caldays(1):t2)';
        date_365_yyyymmdd = datetime2int(tv, 'yyyymmdd');
        date_365_mmdd = mod(date_365_yyyymmdd, 1e4);
    
        ind = find(date_365_mmdd == 228);
        date_366_mmdd = [date_365_mmdd(1:ind); 229; date_365_mmdd(ind+1:end)]; 
    
        %/ Unlike standard calendar, julian considers a year is a leap year if it is divisible by 4, 
        %/ even if it is also divisible by 100.
        %/ So, the goal here is to get those 'excessive' leap days with the given year_list
        output_dates = [];
        for t = 1:length(year_list)
            if mod(year_list(t), 4) == 0  %/ Julian leap year
                dates_eachyear = year_list(t)*1e4 + date_366_mmdd;
            else
                dates_eachyear = year_list(t)*1e4 + date_365_mmdd;
            end
            output_dates = cat(1, output_dates, dates_eachyear);
        end
        return;  
    end

    output_dates = []; flag_break = 0;
    for t = 1:length(years_unique)   

        %/ Retrieve dates based on 'st_year', 'ed_year' directly if given 
        if flag_contin_dates
            date_st_dt = datetime(st_year,   st_month,  st_day,   st_hr, st_min, 0, 'format','yyyyMMdd HH:mm:ss');
            date_ed_dt = datetime(ed_year,   ed_month,  ed_day+1, ed_hr, ed_min, 0, 'format','yyyyMMdd HH:mm:ss'); %/ use eomday() to find the last day of the month
            
            if ~isempty(dt_slct_mo)       dates_dt = (date_st_dt:calmonths(dt_slct_mo):date_ed_dt)';
            elseif ~isempty(dt_slct_min)  dates_dt = (date_st_dt:minutes(dt_slct_min):date_ed_dt)';
            else                          dates_dt = (date_st_dt:hours(dt_slct_hr):date_ed_dt)';   end
            dates_dt(end) = [];
            flag_break = 1;  %/ Break the year loop after converting dates_dt into output_dates
        
        elseif ~isempty(season)
            if ismember(season, [1, 2, 3, 5, 9, 11])  %/ for these seasons do not cross two year_list.
                if season == 1
                    st_month = 3;  ed_month = 5;   %/ MAM
                elseif season == 2
                    st_month = 6;  ed_month = 8;   %/ JJA
                elseif season == 3
                    st_month = 9;  ed_month = 11;  %/ SON
                elseif season == 5
                    st_month = 4;  ed_month = 9;   %/ AMJJAS (Apr-Sep)
                elseif season == 9
                    st_month = 5;  ed_month = 10;  %/ MJJASO (May-Oct)
                elseif season == 11
                    st_month = 5;  ed_month = 9;   %/ MJJAS  (May-Sep)
                else
                    error('something wrong here!');
                end
                st_day = 1; ed_day = eomday(years_unique(t), ed_month);
                
                date_st_dt = datetime(years_unique(t),   st_month,  st_day,   st_hr, st_min, 0, 'format','yyyyMMdd HH:mm:ss');
                date_ed_dt = datetime(years_unique(t),   ed_month,  ed_day+1, ed_hr, ed_min, 0, 'format','yyyyMMdd HH:mm:ss'); %/ use eomday() to find the last day of the month
                
                if ~isempty(dt_slct_mo)       dates_dt = [date_st_dt:calmonths(dt_slct_mo):date_ed_dt]';
                elseif ~isempty(dt_slct_min)  dates_dt = [date_st_dt:minutes(dt_slct_min):date_ed_dt]';
                else                          dates_dt = [date_st_dt:hours(dt_slct_hr):date_ed_dt]';   end
                dates_dt(end) = [];
            
            elseif season == 4   %/ mth = 16; D(0)JF(1)
                if t == length(years_unique) && skip_the_incomplete
                    continue;  %/ then we skip the last year's DJF (incomplete)
                end
                date_st_dt = datetime(years_unique(t),   12,    1, st_hr, st_min, 0, 'format','yyyyMMdd HH:mm:ss');
                date_ed_dt = datetime(years_unique(t),   12, 31+1, ed_hr, ed_min, 0, 'format','yyyyMMdd HH:mm:ss'); %/ use eomday() to find the last day of the month
                if ~isempty(dt_slct_mo)       dates_dt1 = [date_st_dt:calmonths(dt_slct_mo):date_ed_dt]';
                elseif ~isempty(dt_slct_min)  dates_dt1 = [date_st_dt:minutes(dt_slct_min):date_ed_dt]';
                else                          dates_dt1 = [date_st_dt:hours(dt_slct_hr):date_ed_dt]';   end
                dates_dt1(end) = [];
                           
                date_st_dt = datetime(years_unique(t)+1,   1,   1,                              st_hr, st_min, 0, 'format','yyyyMMdd HH:mm:ss');
                date_ed_dt = datetime(years_unique(t)+1,   2,   eomday(years_unique(t)+1, 2)+1, ed_hr, ed_min, 0, 'format','yyyyMMdd HH:mm:ss'); %/ use eomday() to find the last day of the month
                if ~isempty(dt_slct_mo)       dates_dt2 = [date_st_dt:calmonths(dt_slct_mo):date_ed_dt]';
                elseif ~isempty(dt_slct_min)  dates_dt2 = [date_st_dt:minutes(dt_slct_min):date_ed_dt]';
                else                          dates_dt2 = [date_st_dt:hours(dt_slct_hr):date_ed_dt]';   end
                dates_dt2(end) = [];
                dates_dt = [dates_dt1; dates_dt2];
    
            elseif season == 6   %/ mth = 18; Cold season in Asia (Oct(0)-Mar(1))
                if t == length(years_unique) && skip_the_incomplete
                    continue;  %/ then we skip the last year's cold season (incomplete)
                end
                date_st_dt = datetime(years_unique(t),   10,     1, st_hr, st_min, 0, 'format','yyyyMMdd HH:mm:ss');
                date_ed_dt = datetime(years_unique(t),   12,  31+1, ed_hr, ed_min, 0, 'format','yyyyMMdd HH:mm:ss'); %/ use eomday() to find the last day of the month
                if ~isempty(dt_slct_mo)       dates_dt1 = [date_st_dt:calmonths(dt_slct_mo):date_ed_dt]';
                elseif ~isempty(dt_slct_min)  dates_dt1 = [date_st_dt:minutes(dt_slct_min):date_ed_dt]';
                else                          dates_dt1 = [date_st_dt:hours(dt_slct_hr):date_ed_dt]';   end
                dates_dt1(end) = [];
                
                date_st_dt = datetime(years_unique(t)+1,  1,     1, st_hr, st_min, 0, 'format','yyyyMMdd HH:mm:ss');
                date_ed_dt = datetime(years_unique(t)+1,  3,  31+1, ed_hr, ed_min, 0, 'format','yyyyMMdd HH:mm:ss'); %/ use eomday() to find the last day of the month
                if ~isempty(dt_slct_mo)       dates_dt2 = [date_st_dt:calmonths(dt_slct_mo):date_ed_dt]';
                elseif ~isempty(dt_slct_min)  dates_dt2 = [date_st_dt:minutes(dt_slct_min):date_ed_dt]';
                else                          dates_dt2 = [date_st_dt:hours(dt_slct_hr):date_ed_dt]';   end
                dates_dt2(end) = [];
                dates_dt = [dates_dt1; dates_dt2];
                
            elseif season == 7   %/ mth = 19; JFD(0)
                %/ NOTE: This season consider Jan, Feb and Dec in the same year,
                %/       useful when summing up seasonal values to annual value (s.t. no missing data in the final year)
                %/       Use with caution.
                
                date_st_dt = datetime(years_unique(t),   1,   1,                            st_hr, st_min, 0, 'format','yyyyMMdd HH:mm:ss');
                date_ed_dt = datetime(years_unique(t),   2,   eomday(years_unique(t), 2)+1, ed_hr, ed_min, 0, 'format','yyyyMMdd HH:mm:ss'); %/ use eomday() to find the last day of the month
                if ~isempty(dt_slct_mo)       dates_dt1 = [date_st_dt:calmonths(dt_slct_mo):date_ed_dt]';
                elseif ~isempty(dt_slct_min)  dates_dt1 = [date_st_dt:minutes(dt_slct_min):date_ed_dt]';
                else                          dates_dt1 = [date_st_dt:hours(dt_slct_hr):date_ed_dt]';   end
                dates_dt1(end) = [];
                
                date_st_dt = datetime(years_unique(t),   12,        1, st_hr, st_min, 0, 'format', 'yyyyMMdd HH:mm:ss');
                date_ed_dt = datetime(years_unique(t),   12, ed_day+1, ed_hr, ed_min, 0, 'format', 'yyyyMMdd HH:mm:ss'); %/ use eomday() to find the last day of the month
                if ~isempty(dt_slct_mo)       dates_dt2 = [date_st_dt:calmonths(dt_slct_mo):date_ed_dt]';
                elseif ~isempty(dt_slct_min)  dates_dt2 = [date_st_dt:minutes(dt_slct_min):date_ed_dt]';
                else                          dates_dt2 = [date_st_dt:hours(dt_slct_hr):date_ed_dt]';   end
                dates_dt2(end) = [];
                dates_dt = [dates_dt1; dates_dt2];
                
            elseif season == 8   %/ mth = 20; nonJJA (year 0)
                %/ NOTE: This season consider Jan, Feb and Dec in the same year,
                %/       useful when summing up seasonal values to annual value (s.t. no missing data in the final year)
                %/       Use with caution.
                
                date_st_dt = datetime(years_unique(t),   1,      1, st_hr, st_min, 0, 'format','yyyyMMdd HH:mm:ss');
                date_ed_dt = datetime(years_unique(t),   5,   31+1, ed_hr, ed_min, 0, 'format','yyyyMMdd HH:mm:ss'); %/ use eomday() to find the last day of the month
                if ~isempty(dt_slct_mo)       dates_dt1 = [date_st_dt:calmonths(dt_slct_mo):date_ed_dt]';
                elseif ~isempty(dt_slct_min)  dates_dt1 = [date_st_dt:minutes(dt_slct_min):date_ed_dt]';
                else                          dates_dt1 = [date_st_dt:hours(dt_slct_hr):date_ed_dt]';   end
                dates_dt1(end) = [];
                
                date_st_dt = datetime(years_unique(t),    9,        1, st_hr, st_min, 0, 'format', 'yyyyMMdd HH:mm:ss');
                date_ed_dt = datetime(years_unique(t),   12, ed_day+1, ed_hr, ed_min, 0, 'format', 'yyyyMMdd HH:mm:ss'); %/ use eomday() to find the last day of the month
                if ~isempty(dt_slct_mo)       dates_dt2 = [date_st_dt:calmonths(dt_slct_mo):date_ed_dt]';
                elseif ~isempty(dt_slct_min)  dates_dt2 = [date_st_dt:minutes(dt_slct_min):date_ed_dt]';
                else                          dates_dt2 = [date_st_dt:hours(dt_slct_hr):date_ed_dt]';   end
                dates_dt2(end) = [];
                dates_dt = [dates_dt1; dates_dt2];
                
            elseif season == 10   %/ mth = 22; NDJFMA (Nov(0)-Apr(1))
                if t == length(years_unique) && skip_the_incomplete
                    continue;  %/ then we skip the last year's cold season (incomplete)
                end
                date_st_dt = datetime(years_unique(t),   11,     1, st_hr, st_min, 0, 'format','yyyyMMdd HH:mm:ss');
                date_ed_dt = datetime(years_unique(t),   12,  31+1, ed_hr, ed_min, 0, 'format','yyyyMMdd HH:mm:ss'); %/ use eomday() to find the last day of the month
                if ~isempty(dt_slct_mo)       dates_dt1 = [date_st_dt:calmonths(dt_slct_mo):date_ed_dt]';
                elseif ~isempty(dt_slct_min)  dates_dt1 = [date_st_dt:minutes(dt_slct_min):date_ed_dt]';
                else                          dates_dt1 = [date_st_dt:hours(dt_slct_hr):date_ed_dt]';   end
                dates_dt1(end) = [];
                
                date_st_dt = datetime(years_unique(t)+1,  1,     1, st_hr, st_min, 0, 'format','yyyyMMdd HH:mm:ss');
                date_ed_dt = datetime(years_unique(t)+1,  4,  30+1, ed_hr, ed_min, 0, 'format','yyyyMMdd HH:mm:ss'); %/ use eomday() to find the last day of the month
                if ~isempty(dt_slct_mo)       dates_dt2 = [date_st_dt:calmonths(dt_slct_mo):date_ed_dt]';
                elseif ~isempty(dt_slct_min)  dates_dt2 = [date_st_dt:minutes(dt_slct_min):date_ed_dt]';
                else                          dates_dt2 = [date_st_dt:hours(dt_slct_hr):date_ed_dt]';   end
                dates_dt2(end) = [];
                dates_dt = [dates_dt1; dates_dt2];

            else
                error('Invalid input of ''season''!')
            end
        else
            if isempty(ed_day) 
                if ~isempty(ed_month) && ismember(ed_month, [1:12])  
                    ed_day = eomday(years_unique(t),ed_month);  
                end
            end
        
            date_st_dt = datetime(years_unique(t), st_month,  st_day,   st_hr, st_min, 0, 'format','yyyyMMdd HH:mm:ss');
            date_ed_dt = datetime(years_unique(t), ed_month,  ed_day+1, ed_hr, ed_min, 0, 'format','yyyyMMdd HH:mm:ss'); %/ use eomday() to find the last day of the month
            %/ Remove the last date (as added temporarily) to get the full time
            %/ steps of the period.
            if ~isempty(dt_slct_mo)       dates_dt = [date_st_dt:calmonths(dt_slct_mo):date_ed_dt]';
            elseif ~isempty(dt_slct_min)  dates_dt = [date_st_dt:minutes(dt_slct_min):date_ed_dt]';
            else                          dates_dt = [date_st_dt:hours(dt_slct_hr):date_ed_dt]';   end
            dates_dt(end) = [];
        end
        
        if isequal(output_date_format, 'yyyymmdd') || isequal(output_date_format, 'yyyyMMdd')
            output_dates = cat(1, output_dates, yyyymmdd(dates_dt));

        elseif isequal(output_date_format, 'yyyymmddHHMM') || isequal(output_date_format, 'yyyyMMddHHmm')
            output_dates = cat(1, output_dates, dates_dt.Year*1e8 + dates_dt.Month*1e6 + dates_dt.Day*1e4 + dates_dt.Hour*1e2 + dates_dt.Minute);

        elseif isequal(output_date_format, 'yyyymm') || isequal(output_date_format, 'yyyyMM')
            output_dates = cat(1, output_dates, unique(floor(yyyymmdd(dates_dt)./100)));

        elseif isequal(output_date_format, 'yyyyptd')
            %/ recursion
            dates_mmdd          = mod(date_array_gen('year_list', 1979, 'output_date_format', 'yyyymmdd'), 1e4); %/ 365
            dates_ptd_centerday = dates_mmdd(3:5:end);

            %/ Check if the dates cover the start day of any pentads
            ptd = findismember_loop(dates_ptd_centerday, mod(yyyymmdd(dates_dt), 1e4));

            % dates_dt(ptd).Year
            output_dates = cat(1, output_dates, dates_dt(ptd).Year*1e2 + ptd);

        elseif isequal(output_date_format, 'datetime')
            output_dates = cat(1, output_dates, dates_dt); %/ keep datetime as it is  
        else
            error('Invalid ''output_date_format''!');
        end

        %/ Convert double into int64
        if ~isequal(output_date_format, 'datetime')
            output_dates = int64(output_dates);  
        end
        
        if flag_break
            break;
        end
    end

    if noleap
        if isequal(calendar, 'julian')
            error('calendar julian conflicts with noleap = 1!');
        end
        if isequal(output_date_format, 'yyyymmdd') || isequal(output_date_format, 'yyyyMMdd')
            ind_leap = find(mod(output_dates, 1e4) == 229);
    
        elseif isequal(output_date_format, 'yyyymmddHHMM') || isequal(output_date_format, 'yyyyMMddHHmm')
            ind_leap = find(floor(mod(output_dates, 1e8)./1e4) == 229);

        elseif isequal(output_date_format, 'yyyymm') || isequal(output_date_format, 'yyyyMM')
            ind_leap = [];

        elseif isequal(output_date_format, 'yyyyptd') 
            ind_leap = [];

        elseif isequal(output_date_format, 'datetime')
            ind_leap = find(output_dates.Month == 2 & output_dates.Day == 29);
        else
            error('Invalid ''output_date_format''!');
        end
        output_dates(ind_leap) = [];  
    end
end
