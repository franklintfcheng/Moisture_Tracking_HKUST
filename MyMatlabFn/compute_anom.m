function [X, X_mean, X_mean_sig, date] = compute_anom(varargin)
    
    %/ create a set of valid parameters and their default value
    pnames = {'dataset', 'var', 'select_field', 'mth', 'st_month', 'st_day', 'ed_month', 'ed_day', 'year_list', 'compute_anombyAM', 'alpha'};  
    dflts  = cell(1, length(pnames));

    [           dataset,   var,   select_field,   mth,   st_month,  st_day,   ed_month,   ed_day,   year_list,   compute_anombyAM,   alpha] ...
                   = internal.stats.parseArgs(pnames, dflts, varargin{:}); %/ parse function arguments

%%
    %/ If the dim of date is only 1
    dim_time = 3;  %/ Assume
    if size(dataset.(var).(select_field),dim_time) == 1
        X = dataset.(var).(select_field);
        X_mean = dataset.(var).(select_field);
        X_mean_sig = X_mean;
        date = [];
        return;
    end
    
    if compute_anombyAM   
        mth_list = [mth, 0];   
    else  
        mth_list = mth;   
    end

    for m = 1:length(mth_list)
        mth_bc = mth_list(m);
        fprintf('*** mth_bc = %d ***\n', mth_bc);
        
        season = []; skip_the_incomplete = 0;
        if isempty(st_month) && isempty(st_day) && isempty(ed_month) && isempty(ed_day)
            if mth_bc == 0
                st_month = 1; st_day = 1; ed_month = 12; ed_day = 31; 
            elseif mth_bc >= 13       
                %/ 13: MAM, 14: JJA, 15: SON, 16: DJF, 17: Apr-Sep, 18: Oct-Mar
                season = mth_bc - 12;     
                % skip_the_incomplete = 1; %/ necessary?
            else
                st_month  = mth_bc; st_day = 1;
                ed_month  = st_month; 
                ed_day    = eomday(year_bc,ed_month);
            end
        end

        %/ slct_dates_meteo will always be yyyymmdd
        output_date_format = 'yyyymmdd';
        slct_dates_meteo = date_array_gen('year_list', year_list, 'st_month', st_month, 'st_day', st_day, 'ed_month', ed_month, 'ed_day', ed_day,...
                                          'season', season, 'output_date_format', output_date_format, 'skip_the_incomplete', skip_the_incomplete);
        
        if contains(select_field, 'monthly') && isfield(dataset.(var), 'date_yyyymm')
            date_fld = 'date_yyyymm';
            CONV     = 1/100;
        else
            date_fld = 'date_yyyymmdd_AllYr';
            CONV     = 1;
        end

        %/ In case the older dataset does not contain the 'date_yyyymm' field for monthly data
        if isfield(dataset.(var), date_fld)
            ind_date = find(ismember(dataset.(var).(date_fld), slct_dates_meteo*CONV));
        else
            date_fld = 'date_yyyymmdd_AllYr';
            ind_date = find(ismember(dataset.(var).(date_fld)*CONV, slct_dates_meteo*CONV));
        end

        if compute_anombyAM && mth_bc == 0   %/ compute annual mean (AM)
            X_AM = mean(dataset.(var).(select_field)(:,:,ind_date), 3, 'omitnan');
            X = X - X_AM;
        else                                 %/ just compute seasonal / monthly
            X = dataset.(var).(select_field)(:,:,ind_date);
            date = dataset.(var).(date_fld)(ind_date);
        end
        X_mean = mean(X, 3, 'omitnan');
        
        %/ anomaly mean (do t-test if alpha not empty)
        if ~isempty(alpha)
            X_mean_sig = ttest_sig_fn(X, alpha, 3);
        else
            X_mean_sig = [];
        end
    end
end