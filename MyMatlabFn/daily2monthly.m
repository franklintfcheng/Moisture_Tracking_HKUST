function [monthly_data, date_yyyymm] = daily2monthly(varargin)

    %/ create a set of valid parameters and their default value
    pnames = {'daily_data', 'dates', 'ins_or_acc', 'n_MA',  'noleap'}; 
    dflts  = {          [],      [],        'inc',     [],        0};
    [           daily_data,   dates,   ins_or_acc,   n_MA,   noleap] = internal.stats.parseArgs(pnames, dflts, varargin{:});
%%

    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 18 Jun 2024
    %/
    %/ Description: This function converts daily data into monthly.
    %/              Assuming time_dim to be the last dim.
    %/=====================================================================
    
    if numel(num2str(dates(1))) ~= 8                                
        error('Check if ''dates'' is in yyyymmdd (8 digits)!!');              
    end
    if ~isequal(ins_or_acc, 'ins') && ~isequal(ins_or_acc, 'acc')   
        error('''ins_or_acc'' can only be ''ins'' or ''acc''!');    
    end
    
    time_dim = length(size(daily_data));
    if time_dim == 2 && size(daily_data, 1) ~= 1  %/ convert a column to a row vector -> make the time_dim the last dim.
        daily_data     = daily_data';  
        flag_transpose = 1;
    else
        flag_transpose = 0;
    end

    sz = size(daily_data);
    if ~isequal(sz(end), length(dates))
        error('Inconsistent nubmer of dates with daily_data. Asumming the last dimension of daily_data is dates');
    end

    sz          = size(daily_data);            %/ put this line *after* transpose (if so).
    date_yyyymm = unique(floor(dates/1e2));

    %/ Loop to take mean/sum
    monthly_data = nan([sz(1:end-1), length(date_yyyymm)]); %/ Initialize data
    for t = 1:length(date_yyyymm)

        y = floor(date_yyyymm(t)/1e2);
        m = mod(date_yyyymm(t), 1e2);

        ind = find(floor(dates/1e2) == date_yyyymm(t));
        if isempty(ind)
            error('empty ind!');
        end
        
        %/ Check if all dates are available to compute the monthly mean/sum (allow for 'noleap')
        if noleap == 1 && leapyear(y) && m == 2
            DoM = eomday(y,m)-1;
        else
            DoM = eomday(y,m);
        end
        if numel(ind) ~= DoM
            warning('Only %d days are found in the month of %d (noleap == %d)! No calculation will be done.', numel(ind), date_yyyymm(t), noleap)
            continue;
        end

        if isequal(ins_or_acc, 'ins')
            if time_dim == 2
                monthly_data(:,t) = mean(daily_data(:,ind), time_dim, 'omitnan');
            elseif time_dim == 3
                monthly_data(:,:,t) = mean(daily_data(:,:,ind), time_dim, 'omitnan');
            elseif time_dim == 4
                monthly_data(:,:,:,t) = mean(daily_data(:,:,:,ind), time_dim, 'omitnan');
            else
                error('Code not ready for %dD data!', time_dim);
            end
            
        elseif isequal(ins_or_acc, 'acc')
            if time_dim == 2
                monthly_data(:,t) = sum(daily_data(:,ind), time_dim, 'omitnan');
            elseif time_dim == 3
                monthly_data(:,:,t) = sum(daily_data(:,:,ind), time_dim, 'omitnan');
            elseif time_dim == 4
                monthly_data(:,:,:,t) = sum(daily_data(:,:,:,ind), time_dim, 'omitnan');
            else
                error('Code not ready for %dD data!', time_dim);
            end
        end
    end

    % %/ Loop to take mean/sum
    % year_list           = unique(floor(dates/1e4));
    % date_yyyymmdd_AllYr = nan(length(year_list)*12, 1); %/ for dates, make it a column vector.
    % monthly_data = nan([sz(1:end-1), length(year_list)*12]); %/ Initialize data
    % for y = 1:length(year_list)
    %     for m = 1:12
    %         ind = find(floor(dates/1e2) == year_list(y)*1e2 + m);
    %         if isempty(ind)
    %             error('empty ind!');
    %         end
    %         if isequal(ins_or_acc, 'ins')
    %             if time_dim == 2
    %                 monthly_data(:,m+(y-1)*12) = mean(daily_data(:,ind), time_dim, 'omitnan');
    %             elseif time_dim == 3
    %                 monthly_data(:,:,m+(y-1)*12) = mean(daily_data(:,:,ind), time_dim, 'omitnan');
    %             elseif time_dim == 4
    %                 monthly_data(:,:,:,m+(y-1)*12) = mean(daily_data(:,:,:,ind), time_dim, 'omitnan');
    %             else
    %                 error('Code not ready for %dD data!', time_dim);
    %             end
    % 
    %         elseif isequal(ins_or_acc, 'acc')
    %             if time_dim == 2
    %                 monthly_data(:,m+(y-1)*12) = sum(daily_data(:,ind), time_dim, 'omitnan');
    %             elseif time_dim == 3
    %                 monthly_data(:,:,m+(y-1)*12) = sum(daily_data(:,:,ind), time_dim, 'omitnan');
    %             elseif time_dim == 4
    %                 monthly_data(:,:,:,m+(y-1)*12) = sum(daily_data(:,:,:,ind), time_dim, 'omitnan');
    %             else
    %                 error('Code not ready for %dD data!', time_dim);
    %             end
    %         end
    %         date_yyyymmdd_AllYr(m+(y-1)*12) = dates(ind(1));
    %     end   
    % end
    
    if ~isempty(n_MA)  %/ if requesting for a n-month moving average.
        monthly_data = movmean(monthly_data, n_MA, time_dim); 
    end
    
    if flag_transpose
        monthly_data = monthly_data'; %/ transpose back to the original (only for a column vector)
    end
    
end