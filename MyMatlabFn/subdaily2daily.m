function [daily_data, date_yyyymmdd_AllYr]  = subdaily2daily(varargin)

    pnames = {'subdaily_data', 'dates', 'ins_or_acc'};
    dflts  =  cell(1, length(pnames));
    [          subdaily_data,   dates,   ins_or_acc] ...
                            = internal.stats.parseArgs(pnames, dflts, varargin{:});

%% 
    %/=====================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last Update: 26 Jan 2023
    %/
    %/ Description: This function converts subdaily data to daily.
    %/              Assuming dim_time to be the last dim.
    %/             
    %/=====================================================================
    Lr = 12;  %/ yyyymmddHHMM
    if numel(num2str(dates(1))) ~= Lr                               
        error('Check if ''dates'' is in yyyymmddHHMM (%d digits)!!', Lr);              
    end
    if ~isequal(ins_or_acc, 'ins') && ~isequal(ins_or_acc, 'acc')   
        error('''ins_or_acc'' can only be ''ins'' or ''acc''!');    
    end
    
    dim_time = length(size(subdaily_data));
    if dim_time == 2 && size(subdaily_data, 1) ~= 1  %/ convert a column to a row vector -> make the time_dim the last dim.
        subdaily_data  = subdaily_data';  
        flag_transpose = 1;
    else
        flag_transpose = 0;
    end

    sz = size(subdaily_data);            %/ put this line *after* transpose (if so).
    year_list = unique(floor(dates/10^(Lr-4)));

    %/ For double-checking
    date_yyyymmdd_AllYr = date_array_gen('year_list', year_list, 'dt_slct_hr', 24, 'output_date_format', 'yyyymmdd');          
    nt = length(date_yyyymmdd_AllYr);

    %/ Loop to take mean/sum
    daily_data = nan([sz(1:end-1), nt]); %/ Initialize data
    for t = 1:nt
        ind = find(floor(dates/10^(Lr-8)) == date_yyyymmdd_AllYr(t));
        if isempty(ind)
            warning('Missing data for date %d! Will stay as NaN.', date_yyyymmdd_AllYr(t));
        else
            if isequal(ins_or_acc, 'ins')
                if dim_time == 2
                    daily_data(:,t) = mean(subdaily_data(:,ind), dim_time, 'omitnan');
                elseif dim_time == 3
                    daily_data(:,:,t) = mean(subdaily_data(:,:,ind), dim_time, 'omitnan');
                elseif dim_time == 4
                    daily_data(:,:,:,t) = mean(subdaily_data(:,:,:,ind), dim_time, 'omitnan');
                else
                    error('Code not ready for %dD data!', dim_time);
                end

            elseif isequal(ins_or_acc, 'acc')
                if dim_time == 2
                    daily_data(:,t) = sum(subdaily_data(:,ind), dim_time, 'omitnan');
                elseif dim_time == 3
                    daily_data(:,:,t) = sum(subdaily_data(:,:,ind), dim_time, 'omitnan');
                elseif dim_time == 4
                    daily_data(:,:,:,t) = sum(subdaily_data(:,:,:,ind), dim_time, 'omitnan');
                else
                    error('Code not ready for %dD data!', dim_time);
                end
            end
        end
    end
    
    if ~isempty(find(isnan(daily_data), 1))
        warning('[subdaily2daily] The converted daily_data contains Nan (likely due to missing values!)');
    end

    if flag_transpose
        daily_data = daily_data'; %/ transpose back to the original (only for a column vector)
    end
    
end