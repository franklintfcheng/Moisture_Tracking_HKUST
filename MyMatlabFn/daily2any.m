function [output_data, output_SD, output_dates, output_data_prime] = daily2any(varargin)

    %/ create a set of valid parameters and their default value
    pnames = {'data_daily', 'dates',  'mth', 'ins_or_acc'}; 
    dflts  = cell(length(pnames), 1);
    %/ parse function arguments
    [          data_daily,   dates,    mth,   ins_or_acc] = internal.stats.parseArgs(pnames, dflts, varargin{:});

    %%
    %======================================================================
    %/ Author: Franklin Cheng (fandycheng@ust.hk)
    %/ Last update: 9 Mar 2024
    %/
    %/ Description :
    %/      This function converts daily data to monthly/seasonal/yearly data
    %/
    %/ NOTE:
    %/      It works for *monthly* data as well if the input dates are in 
    %/      the format of 'yyymmdd' or 'yyyymm'
    %======================================================================
    L = numel(num2str(dates(1)));

    if L == 8      %/ yyyymmdd
        CONV = 1./1e4;
        output_date_format = 'yyyymmdd';
    elseif L == 6  %/ yyyymm
        CONV = 1./1e2;
        output_date_format = 'yyyymm';
    else
        error('Invalid datetime format with digit = %d. Please check.', L);
    end

    if isempty(ins_or_acc)             
        ins_or_acc = 'ins';
    elseif ~ismember(ins_or_acc, {'ins', 'acc'})   
        error('''ins_or_acc'' can only be ''ins'' or ''acc''!');
    end
    if isempty(data_daily)   
        error('data_daily is empty!');  
    end
    
    sz       = size(data_daily);
    dim_time = length(sz);
    
%     output_data_prime = nan(sz);  %/ compute prime (from Reynolds average) only for 'ins'
    if isequal(ins_or_acc, 'ins')
        output_data_prime = nan(sz);  %/ compute prime (from Reynolds average) only for 'ins'
    else
        % warning('[daily2any]: No prime (anomaly) is computed for ''acc''! Since by def., it is the departure from the mean, not from the sum!')
        output_data_prime = [];
    end
    
    %/ If data is a 1D array, transpose it a row vector s.t. time is at the 2nd dim
    if dim_time == 2 && sz(2) == 1
        data_daily = data_daily';
    end
    
    %/ Get the date indices
    year_list    = unique(floor(dates*CONV));
    new_sz       = [sz(1:end-1), length(year_list)];
    output_data  = nan(new_sz);
%     output_dates = cell(length(years),1);
    flag_set_nan = 0;
    for y = 1:length(year_list)
        if mth == 0
            ind_dates    = find(floor(dates*CONV) == year_list(y));               %/ yearly 
            % subset_dates = dates(ind_dates);
            
        elseif mth >= 1 && mth <= 12
            subset_dates = date_array_gen('year_list', year_list(y), 'mth', mth, 'output_date_format', output_date_format);
            ind_dates    = findismember_loop(dates, subset_dates);           %/ seasonal

        elseif mth >= 13 
            season       = mth - 12;
            subset_dates = date_array_gen('year_list', year_list(y), 'season', season, 'output_date_format', output_date_format);
            ind_dates    = findismember_loop(dates, subset_dates);           %/ seasonal
            
            if length(ind_dates) ~= length(subset_dates)
                if y == length(year_list) && ismember(mth, [16, 18])
                    flag_set_nan = 1;
                    warning('[daily2any]: Note: As mth == %d, the final year of the yearly_data is set NaN', mth);
                else
                    %/ Known Issues: no Jan 1 for all CERA-20C P and E data.
                    warning('[daily2any]: Detected Incomplete date coverage for year = %d, mth = %d!', year_list(y), mth)
                end
            end
        else
            error('Invalid input of mth!')
        end
        
        %/ Consider various dim of data
        if dim_time == 2
            if flag_set_nan
                output_data(:,y) = nan;
            else
                %/ If the field is ins -> take the mean
                if isequal(ins_or_acc, 'ins')
                    output_data(:,y) = mean(data_daily(:,ind_dates), dim_time, 'omitnan');
                    output_data_prime(:,ind_dates) = data_daily(:,ind_dates) - output_data(:,y); %/ departure from the monthly/seasonal/annual mean
                elseif isequal(ins_or_acc, 'acc')
                    output_data(:,y) = sum(data_daily(:,ind_dates), dim_time, 'omitnan');
                end
            end
        elseif dim_time == 3
            if flag_set_nan
                output_data(:,y) = nan;
            else
                %/ If the field is ins -> take the mean
                if isequal(ins_or_acc, 'ins')
                    output_data(:,:,y) = mean(data_daily(:,:,ind_dates), dim_time, 'omitnan');
                    output_data_prime(:,:,ind_dates) = data_daily(:,:,ind_dates) - output_data(:,:,y); %/ departure from the monthly/seasonal/annual mean
                elseif isequal(ins_or_acc, 'acc')
                    output_data(:,:,y) = sum(data_daily(:,:,ind_dates), dim_time, 'omitnan');
                end
            end
        elseif dim_time == 4
            if flag_set_nan
                output_data(:,y) = nan;
            else
                %/ If the field is ins -> take the mean
                if isequal(ins_or_acc, 'ins')
                    output_data(:,:,:,y) = mean(data_daily(:,:,:,ind_dates), dim_time, 'omitnan');
                    output_data_prime(:,:,:,ind_dates) = data_daily(:,:,:,ind_dates) - output_data(:,:,:,y); %/ prime = departure from the monthly/seasonal/annual mean
                elseif isequal(ins_or_acc, 'acc')
                    output_data(:,:,:,y) = sum(data_daily(:,:,:,ind_dates), dim_time, 'omitnan');
                end
            end
        else
            error('Code not set for the dim_time == %d!', dim_time);
        end
    end
%     size(output_data)
    output_SD    = std(output_data,0,dim_time, 'omitnan');
    output_dates = year_list;
%     output_dates = cat(1, output_dates{:});
        
 
end