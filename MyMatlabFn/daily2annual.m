%%
function [annual_data, annual_SD, annual_date] = daily2annual(varargin)

    %/ create a set of valid parameters and their default value
    pnames = {'daily_data', 'dates',  'ins_or_acc'}; 
    dflts  = {          [],      [],             1};
    %/ parse function arguments
    [           daily_data,   dates,   ins_or_acc] = internal.stats.parseArgs(pnames, dflts, varargin{:});

    %%
    if numel(num2str(dates(1))) ~= 8   error('Check if input date is in yyyymmdd!!');  end
    if ~ismember(ins_or_acc, [1, 2])   error('''mean_or_sum'' can only be 1 or 2!');   end
    
    sz = size(daily_data);
    dim_time = length(sz);
    
    %/ If data is a 1D array, transpose it a row vector s.t. time is at the 2nd dim
    if dim_time == 2 && sz(2) == 1
        daily_data = daily_data';
    end
    
    year_list = unique(floor(dates/1e4));
    if dim_time == 2 
        annual_data = nan(sz(1), length(year_list));
        for t = 1:length(year_list)
            ind = find(floor(dates/1e4) == year_list(t));
            
            if ins_or_acc == 1
                annual_data(:,t) = nanmean(daily_data(:,ind), dim_time);
            elseif ins_or_acc == 2
                annual_data(:,t) = nansum(daily_data(:,ind), dim_time);
            end
        end          
        
    elseif dim_time == 3
        annual_data = nan(sz(1), sz(2), length(year_list));
        for t = 1:length(year_list)
            ind = find(floor(dates/1e4) == year_list(t));
            
            if ins_or_acc == 1
                annual_data(:,:,t) = nanmean(daily_data(:,:,ind),dim_time);
            elseif ins_or_acc == 2
                annual_data(:,:,t) = nansum(daily_data(:,:,ind),dim_time);
            end
        end          
        
    else
        error('Code not set for dim_time == %d!', dim_time)
    end
        
    annual_SD = std(annual_data,0,dim_time);
    annual_date = year_list;
end